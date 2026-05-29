#!/usr/bin/env python3
"""
Gene-based (collapsing/burden) association test from a VCF, in pure Python.

    raw per-chromosome VCF (bgzipped + tabix-indexed)
      -> assign each variant to overlapping gene(s)        [pyensembl, GRCh37]
      -> per-sample rare-variant burden score per gene      [cyvcf2 streaming]
      -> phenotype ~ burden + covariates                    [statsmodels OLS/Logit]
      -> per-gene table: beta (with sign), se, p, BH-FDR q

This is the FRONT HALF of the POC. The back half hands the significant genes to
INDRA's discovery engine; the output table's `gene_name` / `gene_id` columns are
what you ground to HGNC for that step. Burden (not SKAT) is used on purpose: it
yields an effect *direction*, which is what INDRA's signed Reverse Causal
Reasoning consumes.

--- CORRECTNESS NOTES (read these) -------------------------------------------
* GENOME BUILD. 1000G Phase 3 (20130502, v5a) is GRCh37/hg19. Gene coordinates
  must come from a GRCh37 Ensembl release -> release 75. Mixing builds silently
  corrupts every variant->gene assignment. This is the #1 way to get nonsense.
* DOSAGE ENCODING. We open the VCF with gts012=True, so variant.gt_types returns
  0/1/2 = alt-allele count and 3 = missing. Verify on your cyvcf2 version by
  printing a few variants and checking the allele frequencies look sane.
* INPUT FILE. A real rare-variant burden test needs the UNFILTERED callset.
  A file named '...bi_maf001...' has already stripped the rare tail this test is
  built to aggregate; use the full ALL.chrN....genotypes.vcf.gz from IGSR/EBI.

Setup:
    pip install cyvcf2 pyensembl statsmodels numpy pandas intervaltree
    pyensembl install --release 75 --species homo_sapiens
"""

from __future__ import annotations

import argparse
import logging
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from cyvcf2 import VCF
from pyensembl import EnsemblRelease
from intervaltree import IntervalTree

log = logging.getLogger("burden")


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def load_phenotype(path, sample_col, pheno_col, covar_cols):
    """Read a CSV/TSV of phenotype + covariates, indexed by sample id.

    Expected columns: <sample_col>, <pheno_col>, and any covariates (e.g. PC1..PCk,
    sex). 1000G has real structure across superpopulations, so include ancestry
    PCs as covariates here -- compute them once (e.g. scikit-allel PCA on a
    pruned genome-wide set) and add them as columns.
    """
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    keep = [sample_col, pheno_col] + list(covar_cols)
    missing = [c for c in keep if c not in df.columns]
    if missing:
        sys.exit(f"phenotype file is missing columns: {missing}")
    df = df[keep].rename(columns={sample_col: "sample.id", pheno_col: "pheno"})
    return df.set_index("sample.id")


def build_gene_tree(release, contig, biotype="protein_coding"):
    """Interval tree of gene intervals on one contig -> {gene_id: (name, tree)}.

    A variant overlapping several genes is counted toward each, matching the
    'a variant may belong to more than one group' convention.
    """
    data = EnsemblRelease(release)
    tree = IntervalTree()
    names = {}
    n = 0
    for g in data.genes():
        if g.contig != contig or (biotype and g.biotype != biotype):
            continue
        if g.end < g.start:
            continue
        tree[g.start : g.end + 1] = g.gene_id  # end-exclusive
        names[g.gene_id] = g.gene_name or g.gene_id
        n += 1
    if n == 0:
        sys.exit(
            f"no {biotype} genes found on contig '{contig}' in release {release}. "
            "Check the contig name (1000G uses '17', not 'chr17') and that "
            "`pyensembl install --release {release}` has run."
        )
    log.info("loaded %d %s genes on contig %s (release %d)", n, biotype, contig, release)
    return names, tree


# --------------------------------------------------------------------------- #
# Genotype -> per-gene burden
# --------------------------------------------------------------------------- #
def alt_dosage(variant):
    """Alt-allele dosage per sample as float; missing -> nan. Assumes gts012=True."""
    gt = variant.gt_types.astype(float)  # 0,1,2 = alt count; 3 = missing
    gt[gt == 3] = np.nan
    return gt


def accumulate_burden(vcf_path, contig, region, tree, samples, maf, weights):
    """Stream the VCF and build {gene_id: burden_vector} and per-gene variant counts.

    `samples` is the ordered list of sample ids to analyze; the VCF is opened
    restricted to them so dosage vectors line up with the regression design.
    Only rare (0 < alt_AF <= maf), biallelic variants are aggregated. Missing
    genotypes are imputed to 0 (no minor alleles); switch to mean-imputation if
    you prefer. `weights`='count' gives a simple collapsing burden; 'beta' applies
    Beta(1,25) up-weighting of rarer variants (needs scipy).
    """
    if weights == "beta":
        from scipy.stats import beta as _beta  # lazy import

    n_samples = len(samples)
    burden, nvar = {}, {}
    vcf = VCF(vcf_path, gts012=True, samples=samples)  # see DOSAGE ENCODING note
    iterator = vcf(region) if region else vcf
    kept = scanned = 0
    for v in iterator:
        scanned += 1
        if v.CHROM != contig or len(v.ALT) != 1:
            continue
        dose = alt_dosage(v)
        nonmiss = ~np.isnan(dose)
        if nonmiss.sum() == 0:
            continue
        af = np.nansum(dose) / (2.0 * nonmiss.sum())
        if not (0.0 < af <= maf):  # keep rare alt alleles only
            continue
        overlaps = tree[v.POS]
        if not overlaps:
            continue
        w = float(_beta.pdf(af, 1, 25)) if weights == "beta" else 1.0
        contrib = np.nan_to_num(dose, nan=0.0) * w
        for iv in overlaps:
            gid = iv.data
            if gid not in burden:
                burden[gid] = np.zeros(n_samples)
                nvar[gid] = 0
            burden[gid] += contrib
            nvar[gid] += 1
        kept += 1
    log.info("scanned %d variants, aggregated %d rare biallelic into %d genes",
             scanned, kept, len(burden))
    return burden, nvar


# --------------------------------------------------------------------------- #
# Per-gene regression
# --------------------------------------------------------------------------- #
def test_genes(burden, nvar, y, covars, names, trait_type, min_variants):
    """Regress phenotype on each gene's burden score + covariates."""
    binary = trait_type == "binary" or (
        trait_type == "auto" and np.unique(y[~np.isnan(y)]).size <= 2
    )
    log.info("trait treated as %s", "binary (Logit)" if binary else "quantitative (OLS)")

    base = np.column_stack([np.ones(len(y)), covars]) if covars.size else np.ones((len(y), 1))
    rows = []
    for gid, score in burden.items():
        if nvar[gid] < min_variants or np.std(score) == 0:
            continue
        X = np.column_stack([base, score])
        try:
            if binary:
                res = sm.Logit(y, X).fit(disp=0)
            else:
                res = sm.OLS(y, X).fit()
            beta = res.params[-1]
            rows.append({
                "gene_id": gid,
                "gene_name": names.get(gid, gid),
                "n_variants": nvar[gid],
                "n_carriers": int((score > 0).sum()),
                "beta": beta,
                "direction": "+" if beta > 0 else "-",
                "se": res.bse[-1],
                "pval": res.pvalues[-1],
            })
        except Exception as e:  # separation / non-convergence on rare genes
            log.debug("skip %s: %s", gid, e)

    out = pd.DataFrame(rows)
    if not out.empty:
        ok = out["pval"].notna()
        out.loc[ok, "qval"] = multipletests(out.loc[ok, "pval"], method="fdr_bh")[1]
        out = out.sort_values("pval", na_position="last").reset_index(drop=True)
    return out


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", help="bgzipped + tabix-indexed VCF (one chromosome)", default='data/GWAS_1kg-genotypes_vcf_ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.bgz')
    ap.add_argument("--pheno",  help="CSV/TSV of phenotype + covariates", default='data/phenotypes.csv')
    ap.add_argument("--contig", default="17", help="contig name as it appears in the VCF (e.g. '17')")
    ap.add_argument("--region", default=None, help="optional region for cyvcf2, e.g. '17' or '17:1-83000000'")
    ap.add_argument("--release", type=int, default=75, help="Ensembl release (75 = GRCh37)")
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--pheno-col", default="ldl")
    ap.add_argument("--covars", default='annotated_sex,age_at_index,population_ASW,population_BEB,population_CDX,population_CEU,population_CHB,population_CHS,population_CLM,population_ESN,population_FIN,population_GBR,population_GIH,population_GWD,population_IBS,population_ITU,population_JPT,population_KHV,population_LWK,population_MSL,population_MXL,population_PEL,population_PJL,population_PUR,population_STU,population_TSI,population_YRI', help="comma-separated covariate column names (e.g. PC1,PC2,PC3,sex)")
    ap.add_argument("--trait-type", choices=["auto", "quantitative", "binary"], default="auto")
    ap.add_argument("--maf", type=float, default=0.01, help="max alt allele frequency to include")
    ap.add_argument("--weights", choices=["count", "beta"], default="count")
    ap.add_argument("--min-variants", type=int, default=2, help="min rare variants for a gene to be tested")
    ap.add_argument("--out", default="data/gene_burden_results.tsv")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    covar_cols = [c for c in args.covars.split(",") if c]
    pheno = load_phenotype(args.pheno, args.sample_col, args.pheno_col, covar_cols)

    # Align phenotype to VCF sample order, then drop incomplete cases.
    vcf_samples = VCF(args.vcf).samples
    shared = [s for s in vcf_samples if s in pheno.index]
    if not shared:
        sys.exit("no overlap between VCF samples and phenotype sample ids")
    pheno = pheno.loc[shared]
    complete = pheno[["pheno"] + covar_cols].notna().all(axis=1)
    pheno = pheno[complete]
    sample_order = pheno.index.tolist()
    log.info("%d samples after intersect + complete-case", len(sample_order))

    names, tree = build_gene_tree(args.release, args.contig)

    # cyvcf2 emits dosage in *its* sample order, which may not match sample_order,
    # so capture that order and remap dosage columns back to sample_order. Note:
    # accumulate_burden opens its own VCF handle; we pass the matching order here.
    vcf_order = VCF(args.vcf, samples=sample_order).samples
    pos = {s: i for i, s in enumerate(vcf_order)}
    reindex = np.array([pos[s] for s in sample_order])

    burden, nvar = accumulate_burden(
        args.vcf, args.contig, args.region, tree, vcf_order, args.maf, args.weights
    )
    burden = {g: vec[reindex] for g, vec in burden.items()}  # -> sample_order

    y = pheno["pheno"].to_numpy(dtype=float)
    covars = pheno[covar_cols].to_numpy(dtype=float) if covar_cols else np.empty((len(y), 0))

    results = test_genes(burden, nvar, y, covars, names, args.trait_type, args.min_variants)
    if results.empty:
        log.warning("no testable genes -- on downsampled/MAF-filtered data this is expected")
    else:
        results.to_csv(args.out, sep="\t", index=False)
        log.info("wrote %d genes -> %s", len(results), args.out)
        print(results.head(15).to_string(index=False))
        # Hand-off to INDRA: ground these to HGNC and post to the discovery engine.
        sig = results.loc[results.get("qval", 1) < 0.05, "gene_name"].tolist()
        log.info("%d genes at FDR<0.05 for the INDRA discovery step: %s", len(sig), sig)


if __name__ == "__main__":
    main()