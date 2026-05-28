import pandas as pd, re

## load in ## 
lab = pd.read_csv("lab_result.tsv", sep="\t")
dem = pd.read_csv("demographic.tsv", sep="\t")

# submitter_id embeds the 1000G sample id: 'HG00325_lab_res' -> 'HG00325'
samp = lambda s: re.match(r"(HG\d+|NA\d+)", str(s)).group(1)
lab["sample_id"] = lab["submitter_id"].map(samp)
dem["sample_id"] = dem["submitter_id"].map(samp)

TRAIT = "ldl"           # pick from lab_result; confirm the exact name in schema.json
pheno = (lab[["sample_id", TRAIT,]]
         .merge(dem[["sample_id", 'annotated_sex']], on="sample_id", how="inner")
        )
        #  .rename(columns={TRAIT: "phenotype"})
        #  .dropna(subset=["phenotype"]))
pheno["annotated_sex"] = pheno["annotated_sex"].map({"male": 1, "female": 0})
pheno.to_csv("pheno_chr17.csv", index=False)