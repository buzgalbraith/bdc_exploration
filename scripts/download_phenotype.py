from gen3helper import gen3Client
import subprocess
import pandas as pd
import re


COMMONS = 'https://gen3.biodatacatalyst.nhlbi.nih.gov'
CRED_FILE = '/Users/buzgalbraith/.gen3/credentials.json'
target_project = 'tutorial-synthetic_data_set_1'
oid = 'dg.4503/7be3c6bb-602c-402a-88d2-50394bf8b433'

if __name__ == "__main__":
    client = gen3Client(endpoint=COMMONS, credential_file=CRED_FILE)
    # client.get_project_files()
    client.download_files(file_ids=[oid], save_directory='./data')
    # subprocess.call(
    #     [
    #         'unzip', 'data/tutorial-synthetic_data_set_1_structured_data.zip', '-d', 'data/'
    #     ]
    # )
    # lab = pd.read_csv("data/lab_result.tsv", sep="\t")
    # dem = pd.read_csv("data/demographic.tsv", sep="\t")
    # # submitter_id embeds the 1000G sample id: 'HG00325_lab_res' -> 'HG00325'
    # samp = lambda s: re.match(r"(HG\d+|NA\d+)", str(s)).group(1)
    # lab["sample_id"] = lab["submitter_id"].map(samp)
    # dem["sample_id"] = dem["submitter_id"].map(samp)

    # TRAIT = "ldl"           # pick from lab_result; confirm the exact name in schema.json
    # pheno = (lab[["sample_id", TRAIT,]]
    #         .merge(dem[["sample_id", 'annotated_sex', 'population', 'age_at_index']], on="sample_id", how="inner")
    #         )
    #         #  .rename(columns={TRAIT: "phenotype"})
    #         #  .dropna(subset=["phenotype"]))
    # pheno["annotated_sex"] = pheno["annotated_sex"].map({"male": 1, "female": 0})
    # pheno = pd.get_dummies(pheno, columns=["population"], drop_first=True, dtype=float)
    # pheno.to_csv("data/phenotypes.csv", index=False)