import pandas as pd

"""
This script processes the gnomAD v4.1 constraint scores present here /hpf/largeprojects/ccmbio/ajain/crg2_hg38/gnomad_v4.1/constraint_metrics/gnomad.v4.1.constraint_metrics.tsv
to retain only Ensemble mane select genes.

Final processed file is located in ~/cre/data/gnomad_scores_v4.1.csv

Script written by Anjali Jain
Date: Aug 27, 2024
"""

# Read in the constraints file that was downloaded from gnomAD: https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
constraints = pd.read_csv(
    "/hpf/largeprojects/ccmbio/ajain/crg2_hg38/gnomad_v4.1/constraint_metrics/gnomad.v4.1.constraint_metrics.tsv",
    sep="\t",
)

# Filtering for the Ensembl canonical transcript
constraints = constraints[
    constraints["gene_id"].str.startswith("ENSG") & constraints["mane_select"] == True
]

# Check if there are any duplicates
if len(constraints["gene_id"]) == len(set(constraints["gene_id"])):
    print(
        f"There are no duplicates in this dataframe. There are {len(constraints['gene_id'])} unique ensembl mane select genes"
    )
else:
    print("Duplicates found!")

# Keeping columns that are relevant
constraints = constraints[
    [
        "gene_id",
        "transcript",
        "lof.oe",
        "lof.oe_ci.lower",
        "lof.oe_ci.upper",
        "lof.pLI",
        "lof.pNull",
        "lof.pRec",
        "mis.oe",
        "mis.z_score",
    ]
]

# renaming the columns
constraints = constraints.rename(
    columns={
        "gene_id": "Ensembl_gene_id",
        "transcript": "Ensembl_transcript_id",
        "lof.oe": "Gnomad_oe_lof_score",
        "lof.oe_ci.lower": "Gnomad_oe_ci_lower",
        "lof.oe_ci.upper": "Gnomad_oe_ci_upper",
        "lof.pLI": "Gnomad_pLI_score",
        "lof.pNull": "Gnomad_pnull_score",
        "lof.pRec": "Gnomad_prec_score",
        "mis.oe": "Gnomad_oe_mis_score",
        "mis.z_score": "Gnomad_mis_z_score",
    }
)

# Write the constraints score to ~/cre/data/gnomad_scores_v4.1.csv
constraints.to_csv("~/cre/data/gnomad_scores_v4.1.csv", index=False)
