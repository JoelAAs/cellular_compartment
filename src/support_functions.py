import pandas as pd

def read_ppi(intact, localisation_annotations, uniprot_gene_name):
    intact_df = pd.read_csv(intact,sep="\t")
    localisation_df = pd.read_csv(localisation_annotations,sep="\t")
    localisation_df = localisation_df[["source", "target_desc"]]
    gene_name_to_uniprot = pd.read_csv(uniprot_gene_name,sep="\t")

    intact_df = intact_df.merge(
        gene_name_to_uniprot,
        left_on="bait",right_on="uniprot_id"
    )

    intact_df = intact_df.merge(
        gene_name_to_uniprot,
        left_on="prey",right_on="uniprot_id",
        suffixes=("_bait", "_prey")
    )
    intact_df = intact_df.merge(
        localisation_df,
        left_on="gene_name_bait",right_on="source",
    )
    intact_df = intact_df.merge(
        localisation_df,
        left_on="gene_name_prey",right_on="source",
        suffixes=("_bait", "_prey")
    )
    intact_df = intact_df[[
        "gene_name_bait",
        "gene_name_prey",
        "detection_method",
        "target_desc_bait",
        "target_desc_prey"]
    ]
    return intact_df