import pandas as pd
from scipy.stats import beta

def get_localisation_data(bait_prey_data, localisation_data, uniprot_to_gene_name_data):
    intact_df = pd.read_csv(bait_prey_data, sep = "\t")
    localisation_df = pd.read_csv(localisation_data, sep = "\t")
    prot_to_gene_df = pd.read_csv(uniprot_to_gene_name_data, sep="\t")
    print(intact_df.columns)

    intact_df = intact_df.merge(prot_to_gene_df, left_on="bait", right_on="uniprot_id")
    del intact_df["uniprot_id"]
    intact_df = intact_df.merge(prot_to_gene_df, left_on="prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
    print(intact_df.columns)

    intact_df = intact_df.merge(localisation_df, left_on="gene_name_bait", right_on="gene_name")
    print(intact_df.columns)
    del intact_df["gene_name"]
    intact_df = intact_df.merge(localisation_df, left_on="gene_name_prey", right_on="gene_name", suffixes=("_bait", "_prey"))
    intact_df["match"] = intact_df["localisation_bait"] == intact_df["localisation_prey"]

    return intact_df


def get_beta_posterior_values(test_df, alpha_prior, beta_prior):
    # TODO: Add 
    n = test_df.shape[0]
    x = test_df["match"].sum()

    alpha_post = alpha_prior + x
    beta_post = beta_prior + n - x
    p_est = alpha_post/(alpha_post+beta_post)
    cred_interval = beta.ppf([0.025, 0.975], alpha_post, beta_post)

    return p_est, cred_interval


def get_prior_information(localisation_df, current_localisation, pseudo_n=1000):
    alpha_prior = localisation_df.groupby("localisation").size()/localisation_df.shape[0]*pseudo_n
    alpha_prior = round(alpha_prior[current_localisation])
    beta_prior = pseudo_n - alpha_prior
    return alpha_prior, beta_prior


