import pandas as pd
from scipy.stats import beta
import math

def get_localisation_data(bait_prey_data, localisation_data, uniprot_to_gene_name_data, detection_methods):
    intact_df = pd.read_csv(bait_prey_data, sep = "\t")
    localisation_df = pd.read_csv(localisation_data, sep = "\t")
    prot_to_gene_df = pd.read_csv(uniprot_to_gene_name_data, sep="\t")

    intact_df = intact_df[intact_df["detection_method"].isin(detection_methods)]

    intact_df = intact_df.merge(prot_to_gene_df, left_on="bait", right_on="uniprot_id")
    del intact_df["uniprot_id"]
    intact_df = intact_df.merge(prot_to_gene_df, left_on="prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))

    intact_df = intact_df.merge(localisation_df, left_on="gene_name_bait", right_on="gene_name")
    del intact_df["gene_name"]
    intact_df = intact_df.merge(localisation_df, left_on="gene_name_prey", right_on="gene_name", suffixes=("_bait", "_prey"))
    intact_df["match"] = intact_df["localisation_bait"] == intact_df["localisation_prey"]

    return intact_df


def get_beta_posterior_values(test_df, alpha_prior, beta_prior):
    n = test_df.shape[0]
    x = test_df["match"].sum()

    p = test_df.groupby("gene_name_bait").size()/n
    k = len(p)
    if k > 1:
        h_norm = -sum(map(lambda x: x*math.log2(x), p))/math.log2(k)  # normalised entropy
    else:
        h_norm = math.nan

    alpha_post = alpha_prior + x
    beta_post = beta_prior + n - x
    p_est = alpha_post/(alpha_post+beta_post)
    cred_interval = beta.ppf([0.025, 0.975], alpha_post, beta_post)

    return p_est, cred_interval, n, k, h_norm


def get_prior_information(localisation_df, current_localisation, pseudo_n=10):
    # TODO check_overlap
    alpha_min = 1/localisation_df.shape[0]
    alpha_prior = localisation_df.groupby("localisation").size()/localisation_df.shape[0]*pseudo_n
    alpha_prior = alpha_prior[current_localisation]
    if alpha_prior < alpha_min:
        alpha_prior = alpha_min
    beta_prior = pseudo_n - alpha_prior
    return alpha_prior, beta_prior
