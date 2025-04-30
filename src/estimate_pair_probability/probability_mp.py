import multiprocessing as mp
import pandas as pd
import numpy as np
import sys

def ppi_pair_binom_p(row, pseudo_n):
    """
    We assume the same localisation probability of bait-prey as the prior
    """
    ms_p = row["p_estimation_ms"]
    alpha_prior_ms = ms_p*pseudo_n
    beta_prior_ms  = (1 - ms_p)*pseudo_n
    alpha_post_ms  = alpha_prior_ms + row["n_observed_ms"]
    beta_post_ms   = beta_prior_ms +  row["n_tests_ms"] - row["n_observed_ms"]
    p_est_ms = alpha_post_ms/(alpha_post_ms + beta_post_ms)

    y2h_p = row["p_estimation_y2h"]
    alpha_prior_y2h = y2h_p*pseudo_n
    beta_prior_y2h  = (1 - y2h_p)*pseudo_n
    alpha_post_y2h  = alpha_prior_y2h + row["n_observed_y2h"]
    beta_post_y2h   = beta_prior_y2h +  row["n_tests_y2h"] - row["n_observed_y2h"]
    p_est_y2h = alpha_post_y2h/(alpha_post_y2h + beta_post_y2h)

    p_est_weighted = row["n_tests_ms"]/(row["n_tests_ms"] + row["n_tests_y2h"])*p_est_ms + row["n_tests_y2h"]/(row["n_tests_ms"] + row["n_tests_y2h"])*p_est_y2h

    return pd.concat([
        row,
        pd.Series({
            "weighted_ppi_p": p_est_weighted,
            "ms_ppi_p": p_est_ms,
            "y2h_ppi_p": p_est_y2h
        })])


def _process_chunk(df):
    return df.apply(
        ppi_pair_binom_p,axis=1,args=(1,))

def parallel_process_df(df_chunks):
    with mp.Pool() as pool:
        results = pool.map(_process_chunk, df_chunks)

    return pd.concat(results)

if __name__ == '__main__':
    args = sys.argv

    inferred_data_per_bait_full   = args[1]
    ms_localisation_probability   = args[2]
    y2h_localisation_probability  = args[3]

    pseudo_n  = int(args[4])
    n_cores   = int(args[5])

    pair_probability  = args[6]

    study_count_df = pd.read_csv(inferred_data_per_bait_full, sep="\t")
    ms_localisation = pd.read_csv(ms_localisation_probability, sep="\t")
    ms_localisation = ms_localisation[
        ms_localisation["localisation_bait"] == ms_localisation["localisation_prey"]]
    y2h_localisation = pd.read_csv(y2h_localisation_probability, sep="\t")
    y2h_localisation = y2h_localisation[
        y2h_localisation["localisation_bait"] == y2h_localisation["localisation_prey"]]

    full_localisation = ms_localisation.merge(
        y2h_localisation,
        on=["localisation_bait", "localisation_prey"],
        suffixes=("_ms", "_y2h"))
    full_localisation["localisation"] = full_localisation["localisation_bait"]
    del full_localisation["localisation_bait"]
    del full_localisation["localisation_prey"]

    study_count_df = study_count_df.merge(full_localisation, on="localisation")
    partitions = np.array_split(study_count_df, n_cores)

    probability_estimate_df = parallel_process_df(partitions)

    probability_estimate_df.to_csv(
        pair_probability,
        sep="\t",
        index=False
    )