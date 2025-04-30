from get_localisation_binom_per_method import get_prior_information, get_beta_posterior_values

import pandas as pd

rule get_method_localisation_subsets:
    input:
        bait_prey_data = "work_folder_{project}/cell_line_subset/{cell_line}_bait_prey_publications.csv",
    output:
        expand(
            "work_folder_{{project}}/method_localisation/{{cell_line}}/{method}_{localisation}_ppi.csv",
            method = config["methods_considered"],
            localisation = config["localisations_considered"]
        )
    run:
        ppi_localisation_df = pd.read_csv(input.bait_prey_data,sep="\t")
        for current_method in config["methods_considered"]:
            for current_localisation in config["localisations_considered"]:
                ss_localisation_df = ppi_localisation_df[
                    (ppi_localisation_df["detection_method"] == current_method) &
                    (ppi_localisation_df["localisation_bait"] == current_localisation)
                ]
                ss_localisation_df.to_csv(
                    f"work_folder_{wildcards.project}/method_localisation/{wildcards.cell_line}/{current_method}_{current_localisation}_ppi.csv",
                    sep="\t", index=False)

rule get_probability_of_match:
    params:
        pseudo_n = 1000
    input:
        mehtod_local_ss = "work_folder_{project}/method_localisation/{cell_line}/{method}_{localisation}_ppi.csv",
        localisation_data = "data/localisation/gene_to_location.csv"
    output:
        method_local_pval = "work_folder_{project}/method_localisation/{cell_line}/beta_estimation/{method}_{localisation}.csv"
    run:
        method_localisation_ss = pd.read_csv(input.mehtod_local_ss, sep="\t")
        localisation_df = pd.read_csv(input.localisation_data, sep="\t")

        alpha_prior, beta_prior = get_prior_information(localisation_df, wildcards.localisation, params.pseudo_n)
        p_est, cred_int, n, k, h_norm = get_beta_posterior_values(method_localisation_ss, alpha_prior, beta_prior)

        with open(output.method_local_pval, "w") as w:
            w.write(f"{wildcards.method}\t"
                    f"{wildcards.localisation}\t"
                    f"{p_est}\t"
                    f"{cred_int[0]}\t"
                    f"{cred_int[1]}\t"
                    f"{n}\t"
                    f"{k}\t"
                    f"{h_norm}\n")

rule aggregate_p_estimations:
    input:
        p_estimations = expand(
            "work_folder_{{project}}/method_localisation/{{cell_line}}/beta_estimation/{method}_{localisation}.csv",
            method = config["methods_considered"],
            localisation = config["localisations_considered"]
        )
    output:
        aggregated_values = "work_folder_{project}/method_localisation/{cell_line}/beta_estimation/method_localisation_p_estimation.csv"
    run:
        with open(output.aggregated_values, "w") as w:
            w.write("\t".join([
                "detection_method", "localisation_bait", "p_estimation",
                "low_0.95", "high_0.95", "n_ppis",
                "n_unique_bait", "normalised_entropy"
            ]) + "\n")
        for estimation_file in input.p_estimations:
            shell("cat {estimation_file} >> {output.aggregated_values}")



