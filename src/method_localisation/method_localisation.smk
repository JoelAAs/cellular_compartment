from get_localisation_binom_per_method import get_localisation_data, get_prior_information, get_beta_posterior_values
import pandas as pd

rule get_method_localisation_subsets:
    input:
        bait_prey_data = "data/bait_prey_publications.csv",
        localisation_data = "data/localisation/gene_to_location.csv",
        uniprot_to_gene_name_data = "data/uniprot_to_gene_name.csv"
    output:
        expand(
            "work_folder/method_localisation/{method}_{localisation}_ppi.csv",
            method = config["methods_considered"],
            localisation = config["localisations_considered"]
        )
    run:
        localisation_df = get_localisation_data(
            input.bait_prey_data,
            input.localisation_data,
            input.uniprot_to_gene_name_data)

        for current_method in config["methods_considered"]:
            for current_localisation in config["localisations_considered"]:
                ss_localisation_df = localisation_df[
                    (localisation_df["detection_method"] == current_method) &
                    (localisation_df["localisation_bait"] == current_localisation)
                ]
                ss_localisation_df.to_csv(
                    f"work_folder/method_localisation/{current_method}_{current_localisation}_ppi.csv",
                    sep="\t", index=False)

rule get_probability_of_match:
    params:
        pseudo_n = 1000
    input:
        mehtod_local_ss = "work_folder/method_localisation/{method}_{localisation}_ppi.csv",
        localisation_data = "data/localisation/gene_to_location.csv"
    output:
        method_local_pval = "work_folder/method_localisation/beta_estimation/{method}_{localisation}.csv"
    run:
        method_localisation_ss = pd.read_csv(input.mehtod_local_ss, sep="\t")
        localisation_df = pd.read_csv(input.localisation_data, sep="\t")

        alpha_prior, beta_prior = get_prior_information(localisation_df, wildcards.localisation, params.pseudo_n)
        p_est, cred_int = get_beta_posterior_values(method_localisation_ss, alpha_prior, beta_prior)

        with open(output.method_local_pval, "w") as w:
            w.write(f"{wildcards.method}\t{wildcards.localisation}\t{p_est}\t{cred_int[0]}\t{cred_int[1]}\n")

rule aggregate_p_estimations:
    input:
        p_estimations = expand(
            "work_folder/method_localisation/beta_estimation/{method}_{localisation}.csv",
            method = config["methods_considered"],
            localisation = config["localisations_considered"]
        )
    output:
        aggregated_values = "work_folder/method_localisation/beta_estimation/method_localisation_p_estimation.csv"
    run:
        with open(output.aggregated_values, "w") as w:
            w.write("\t".join(["detection_method", "localisation_bait", "p_estimation", "low_0.95", "high_0.95"]) + "\n")
        for estimation_file in input.p_estimations:
            shell("cat {estimation_file} >> {output.aggregated_values}")



