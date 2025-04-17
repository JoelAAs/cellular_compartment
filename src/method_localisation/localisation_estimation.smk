from scipy.stats import fisher_exact

from get_localisation_binom_per_method import get_localisation_data, get_prior_information, get_beta_posterior_values
import pandas as pd

rule get_localisation_beta_estimate:
    params:
        pseudo_n = 10
    input:
        bait_prey_data = "data/bait_prey_publications.csv",
        localisation_data = "data/localisation/gene_to_location.csv",
        uniprot_to_gene_name_data = "data/uniprot_to_gene_name.csv"
    output:
        localisation = "work_folder/localisation/beta_estimation.csv"
    run:
        ppi_localisation_df = get_localisation_data(
            input.bait_prey_data,
            input.localisation_data,
            input.uniprot_to_gene_name_data,
            config["methods_considered"]
        )

        localisation_df = pd.read_csv(input.localisation_data, sep="\t")
        localisations_baits = ppi_localisation_df["localisation_bait"].unique()
        localisations_prey = ppi_localisation_df["localisation_prey"].unique()

        with open(output.localisation, "w") as w:
            w.write("\t".join([
                "localisation_bait", "localisation_prey", "p_estimation",
                "low_0.95", "high_0.95", "n_ppis",
                "n_unique_bait", "normalised_entropy"
            ]) + "\n")
            for bait_localisation in localisations_baits:
                print(bait_localisation)
                ppi_bait_ss = ppi_localisation_df[
                    ppi_localisation_df["localisation_bait"] == bait_localisation
                ].copy()
                for prey_localisation in localisations_prey:
                    ppi_bait_ss["match"] = ppi_bait_ss["localisation_prey"] == prey_localisation

                    alpha_prior, beta_prior = get_prior_information(localisation_df, prey_localisation, params.pseudo_n)
                    p_est, cred_int, n, k, h_norm = get_beta_posterior_values(ppi_bait_ss, alpha_prior, beta_prior)

                    w.write(f"{bait_localisation}\t"
                        f"{prey_localisation}\t"
                        f"{p_est}\t"
                        f"{cred_int[0]}\t"
                        f"{cred_int[1]}\t"
                        f"{n}\t"
                        f"{k}\t"
                        f"{h_norm}\n")


rule get_localisation_overlap:
    input:
        localisation_data = "data/localisation/gene_to_location.csv"
    output:
        localisation_overlap = "work_folder/localisation_overlap.csv"
    run:
        localisation_df = pd.read_csv(input.localisation_data, sep="\t")

        all_localisations = localisation_df["localisation"].unique().tolist()
        with open(output.localisation_overlap, "w") as w:
            w.write("localisation_a\tlocalisation_b\tiou\n")
            for from_localisation in all_localisations:
                for to_localisation in all_localisations:
                    genes_from = set(localisation_df[
                        localisation_df["localisation"] == from_localisation
                    ]["gene_name"].tolist())

                    genes_to = set(localisation_df[
                        localisation_df["localisation"] == to_localisation
                    ]["gene_name"].tolist())

                    iou = len(genes_from & genes_to)/len(genes_from | genes_to)
                    w.write(f"{from_localisation}\t{to_localisation}\t{iou}\n")




rule get_tests_against_background:
    input:
        bait_prey_data = "data/bait_prey_publications.csv",
        localisation_data = "data/localisation/gene_to_location.csv",
        uniprot_to_gene_name_data = "data/uniprot_to_gene_name.csv"
    output:
        fisher_exact = "work_folder/localisation/fisher_exact.csv"
    run:
        ppi_localisation_df = get_localisation_data(
            input.bait_prey_data,
            input.localisation_data,
            input.uniprot_to_gene_name_data,
            config["methods_considered"]
        )

        localisation_df = pd.read_csv(input.localisation_data,sep="\t")
        localisations_baits = ppi_localisation_df["localisation_bait"].unique()
        localisations_prey = ppi_localisation_df["localisation_prey"].unique()

        with open(output.fisher_exact,"w") as w:
            w.write("\t".join([
                "localisation_bait", "localisation_prey", "p_value", "statistic"
            ]) + "\n")
            for bait_localisation in localisations_baits:
                print(bait_localisation)
                ppi_localisation_df_other = ppi_localisation_df[
                    ppi_localisation_df["localisation_bait"] != bait_localisation
                    ].copy()
                ppi_localisation_df_local = ppi_localisation_df[
                    ppi_localisation_df["localisation_bait"] == bait_localisation
                    ].copy()

                for prey_localisation in localisations_prey:
                    ppi_localisation_df_other["match"] = ppi_localisation_df_other["localisation_prey"] == prey_localisation
                    ppi_localisation_df_local["match"] = ppi_localisation_df_local["localisation_prey"] == prey_localisation

                    prey_localisation_hits_other = ppi_localisation_df_other["match"].sum()
                    prey_localisation_hits_local = ppi_localisation_df_local["match"].sum()


                    table = [
                        [prey_localisation_hits_other, len(ppi_localisation_df_other)-prey_localisation_hits_other],
                        [prey_localisation_hits_local, len(ppi_localisation_df_local)-prey_localisation_hits_local]
                        ]
                    print(table)
                    res = fisher_exact(table)
                    w.write(f"{bait_localisation}\t"
                            f"{prey_localisation}\t"
                            f"{res.pvalue}\t"
                            f"{res.statistic}\n")