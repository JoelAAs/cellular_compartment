from scipy.stats import fisher_exact
from get_localisation_binom_per_method import get_localisation_data, get_prior_information, get_beta_posterior_values
import pandas as pd

rule get_yeast_detection:
    input:
        bait_prey_data = "data/bait_prey_publications.csv",
        localisation_data = "data/localisation/gene_to_location.csv",
        uniprot_to_gene_name_data = "data/uniprot_to_gene_name.csv"
    output:
        localisation = "work_folder_Y2H/cell_line_subset/yeast_bait_prey_publications.csv"
    run:
        intact_df = pd.read_csv(input.bait_prey_data,sep="\t")
        localisation_df = pd.read_csv(input.localisation_data,sep="\t")
        prot_to_gene_df = pd.read_csv(input.uniprot_to_gene_name_data, sep="\t")

        intact_df = intact_df[intact_df["detection_method"].isin(config["methods_considered"])]

        intact_df = intact_df.merge(prot_to_gene_df, left_on="bait", right_on="uniprot_id")
        del intact_df["uniprot_id"]
        intact_df = intact_df.merge(prot_to_gene_df, left_on="prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))

        intact_df = intact_df.merge(localisation_df,left_on="gene_name_bait",right_on="gene_name")
        del intact_df["gene_name"]
        intact_df = intact_df.merge(localisation_df,left_on="gene_name_prey",right_on="gene_name",suffixes=(
        "_bait", "_prey"))
        intact_df["match"] = intact_df["localisation_bait"] == intact_df["localisation_prey"]
        intact_df["cl_id"] = "yeast"
        intact_df.to_csv(
            output.localisation,
            sep="\t", index=False
        )

rule get_any_MS_detection:
    input:
        bait_prey_data = "data/bait_prey_publications.csv",
        localisation_data = "data/localisation/gene_to_location.csv",
        uniprot_to_gene_name_data = "data/uniprot_to_gene_name.csv"
    output:
        localisation = "work_folder_MS/cell_line_subset/any_bait_prey_publications.csv"
    run:
        intact_df = pd.read_csv(input.bait_prey_data,sep="\t")
        localisation_df = pd.read_csv(input.localisation_data,sep="\t")
        prot_to_gene_df = pd.read_csv(input.uniprot_to_gene_name_data, sep="\t")

        intact_df = intact_df[intact_df["detection_method"].isin(config["methods_considered"])]

        intact_df = intact_df.merge(prot_to_gene_df, left_on="bait", right_on="uniprot_id")
        del intact_df["uniprot_id"]
        intact_df = intact_df.merge(prot_to_gene_df, left_on="prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))

        intact_df = intact_df.merge(localisation_df,left_on="gene_name_bait",right_on="gene_name")
        del intact_df["gene_name"]
        intact_df = intact_df.merge(localisation_df,left_on="gene_name_prey",right_on="gene_name",suffixes=(
        "_bait", "_prey"))
        intact_df["match"] = intact_df["localisation_bait"] == intact_df["localisation_prey"]
        intact_df["cl_id"] = "any"
        intact_df.to_csv(
            output.localisation,
            sep="\t", index=False
        )


rule get_CL_specific:
    input:
        bait_prey_data = "data/CL_annotated_bait_prey.csv",
        localisation_data = "data/localisation/gene_to_location.csv",
    output:
        localisation = "work_folder_{project}/cell_line_subset/{cell_line}_bait_prey_publications.csv"
    run:
        ppi_localisation_df = get_localisation_data(
            input.bait_prey_data,
            input.localisation_data,
            config["methods_considered"]
        )
        ppi_localisation_df[ppi_localisation_df["cl_id"] == wildcards.cell_line].to_csv(
            output.localisation,
            sep="\t", index=False
        )




rule get_localisation_beta_estimate:
    params:
        pseudo_n = 10
    input:
        bait_prey_data = "work_folder_{project}/cell_line_subset/{cell_line}_bait_prey_publications.csv",
        localisation_data = "data/localisation/gene_to_location.csv"
    output:
        localisation = "work_folder_{project}/localisation/{cell_line}/beta_estimation.csv"
    run:
        ppi_localisation_df = pd.read_csv(input.bait_prey_data, sep="\t")
        localisation_df = pd.read_csv(input.localisation_data, sep="\t")
        localisations_baits = ppi_localisation_df["localisation_bait"].unique()
        localisations_prey = ppi_localisation_df["localisation_prey"].unique()

        with open(output.localisation, "w") as w:
            w.write("\t".join([
                "localisation_bait", "localisation_prey", "p_estimation",
                "low_0.95", "high_0.95", "n_ppis",
                "n_unique_bait", "normalised_entropy",
                "cell_line"
            ]) + "\n")
            for bait_localisation in localisations_baits:
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
                        f"{h_norm}\t"
                        f"{wildcards.cell_line}\n")


rule aggregate_beta_estimates:
    input:
        p_estimations = expand(
            "work_folder_{{project}}/localisation/{cell_line}/beta_estimation.csv", cell_line = config["cell_lines"]
        )
    output:
        all_localisation_beta = "work_folder_{project}/localisation/beta_estimation.csv"
    run:
        with open(output.all_localisation_beta, "w") as w:
            w.write("\t".join([
                "localisation_bait", "localisation_prey", "p_estimation",
                "low_0.95", "high_0.95", "n_ppis",
                "n_unique_bait", "normalised_entropy", "cell_line"
            ]) + "\n")
        for estimation_file in input.p_estimations:
            shell("tail -n+2 {estimation_file} >> {output.all_localisation_beta}")


rule get_localisation_overlap:
    input:
        localisation_data = "data/localisation/gene_to_location.csv"
    output:
        localisation_overlap = "work_folder_{project}/localisation_overlap.csv"
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
        bait_prey_data = "work_folder_{project}/cell_line_subset/{cell_line}_bait_prey_publications.csv"
    output:
        fisher_exact = "work_folder_{project}/localisation/{cell_line}/fisher_exact.csv"
    run:
        ppi_localisation_df = pd.read_csv(input.bait_prey_data,sep="\t")
        localisations_baits = ppi_localisation_df["localisation_bait"].unique()
        localisations_prey = ppi_localisation_df["localisation_prey"].unique()

        with open(output.fisher_exact,"w") as w:
            w.write("\t".join([
                "localisation_bait", "localisation_prey", "p_value", "n_ppi", "n_ppi_other", "statistic", "cell_line"
            ]) + "\n")
            for bait_localisation in localisations_baits:
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

                    n_ppi_other = len(ppi_localisation_df_other)
                    n_ppi_local = len(ppi_localisation_df_local)


                    table = [
                        [prey_localisation_hits_other, n_ppi_other-prey_localisation_hits_other],
                        [prey_localisation_hits_local, n_ppi_local-prey_localisation_hits_local]
                        ]
                    res = fisher_exact(table)
                    w.write(f"{bait_localisation}\t"
                            f"{prey_localisation}\t"
                            f"{res.pvalue}\t"
                            f"{n_ppi_local}\t"
                            f"{n_ppi_other}\t"
                            f"{res.statistic}\t"
                            f"{wildcards.cell_line}\n")


rule aggregate_fisher:
    input:
        p_estimations = expand(
            "work_folder_{{project}}/localisation/{cell_line}/fisher_exact.csv", cell_line = config["cell_lines"]
        )
    output:
        all_localisation_beta = "work_folder_{project}/localisation/fisher_exact.csv"
    run:
        with open(output.all_localisation_beta, "w") as w:
            w.write("\t".join([
                "localisation_bait", "localisation_prey", "p_value", "n_ppi", "n_ppi_other", "statistic", "cell_line"
            ]) + "\n")
        for estimation_file in input.p_estimations:
            shell("tail -n+2 {estimation_file} >> {output.all_localisation_beta}")