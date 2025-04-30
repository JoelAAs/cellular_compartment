import numpy as np
import pandas as pd
import multiprocessing as mp

def get_localisation_data(bait_prey_df, localisation_df, prot_to_gene_df):
    bait_prey_df = bait_prey_df.merge(prot_to_gene_df, left_on="bait", right_on="uniprot_id")
    del bait_prey_df["uniprot_id"]
    bait_prey_df = bait_prey_df.merge(prot_to_gene_df, left_on="prey", right_on="uniprot_id", suffixes=("_bait", "_prey"))
    del bait_prey_df["uniprot_id"]

    bait_prey_df = bait_prey_df.merge(localisation_df, left_on="gene_name_bait", right_on="gene_name")
    del bait_prey_df["gene_name"]
    bait_prey_df = bait_prey_df.merge(localisation_df, left_on="gene_name_prey", right_on="gene_name", suffixes=("_bait", "_prey"))
    del bait_prey_df["gene_name"]

    return bait_prey_df


def create_or_update(c_dict, key, value):
    if key in c_dict:
        if value != dict():
            c_dict[key] += value
    else:
        c_dict[key] = value


def infer_n_tests_per_bait_ms(ppi_df, localisation_genes, current_localisation, output_file):
    bait_prey_pair_tests = dict()
    bait_prey_pair_observations = dict()

    localisation_ss = ppi_df[ppi_df["gene_name_bait"].isin(localisation_genes)]
    localisation_ss = localisation_ss[
        ~localisation_ss[["gene_name_bait", "gene_name_prey", "pubmed_id"]].duplicated()
    ]  # Remove multi method as it's a bit inconsistent
    localisation_n_studies = localisation_ss.groupby(
        "gene_name_bait",
        as_index=False)["pubmed_id"].nunique()
    localisation_n_studies = localisation_n_studies.rename({"pubmed_id": "n_studies"}, axis=1)

    for _, bait, n in localisation_n_studies.itertuples():
        create_or_update(bait_prey_pair_tests, bait, dict())
        _ = [create_or_update(
            bait_prey_pair_tests[bait], prey, n) for prey in localisation_genes
        ]

    localisation_ss = localisation_ss[localisation_ss["gene_name_prey"].isin(localisation_genes)]
    localisation_pair_obs = localisation_ss.groupby(["gene_name_bait", "gene_name_prey"], as_index=False).size()
    for _, bait, prey, obs_value in localisation_pair_obs.itertuples():
        create_or_update(bait_prey_pair_observations, bait, dict())
        create_or_update(bait_prey_pair_observations[bait], prey, obs_value)

    with open(output_file, "a+") as w:
        for bait, prey_dict in bait_prey_pair_tests.items():
            for prey, n_tests in prey_dict.items():
                try:
                    observed = bait_prey_pair_observations[bait][prey]
                except KeyError:
                    observed = 0
                lineout = (f"{bait}\t"
                           f"{prey}\t"
                           f"{current_localisation}\t"
                           f"{n_tests}\t"
                           f"{observed}\t"
                           f"MS\n")
                w.write(lineout)


def infer_n_tests_per_bait_y2h(ppi_df, localisation_genes, current_localisation, output_file):
    # NOTE: Could introduce a pid: prey-pool dict:
    bait_prey_pair_tests = dict()
    bait_prey_pair_observations = dict()

    localisation_ss = ppi_df[ppi_df["gene_name_bait"].isin(localisation_genes)]
    localisation_ss = localisation_ss[
        ~localisation_ss[["gene_name_bait", "gene_name_prey", "pubmed_id"]].duplicated()
    ]  # Remove multi method as it's a bit inconsistent

    for pid in localisation_ss["pubmed_id"].unique():
        study_ss = localisation_ss[localisation_ss["pubmed_id"] == pid]
        preys = study_ss["gene_name_prey"].unique()
        baits = study_ss["gene_name_bait"].unique()
        for bait in baits:
            create_or_update(bait_prey_pair_tests, bait, dict())
            _ = [create_or_update(
                bait_prey_pair_tests[bait], prey, 1) for prey in preys
            ]

    observed_localisation_ss = localisation_ss[
        localisation_ss["gene_name_prey"].isin(localisation_genes)
    ]
    for _, bait, prey in observed_localisation_ss[["gene_name_bait", "gene_name_prey"]].itertuples():
        create_or_update(bait_prey_pair_observations, bait, dict())
        create_or_update(bait_prey_pair_observations[bait], prey, 1)

    with open(output_file, "a+") as w:
        for bait in localisation_genes:
            for prey in localisation_genes:
                try:
                    observed = bait_prey_pair_observations[bait][prey]
                except KeyError:
                    observed = 0

                try:
                    n_tests = bait_prey_pair_tests[bait][prey]
                except KeyError:
                    continue # if not tested don't write it out

                lineout = (
                    f"{bait}\t"
                    f"{prey}\t"
                    f"{current_localisation}\t"
                    f"{n_tests}\t"
                    f"{observed}\t"
                    f"Y2H\n")
                w.write(lineout)

rule n_tests_per_baits:
    params:
        min_annotated_genes = 200, # arbitrary, but if the assumption of inter localisation interactions being "more" common a higher number is needed
        ms_methods = [
            "MI-0006",
            "MI-0007",
            "MI-0096",
            "MI-0676",
            "MI-1314"
        ],
        y2h_methods = [
            "MI-0397",  # (two hybrid array)
            "MI-1112",  # (two hybrid prey pooling approach)
            "MI-0398",  # (two hybrid pooling approach)
            "MI-0018",  # (two hybrid)
            "MI-0399",  # (two hybrid fragment pooling approach) NOTE: One-to-one
            "MI-1356",  # (validated two hybrid)
            "MI-2215",  # (barcode fusion genetics two hybrid)
            "MI-1113"  # (two hybrid bait and prey pooling approach)
        ]
    input:
        bait_preys = "data/bait_prey_publications.csv",
        localisation_data = "data/localisation/gene_to_location.csv",
        uniprot_to_gene = "data/uniprot_to_gene_name.csv"
    output:
        inferred_data_per_bait_ms   = "work_folder/inferred_negative_localisation_ms.csv",
        inferred_data_per_bait_y2h  = "work_folder/inferred_negative_localisation_y2h.csv",
        inferred_data_per_bait_full = "work_folder/inferred_negative_localisation.csv"
    run:
        bait_prey_df = pd.read_csv(input.bait_preys,sep="\t")
        localisation_df = pd.read_csv(input.localisation_data,sep="\t")
        prot_to_gene_df = pd.read_csv(input.uniprot_to_gene,sep="\t")

        bait_prey_df = get_localisation_data(
            bait_prey_df,
            localisation_df,
            prot_to_gene_df
        )

        localisations_considered = localisation_df.groupby("localisation", as_index=False).size()
        localisations_considered = localisations_considered[
            localisations_considered["size"] > params.min_annotated_genes]["localisation"].tolist()

        bait_prey_df_ms  = bait_prey_df[bait_prey_df["detection_method"].isin(params.ms_methods)]
        bait_prey_df_y2h = bait_prey_df[bait_prey_df["detection_method"].isin(params.y2h_methods)]

        header = "\t".join([
            "gene_name_bait",
            "gene_name_prey",
            "localisation",
            "n_tests",
            "n_observed",
            "detection_method"
        ]) + "\n"
        with open(output.inferred_data_per_bait_ms, "w") as w:
            w.write(header)
        with open(output.inferred_data_per_bait_y2h, "w") as w:
            w.write(header)

        for current_localisation in localisations_considered:
            localisation_genes = localisation_df[
                localisation_df["localisation"] == current_localisation]["gene_name"].tolist()
            localisation_ms_ss  = bait_prey_df_ms[bait_prey_df_ms["localisation_bait"] == current_localisation].copy()
            localisation_y2h_ss = bait_prey_df_y2h[bait_prey_df_y2h["localisation_bait"] == current_localisation].copy()

            infer_n_tests_per_bait_ms(localisation_ms_ss, localisation_genes, current_localisation, output.inferred_data_per_bait_ms)
            infer_n_tests_per_bait_y2h(localisation_y2h_ss, localisation_genes, current_localisation, output.inferred_data_per_bait_y2h)

        ms_df = pd.read_csv(output.inferred_data_per_bait_ms, sep="\t")
        y2h_df = pd.read_csv(output.inferred_data_per_bait_y2h,sep="\t")

        full_df = ms_df.merge(
            y2h_df, on = ["gene_name_bait", "gene_name_prey", "localisation"], how="outer", suffixes=("_ms", "_y2h")
        ).fillna(0)
        full_df.to_csv(output.inferred_data_per_bait_full, sep="\t", index=False)


rule estimate_probability_of_prey:
    params:
        pseudo_n = 5,
        n_cores = 20
    input:
        inferred_data_per_bait_full = "work_folder/inferred_negative_localisation.csv",
        ms_localisation_probability = "work_folder_MS/localisation/any/beta_estimation.csv",
        y2h_localisation_probability = "work_folder_Y2H/localisation/yeast/beta_estimation.csv"
    output:
        pair_probability = "work_folder/ppi_localisation_prior_probability.csv"
    shell:
        """
        python src/estimate_pair_probability/probability_mp.py \
            {input.inferred_data_per_bait_full} \
            {input.ms_localisation_probability} \
            {input.y2h_localisation_probability} \
            {params.pseudo_n} \
            {params.n_cores} \
            {output.pair_probability}
        """

