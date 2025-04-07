import pandas as pd
import numpy as np
from scipy.stats import norm
from support_functions import read_ppi
n_permutations = config["n_permutations"]


def get_prey_localisation_likelihood(bait_prey_df):
    bait_prey_counts = bait_prey_df.groupby(["localisation_bait", "localisation_prey"], as_index=False).size()
    bait_sum = bait_prey_counts.groupby("localisation_bait", as_index=False)["size"].sum()
    bait_prey_counts = bait_prey_counts.merge(bait_sum, on="localisation_bait", suffixes=("_bait_prey_count", "_bait_count"))
    bait_prey_counts["likelihood_prey"] = bait_prey_counts["size_bait_prey_count"].div(bait_prey_counts["size_bait_count"])
    return bait_prey_counts



rule bioid_per_bait_localisation:
    input:
        selected_df="work_folder/selected_baits.csv"
    output:
        bioid_file = "work_folder/localisation_probability/bioID_localisation.csv"
    run:
        intact_df = pd.read_csv(input.selected_df,sep="\t")

        bioid_ms_ss = intact_df[intact_df["detection_method"] == "MI-1314"]

        permutation_prey_counts = get_prey_localisation_likelihood(bioid_ms_ss)
        permutation_prey_counts.to_csv(
            output.bioid_file,
            sep = "\t", index = False
        )

rule permute_per_bait_localisation:
    params:
        frac = 0.9,
        other_ms_methods = [
            "MI-0006",
            "MI-0007",
            "MI-0096",
            "MI-0004",
            "MI-0019"
            ]
    input:
        selected_df="work_folder/selected_baits.csv"
    output:
        permutation_sets = expand(
            "work_folder/localisation_probability/permutation/{{localisation}}_0.9_set_{n}.csv",
            n = range(n_permutations)
        )
    run:
        intact_df = pd.read_csv(input.selected_df,sep="\t")

        other_ms_ss = intact_df[intact_df["detection_method"].isin(params.other_ms_methods)]
        other_ms_localisation_ss = other_ms_ss[other_ms_ss["localisation_bait"] == wildcards.localisation]

        for i, permutation_file in enumerate(output.permutation_sets):
            sample_df = other_ms_localisation_ss.sample(frac=params.frac)
            permutation_prey_counts = get_prey_localisation_likelihood(sample_df)
            permutation_prey_counts["permutation"] = i
            permutation_prey_counts.to_csv(
                permutation_file,
                sep = "\t", index = False
            )


rule estimate_quant_localisation:
    input:
        biotin_file = "work_folder/localisation_probability/bioID_localisation.csv",
        localisation_bait_probability = expand(
            "work_folder/localisation_probability/permutation/{{localisation}}_0.9_set_{n}.csv",
            n = range(n_permutations)
        )
    output:
        biotid_cdf_quant = "work_folder/localisation_probability/bioid_quantile/bait_{localisation}.csv"
    run:
        localisation_df_list = [
            pd.read_csv(
                file,
                sep="\t"
            ) for file in input.localisation_bait_probability
        ]
        localisation_df = pd.concat(localisation_df_list)

        biotin_df = pd.read_csv(
            input.biotin_file,
            sep="\t"
        )
        with open(output.biotid_cdf_quant, "w") as w:
            w.write("localisation_bait\tlocalisation_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tfrac_in_permutation\tin_bioid\n")
            for current_localisation in config["localisation"]:
                biotin_observed = True
                all_probabilities = np.zeros(n_permutations)
                ss_df = localisation_df[localisation_df["localisation_prey"] == current_localisation]
                permuted = np.array(ss_df["likelihood_prey"])

                frac_observed = len(permuted)/n_permutations
                if len(permuted) != 0:
                    all_probabilities[:len(permuted)] = permuted

                observed_value = biotin_df.loc[
                                 (biotin_df["localisation_prey"] == current_localisation) &
                                 (biotin_df["localisation_bait"] == wildcards.localisation)
                ]["likelihood_prey"].values
                if len(observed_value) > 0:
                    observed_value=observed_value[0]
                else:
                    observed_value = 0
                    biotin_observed = False


                mu = np.mean(all_probabilities)
                std = np.std(all_probabilities)
                quant_value = norm.cdf(observed_value, mu, std)
                w.write(
                    f"{wildcards.localisation}\t"
                    f"{current_localisation}\t"
                    f"{mu}\t"
                    f"{std}\t"
                    f"{quant_value}\t"
                    f"{observed_value}\t"
                    f"{frac_observed}\t"
                    f"{biotin_observed}\n"
                )


rule aggregate_quant_data:
    input:
        all_baits_probs = expand(
            "work_folder/localisation_probability/bioid_quantile/bait_{bait_localisation}.csv",
            bait_localisation = config["localisation"]
        )
    output:
        biotid_all = "work_folder/localisation_probability/bait_all.csv"
    run:
        with open(output.biotid_all, "w") as w:
            w.write("localisation_bait\tlocalisation_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tin_permutation\tin_bioid\n")
            for localisation_file in input.all_baits_probs:
                lines = [line for line in open(localisation_file, "r")]
                for line in lines[1:]:
                    w.write(line)

        shell("awk '(NR == 1) || (FNR > 1)' {input.all_baits_probs} > combined.csv")





