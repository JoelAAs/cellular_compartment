import pandas as pd
import numpy as np
from scipy.stats import norm
from support_functions import read_ppi
n_permutations = config["n_permutations"]


def get_prey_likelihood(bait_prey_df):
    bait_prey_counts = bait_prey_df.groupby(["target_desc_bait", "gene_name_prey"], as_index=False).size()
    bait_sum = bait_prey_counts.groupby(["target_desc_bait"], as_index=False)["size"].sum()
    bait_prey_counts = bait_prey_counts.merge(bait_sum, on="target_desc_bait", suffixes=("_prey_count", "_bait_localisation_count"))
    bait_prey_counts["likelihood_prey"] = bait_prey_counts["size_prey_count"].div(bait_prey_counts["size_bait_localisation_count"])
    return bait_prey_counts


rule bioid_per_prey_localisation:
    input:
        intact = "data/bait_prey_publications.csv",
        localisation_annotations = "data/gene_attribute_edges.txt",
        uniprot_gene_name = "data/uniprot_to_gene_name.csv",
        common_baits = "work_folder/enrichment_analysis/bait_list_common.csv"
    output:
        bioid_file = "work_folder/prey_probability/bioID_localisation.csv"
    run:
        intact_df = read_ppi(input.intact, input.localisation_annotations, input.uniprot_gene_name)
        common_baits = pd.read_csv(input.common_baits, sep="\t")

        bioid_ms_ss = intact_df[intact_df["detection_method"] == "MI-1314"]
        bioid_ms_ss = bioid_ms_ss[bioid_ms_ss["gene_name_bait"].isin(common_baits["gene_name"])]
        permutation_prey_counts = get_prey_likelihood(bioid_ms_ss)
        permutation_prey_counts.to_csv(
            output.bioid_file,
            sep = "\t", index = False
        )

rule permute_localisation_to_prey:
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
        intact = "data/bait_prey_publications.csv",
        localisation_annotations = "data/gene_attribute_edges.txt",
        uniprot_gene_name = "data/uniprot_to_gene_name.csv",
        common_baits= "work_folder/enrichment_analysis/bait_list_common.csv"
    output:
        permutation_sets = expand(
            "work_folder/prey_probability/permutations/{{localisation}}_0.9_set_{n}.csv",
            n = range(n_permutations)
        )
    run:
        intact_df = read_ppi(input.intact, input.localisation_annotations, input.uniprot_gene_name)
        common_baits = pd.read_csv(input.common_baits,sep="\t")

        ms_df = intact_df.loc[
            intact_df["detection_method"].isin(params.other_ms_methods)
        ]
        ms_df = ms_df[ms_df["gene_name_bait"].isin(common_baits["gene_name"])]

        ms_localisation_df = ms_df[ms_df["target_desc_bait"] == wildcards.localisation]
        for i, permutation_file in enumerate(output.permutation_sets):
            sample_df = ms_localisation_df.sample(frac=params.frac)
            permutation_prey_counts = get_prey_likelihood(sample_df)
            permutation_prey_counts["permutation"] = i
            permutation_prey_counts.to_csv(
                permutation_file,
                sep = "\t", index = False
            )


rule estimate_quant_prey:
    input:
        permutation_sets = expand(
            "work_folder/prey_probability/permutations/{{localisation}}_0.9_set_{n}.csv",
            n = range(n_permutations)
        ),
        biotin_file= "work_folder/prey_probability/bioID_localisation.csv"
    output:
        biotid_cdf_quant = "work_folder/prey_probability/bioid_prey_quantile/bait_{localisation}.csv"
    run:
        localisation_df_list = [
            pd.read_csv(
                file,
                sep="\t"
            ) for file in input.permutation_sets
        ]
        localisation_df = pd.concat(localisation_df_list)

        biotin_df = pd.read_csv(
            input.biotin_file,
            sep="\t"
        )
        with open(output.biotid_cdf_quant, "w") as w:
            w.write("target_desc_bait\tgene_name_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tfrac_in_permutation\tin_bioid\n")
            for current_gene_prey in localisation_df["gene_name_prey"].unique():
                biotin_observed = True
                all_probabilities = np.zeros(n_permutations)
                ss_df = localisation_df[localisation_df["gene_name_prey"] == current_gene_prey]
                permuted = np.array(ss_df["likelihood_prey"])

                frac_observed = len(permuted)/n_permutations
                if len(permuted) != 0:
                    all_probabilities[:len(permuted)] = permuted

                observed_value = biotin_df.loc[
                                 (biotin_df["gene_name_prey"] == current_gene_prey) &
                                 (biotin_df["target_desc_bait"] == wildcards.localisation)
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
                    f"{current_gene_prey}\t"
                    f"{mu}\t"
                    f"{std}\t"
                    f"{quant_value}\t"
                    f"{observed_value}\t"
                    f"{frac_observed}\t"
                    f"{biotin_observed}\n"
                )

rule aggregate_prey_data:
    input:
        all_baits_probs = expand(
            "work_folder/prey_probability/bioid_prey_quantile/bait_{bait_localisation}.csv",
            bait_localisation = config["localisation"]
        )
    output:
        biotid_all = "work_folder/prey_probability/bait_all.csv"
    run:
        with open(output.biotid_all, "w") as w:
            w.write("target_desc_bait\tgene_name_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tin_permutation\tin_bioid\n")
            for localisation_file in input.all_baits_probs:
                lines = [line for line in open(localisation_file, "r")]
                for line in lines[1:]:
                    w.write(line)

        shell("awk '(NR == 1) || (FNR > 1)' {input.all_baits_probs} > combined.csv")
