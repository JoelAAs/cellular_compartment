import pandas as pd
import numpy as np
from scipy.stats import norm

n_permutations = 100000
localisation_classification = [
    "GO:0000137", "GO:0005856", "GO:0042470",
    "GO:0000138", "GO:0005886", "GO:0042579",
    "GO:0000323", "GO:0005911", "GO:0043226",
    "GO:0005575", "GO:0005923", "GO:0043227",
    "GO:0005576", "GO:0008021", "GO:0043228",
    "GO:0005634", "GO:0016020", "GO:0043229",
    "GO:0005635", "GO:0016023", "GO:0043231",
    "GO:0005730", "GO:0016323", "GO:0043232",
    "GO:0005737", "GO:0016324", "GO:0044422",
    "GO:0005739", "GO:0016604", "GO:0044424",
    "GO:0005741", "GO:0016607", "GO:0044425",
    "GO:0005743", "GO:0019866", "GO:0044428",
    "GO:0005764", "GO:0019867", "GO:0044429",
    "GO:0005768", "GO:0030054", "GO:0044430",
    "GO:0005769", "GO:0030133", "GO:0044431",
    "GO:0005770", "GO:0030141", "GO:0044444",
    "GO:0005773", "GO:0031090", "GO:0044446",
    "GO:0005777", "GO:0031410", "GO:0044451",
    "GO:0005783", "GO:0031966", "GO:0044456",
    "GO:0005793", "GO:0031967", "GO:0044459",
    "GO:0005794", "GO:0031968", "GO:0044464",
    "GO:0005797", "GO:0031975", "GO:0048770",
    "GO:0005802", "GO:0031982", "GO:0070160",
    "GO:0005811", "GO:0031984", "GO:0098588",
    "GO:0005813", "GO:0031985", "GO:0098589",
    "GO:0005815", "GO:0031988", "GO:0098590"

]

def bin_it(size, n_permutations, start, batch_list):
    if start + size >= n_permutations:
        batch_list.append(range(start, n_permutations))
        return batch_list
    else:
        end = start + size
        batch_list.append(range(start, end))
        return bin_it(size, n_permutations, end, batch_list)

batches = bin_it(5000, n_permutations, 0, [])
n_batches = len(batches)
config["batches"] = batches

def get_localisation_set(wc):
    expected_input = [
        f"work_folder/localization_permutations/set_frac_0.9_set_{i}.csv"
        for i in config["batches"][int(wc.n)]
    ]
    return expected_input



rule get_permutations_localisation:
    input:
        permut_sets = lambda wc: get_localisation_set(wc)
    output:
        expand(
            "work_folder/localization_permutations/per_localisation/{localisation}_permutations_0.9_batch_{{n}}.csv",
            localisation=localisation_classification
        )
    run:
        file_dict = dict()
        for permutation_file in input.permut_sets:
            with open(permutation_file, "r") as f:
                header = True
                for line in f:
                    if header:
                        header = False
                    else:
                       row_values = line.strip().split("\t")
                       try:
                           file_dict[row_values[0]].write(line)
                       except KeyError:
                           output_name = f"work_folder/localization_permutations/per_localisation/{row_values[0]}_permutations_0.9_batch_{wildcards.n}.csv"
                           file_dict[row_values[0]] = open(output_name, "w")
                           file_dict[row_values[0]].write("target_desc_bait\ttarget_desc_prey\tsize_bait_prey_count\tsize_bait_count\tlikelihood_prey\tpermutation\n")
                           file_dict[row_values[0]].write(line)
        for key in file_dict:
            file_dict[key].close()


rule estimate_quant:
    params:
        n_permuts = n_permutations
    input:
        biotin_file = "work_folder/bioID_localisation.csv",
        localisation_bait_probability = expand(
            "work_folder/localization_permutations/per_localisation/{bait_localisation}_permutations_0.9_batch_{n}.csv",
            n = range(n_batches)
        )
    output:
        biotid_cdf_quant = "work_folder/bioid_quantile/bait_{bait_localisation}.csv"
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
            w.write("target_desc_bait\ttarget_desc_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tin_permutation\tin_bioid\n")
            for current_localisation in localisation_classification:
                observed = True
                biotin_observed = True
                all_probabilities = np.zeros(n_permutations)
                ss_df = localisation_df[localisation_df["target_desc_prey"] == current_localisation]
                permuted = np.array(ss_df["likelihood_prey"])
                if len(permuted) == 0:
                    observed = False
                all_probabilities[:len(permuted)] = permuted

                observed_value = biotin_df.loc[
                    (biotin_df["target_desc_prey"] == current_localisation) &
                    (biotin_df["target_desc_bait"] == wildcards.bait_localisation)]["likelihood_prey"].values
                if len(observed_value) > 0:
                    observed_value=observed_value[0]
                else:
                    observed_value = 0
                    biotin_observed = False


                mu = np.mean(all_probabilities)
                std = np.std(all_probabilities)
                quant_value = norm.cdf(observed_value, mu, std)
                w.write(
                    f"{wildcards.bait_localisation}\t"
                    f"{current_localisation}\t"
                    f"{mu}\t"
                    f"{std}\t"
                    f"{quant_value}\t"
                    f"{observed_value}\t"
                    f"{observed}\t"
                    f"{biotin_observed}\n"
                )


rule aggregate_quant_data:
    input:
        expand(
            "work_folder/bioid_quantile/bait_{bait_localisation}.csv",
            bait_localisation = localisation_classification)
    output:
        biotid_all = "work_folder/bioid_quantile/bait_all.csv"
    run:
        with open(output.biotid_all, "w") as w:
            w.write("target_desc_bait\ttarget_desc_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tin_permutation\tin_bioid\n")

        shell("awk '(NR == 1) || (FNR > 1)' {input.bait_localisation} > combined.csv")





