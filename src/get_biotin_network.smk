import pandas as pd
import numpy as np
from scipy.stats import norm

n_permutations = 10000
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

def get_prey_likelihood(bait_prey_df):
    bait_prey_counts = bait_prey_df.groupby(["target_desc_bait", "target_desc_prey"], as_index=False).size()
    bait_sum = bait_prey_counts.groupby("target_desc_bait", as_index=False)["size"].sum()
    bait_prey_counts = bait_prey_counts.merge(bait_sum, on="target_desc_bait", suffixes=("_bait_prey_count", "_bait_count"))
    bait_prey_counts["likelihood_prey"] = bait_prey_counts["size_bait_prey_count"].div(bait_prey_counts["size_bait_count"])
    return bait_prey_counts

def read_ppi(intact, localisation_annotations, uniprot_gene_name):
    intact_df = pd.read_csv(intact,sep="\t")
    localisation_df = pd.read_csv(localisation_annotations,sep="\t")
    localisation_df = localisation_df[["source", "target_desc"]]
    gene_name_to_uniprot = pd.read_csv(uniprot_gene_name,sep="\t")

    intact_df = intact_df.merge(
        gene_name_to_uniprot,
        left_on="bait",right_on="uniprot_id"
    )

    intact_df = intact_df.merge(
        gene_name_to_uniprot,
        left_on="prey",right_on="uniprot_id",
        suffixes=("_bait", "_prey")
    )
    intact_df = intact_df.merge(
        localisation_df,
        left_on="gene_name_bait",right_on="source",
    )
    intact_df = intact_df.merge(
        localisation_df,
        left_on="gene_name_prey",right_on="source",
        suffixes=("_bait", "_prey")
    )
    intact_df = intact_df[[
        "gene_name_bait",
        "gene_name_prey",
        "detection_method",
        "target_desc_bait",
        "target_desc_prey"]
    ]
    return intact_df


rule bioid_per_bait_localisation:
    input:
        intact = "data/bait_prey_publications.csv",
        localisation_annotations = "data/gene_attribute_edges.txt",
        uniprot_gene_name = "data/uniprot_to_gene_name.csv"
    output:
        bioid_file = "work_folder/bioID_localisation.csv"
    run:
        intact_df = read_ppi(input.intact, input.localisation_annotations, input.uniprot_gene_name)

        bioid_ms_ss = intact_df[intact_df["detection_method"] == "MI-1314"]
        permutation_prey_counts = get_prey_likelihood(bioid_ms_ss)
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
        intact = "data/bait_prey_publications.csv",
        localisation_annotations = "data/gene_attribute_edges.txt",
        uniprot_gene_name = "data/uniprot_to_gene_name.csv"
    output:
        permutation_sets = expand(
            "work_folder/bait_localisation_permutation/{{localisation}}_0.9_set_{n}.csv",
            n = range(n_permutations)
        )
    run:
        intact_df = read_ppi(input.intact, input.localisation_annotations, input.uniprot_gene_name)

        other_ms_ss = intact_df[intact_df["detection_method"].isin(params.other_ms_methods)]
        other_ms_localisation_ss = other_ms_ss[other_ms_ss["target_desc_bait"] == wildcards.localisation]
        for i, permutation_file in enumerate(output.permutation_sets):
            sample_df = other_ms_localisation_ss.sample(frac=params.frac)
            permutation_prey_counts = get_prey_likelihood(sample_df)
            permutation_prey_counts["permutation"] = i
            permutation_prey_counts.to_csv(
                permutation_file,
                sep = "\t", index = False
            )


rule estimate_quant:
    params:
        n_permuts = n_permutations
    input:
        biotin_file = "work_folder/bioID_localisation.csv",
        localisation_bait_probability = expand(
            "work_folder/bait_localisation_permutation/{{localisation}}_0.9_set_{n}.csv",
            n = range(n_batches)
        )
    output:
        biotid_cdf_quant = "work_folder/bioid_quantile/bait_{localisation}.csv"
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
            w.write("target_desc_bait\ttarget_desc_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tfrac_in_permutation\tin_bioid\n")
            for current_localisation in localisation_classification:
                biotin_observed = True
                all_probabilities = np.zeros(n_permutations)
                ss_df = localisation_df[localisation_df["target_desc_prey"] == current_localisation]
                permuted = np.array(ss_df["likelihood_prey"])

                frac_observed = len(permuted)/n_permutations
                if len(permuted) != 0:
                    all_probabilities[:len(permuted)] = permuted

                observed_value = biotin_df.loc[
                    biotin_df["target_desc_prey"] == current_localisation
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
            "work_folder/bioid_quantile/bait_{bait_localisation}.csv",
            bait_localisation = localisation_classification
        )
    output:
        biotid_all = "work_folder/bioid_quantile/bait_all.csv"
    run:
        with open(output.biotid_all, "w") as w:
            w.write("target_desc_bait\ttarget_desc_prey\tprobability_mean\tprobability_std\tquantile_value\tobserved_value\tin_permutation\tin_bioid\n")
            for localisation_file in input.all_baits_probs:
                lines = [line for line in open(localisation_file, "r")]
                for line in lines[1:]:
                    w.write(line)

        shell("awk '(NR == 1) || (FNR > 1)' {input.all_baits_probs} > combined.csv")





