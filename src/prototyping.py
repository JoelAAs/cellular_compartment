##
#import dask.dataframe as dd
import time

import pandas as pd
from dask import delayed
from multiprocessing import Pool
### Functions

def get_prey_likelihood(bait_prey_df):
    bait_prey_counts = bait_prey_df.groupby(["target_desc_bait", "target_desc_prey"], as_index=False).size()
    bait_sum = bait_prey_counts.groupby("target_desc_bait", as_index=False)["size"].sum()
    bait_prey_counts = bait_prey_counts.merge(bait_sum, on="target_desc_bait", suffixes=("_bait_prey_count", "_bait_count"))
    bait_prey_counts["likelihood_prey"] = bait_prey_counts["size_bait_prey_count"].div(bait_prey_counts["size_bait_count"])
    return bait_prey_counts


def worker(args):
    i, frac, other_ms_ss, get_prey_likelihood = args
    print(f"running set {i}")
    permutation_df = other_ms_ss.sample(frac=frac)
    permutation_prey_counts = get_prey_likelihood(permutation_df)
    permutation_prey_counts["permutation"] = i
    permutation_prey_counts.to_csv(f"work_folder/localization_permutations/set_frac_{frac}_set_{i}.csv", sep = "\t", index=False)
    return i

### Input config
pwd = "~/Projects/PPI-bias"
intact_df = pd.read_csv("data/bait_prey_publications.csv", sep = "\t")
localisation_df = pd.read_csv("data/gene_attribute_edges.txt", sep = "\t")
localisation_df = localisation_df[["source", "target_desc"]]
gene_name_to_uniprot = pd.read_csv("data/uniprot_to_gene_name.csv", sep = "\t")

intact_df = intact_df.merge(
    gene_name_to_uniprot,
    left_on="bait", right_on= "uniprot_id"
)

intact_df = intact_df.merge(
    gene_name_to_uniprot,
    left_on="prey", right_on= "uniprot_id",
    suffixes=("_bait", "_prey")
)
intact_df = intact_df.merge(
    localisation_df,
    left_on="gene_name_bait", right_on= "source",
)
intact_df = intact_df.merge(
    localisation_df,
    left_on="gene_name_prey", right_on="source",
    suffixes=("_bait", "_prey")
)
intact_df = intact_df[[
    "gene_name_bait",
    "gene_name_prey",
    "detection_method",
    "target_desc_bait",
    "target_desc_prey"]
]
other_ms_methods = [
    "MI-0006",
    "MI-0007",
    "MI-0096",
    "MI-0004",
    "MI-0019"
]


### Permutations
bioID_ss = intact_df[intact_df["detection_method"] == "MI-1314"]
other_ms_ss = intact_df[intact_df["detection_method"].isin(other_ms_methods)]
print("other_ms_ss computed")
n_permutations = 50000
n_cores = 6
frac = .9

time_start = time.time()
with Pool(n_cores) as pool:
    results = pool.map(
        worker,
        [(i, frac, other_ms_ss, get_prey_likelihood) for i in range(n_permutations)]
    )
print(f"{time.time() - time_start} seconds spent on {n_permutations} using {n_cores} cores.")
###
bioid_counts = get_prey_likelihood(bioID_ss)
bioid_counts.to_csv("work_folder/bioID_localisation.csv", sep="\t", index=False)


