import pandas as pd
from scipy.stats import fisher_exact

rule get_shared_distribution:
    params:
        limit_pval = 0.05
    input:
        bioid_baits="work_folder/enrichment_analysis/bait_lists/bioid_baits.csv",
        ms_baits="work_folder/enrichment_analysis/bait_lists/ms_baits.csv"
    output:
        shared = "work_folder/enrichment_analysis/bait_list_common.csv"
    run:
        bioid_bait_df = pd.read_csv(
            input.bioid_baits
        )
        bioid_bait_df = bioid_bait_df.groupby("gene_name", as_index=False).size()
        bioid_bait_df["n"] = bioid_bait_df["size"].sum()
        bioid_bait_df["p"] = bioid_bait_df["size"].div(bioid_bait_df["n"])

        ms_bait_df = pd.read_csv(
            input.ms_baits
        )
        ms_bait_df = ms_bait_df.groupby("gene_name",as_index=False).size()
        ms_bait_df["n"] = ms_bait_df["size"].sum()
        ms_bait_df["p"] = ms_bait_df["size"].div(ms_bait_df["n"])
        baits_df = bioid_bait_df.merge(ms_bait_df, on="gene_name", suffixes=["_bioid", "_ms"])

        def fisher_test(row):
            table = [
                [row["size_bioid"], row["n_bioid"] - row["size_bioid"]],
                [row["size_ms"], row["n_bioid"] - row["size_ms"]]
            ]
            res = fisher_exact(table)
            return res.pvalue

        baits_df["p"] = baits_df.apply(fisher_test, axis = 1)
        common_baits = baits_df[baits_df["p"] > params.limit_pval]
        common_baits.to_csv(output.shared, sep ="\t", index=False)

