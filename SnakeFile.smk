
include: "src/get_biotin_network.smk"


rule all:
    input:
        "work_folder/bioid_quantile/bait_all.csv"