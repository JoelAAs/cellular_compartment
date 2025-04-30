

wildcard_constraints:
    method="MI-[0-9]+"

# include: "src/get_biotin_network.smk"
# include: "src/get_prey_localisation.smk"
# include: "src/misc_analysis.smk"
# include: "src/enrichment_analysis.smk"
# include: "src/bait_distribution.smk"
# include: "src/method_localisation/method_localisation.smk"
include: "src/method_localisation/localisation_estimation.smk"
#include: "src/check_mitochondria/infer_negative_data.smk"
include: "src/estimate_pair_probability/Infer_negative_data.smk"


rule all:
    input:
        "work_folder/inferred_negative_localisation.csv",
        f"work_folder_{config['project']}/localisation/beta_estimation.csv",
        # expand(
        #     "work_folder_{project}/method_localisation/{cell_line}/beta_estimation/method_localisation_p_estimation.csv",
        #     cell_line=config["cell_lines"], project=config["project"]
        # ),
        # f"work_folder_{config['project']}/localisation_overlap.csv",
        # f"work_folder_{config['project']}/localisation/fisher_exact.csv",
        # "work_folder/Mitochondria/interactions_observed_vs_theoretical.csv",
        # "work_folder/Endoplasmic-reticulum/binomial_ppi_observed.csv",
        # "work_folder/Nucleoplasm/binomial_ppi_observed.csv",
        # "work_folder/Mitochondria/binomial_ppi_observed.csv"
        # "work_folder/enrichment_analysis/plots/venn_diagram_doid.png"
        # "work_folder/localisation_probability/bait_all.csv",
        # "work_folder/prey_probability/bait_all.csv",
        #"work_folder/localisation_data/iou.csv",
        #"work_folder/localisation_data/shared.csv"
