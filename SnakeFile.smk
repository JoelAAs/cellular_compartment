#config["n_permutations"] = 1000

config["localisations_considered"] = [
    'Acrosome', 'Actin-filaments', 'Aggresome','Annulus',
    'Basal-body', 'Calyx', 'Cell-Junctions', 'Centriolar-satellite',
    'Centrosome', 'Cleavage-furrow', 'Connecting-piece', 'Cytokinetic-bridge',
    'Cytoplasmic-bodies', 'Cytosol', 'End-piece', 'Endoplasmic-reticulum',
    'Endosomes', 'Equatorial-segment', 'Flagellar-centriole', 'Focal-adhesion-sites',
    'Golgi-apparatus', 'Intermediate-filaments', 'Kinetochore', 'Lipid-droplets',
    'Lysosomes', 'Microtubule-ends', 'Microtubules', 'Mid-piece',
    'Midbody', 'Midbody-ring', 'Mitochondria', 'Mitotic-chromosome',
    'Mitotic-spindle', 'Nuclear-bodies', 'Nuclear-membrane', 'Nuclear-speckles',
    'Nucleoli', 'Nucleoli-fibrillar-center', 'Nucleoli-rim',
    'Nucleoplasm', 'Perinuclear-theca', 'Peroxisomes', 'Plasma-membrane',
    'Primary-cilium', 'Primary-cilium-tip', 'Primary-cilium-transition-zone', 'Principal-piece',
    'Rods-and-Rings', 'Vesicles'
]

config["methods_considered"] = [
    "MI-0006",
    "MI-0007",
    #"MI-0018",
    "MI-0096",
    #"MI-0397",
    #"MI-0398",
    "MI-0676",
    #"MI-1112",
    "MI-1314"
]

wildcard_constraints:
    method="MI-[0-9]+"

# include: "src/get_biotin_network.smk"
# include: "src/get_prey_localisation.smk"
# include: "src/misc_analysis.smk"
# include: "src/enrichment_analysis.smk"
# include: "src/bait_distribution.smk"
include: "src/method_localisation/method_localisation.smk"
include: "src/method_localisation/localisation_estimation.smk"

rule all:
    input:
        "work_folder/localisation/beta_estimation.csv",
        "work_folder/method_localisation/beta_estimation/method_localisation_p_estimation.csv",
        "work_folder/localisation_overlap.csv",
        "work_folder/localisation/fisher_exact.csv"
        # "work_folder/enrichment_analysis/plots/venn_diagram_doid.png",
        # "work_folder/localisation_probability/bait_all.csv",
        # "work_folder/prey_probability/bait_all.csv",
        #"work_folder/localisation_data/iou.csv",
        #"work_folder/localisation_data/shared.csv"
