config["n_permutations"] = 1000
# config["localisation"] = [
#     "GO:0000137", "GO:0005856", "GO:0042470",
#     "GO:0000138", "GO:0005886", "GO:0042579",
#     "GO:0000323", "GO:0005911", "GO:0043226",
#     "GO:0005575", "GO:0005923", "GO:0043227",
#     "GO:0005576", "GO:0008021", "GO:0043228",
#     "GO:0005634", "GO:0016020", "GO:0043229",
#     "GO:0005635", "GO:0016023", "GO:0043231",
#     "GO:0005730", "GO:0016323", "GO:0043232",
#     "GO:0005737", "GO:0016324", "GO:0044422",
#     "GO:0005739", "GO:0016604", "GO:0044424",
#     "GO:0005741", "GO:0016607", "GO:0044425",
#     "GO:0005743", "GO:0019866", "GO:0044428",
#     "GO:0005764", "GO:0019867", "GO:0044429",
#     "GO:0005768", "GO:0030054", "GO:0044430",
#     "GO:0005769", "GO:0030133", "GO:0044431",
#     "GO:0005770", "GO:0030141", "GO:0044444",
#     "GO:0005773", "GO:0031090", "GO:0044446",
#     "GO:0005777", "GO:0031410", "GO:0044451",
#     "GO:0005783", "GO:0031966", "GO:0044456",
#     "GO:0005793", "GO:0031967", "GO:0044459",
#     "GO:0005794", "GO:0031968", "GO:0044464",
#     "GO:0005797", "GO:0031975", "GO:0048770",
#     "GO:0005802", "GO:0031982", "GO:0070160",
#     "GO:0005811", "GO:0031984", "GO:0098588",
#     "GO:0005813", "GO:0031985", "GO:0098589",
#     "GO:0005815", "GO:0031988", "GO:0098590"
# ]

config["localisation"] = [
    'Acrosome', 'Actin_filaments', 'Aggresome','Annulus',
    'Basal_body', 'Calyx', 'Cell_Junctions', 'Centriolar_satellite',
    'Centrosome', 'Cleavage_furrow', 'Connecting_piece', 'Cytokinetic_bridge',
    'Cytoplasmic_bodies', 'Cytosol', 'End_piece', 'Endoplasmic_reticulum',
    'Endosomes', 'Equatorial_segment', 'Flagellar_centriole', 'Focal_adhesion_sites',
    'Golgi_apparatus', 'Intermediate_filaments', 'Kinetochore', 'Lipid_droplets',
    'Lysosomes', 'Microtubule_ends', 'Microtubules', 'Mid_piece',
    'Midbody', 'Midbody_ring', 'Mitochondria', 'Mitotic_chromosome',
    'Mitotic_spindle', 'Nuclear_bodies', 'Nuclear_membrane', 'Nuclear_speckles',
    'Nucleoli', 'Nucleoli_fibrillar_center', 'Nucleoli_rim',
    'Nucleoplasm', 'Perinuclear_theca', 'Peroxisomes', 'Plasma_membrane',
    'Primary_cilium', 'Primary_cilium_tip', 'Primary_cilium_transition_zone', 'Principal_piece',
    'Rods_and_Rings', 'Vesicles'
]

include: "src/get_biotin_network.smk"
include: "src/get_prey_localisation.smk"
include: "src/misc_analysis.smk"
include: "src/enrichment_analysis.smk"
include: "src/bait_distribution.smk"

rule all:
    input:
        #"work_folder/enrichment_analysis/plots/venn_diagram_doid.png",
        "work_folder/localisation_probability/bait_all.csv",
        "work_folder/prey_probability/bait_all.csv"
        # "work_folder/localisation_data/iou.csv",
        # "work_folder/localisation_data/shared.csv"
