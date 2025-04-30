import pandas as pd
from support_functions import read_ppi
rule get_selected_subset:
    input:
        intact = "data/bait_prey_publications.csv",
        localisation_annotations = "data/localisation/gene_to_location.csv",
        uniprot_gene_name = "data/uniprot_to_gene_name.csv",
        common_baits = "work_folder_{project}/enrichment_analysis/bait_list_common.csv"
    output:
        selected_df = "work_folder_{project}/selected_baits.csv"
    run:
        intact_df = read_ppi(input.intact, input.localisation_annotations, input.uniprot_gene_name)
        common_baits = pd.read_csv(input.common_baits,sep="\t")
        intact_df = intact_df[intact_df["gene_name_bait"].isin(common_baits["gene_name"])]
        intact_df.to_csv(output.selected_df, sep="\t", index=False)




rule get_shared_and_iou_localisation:
    input:
        local_data = "data/gene_attribute_edges.txt"
    output:
        iou = "work_folder_{project}/localisation_data/iou.csv",
        shared = "work_folder_{project}/localisation_data/shared.csv"
    run:
        localisation_data = pd.read_csv(input.local_data, sep="\t")
        localisation_data = localisation_data.iloc[1:] # weird second header
        local_annotation = localisation_data["target_desc"].unique()
        with open(output.iou, "w") as w_iou:
            with open(output.shared, "w") as w_shared:
                w_shared.write("localisation_a\tlocalisation_b\tn_shared\n")
                w_iou.write("localisation_a\tlocalisation_b\tiou\n")
                for localisation_GO_from in local_annotation:
                    for localisation_GO_to in local_annotation:
                        genes_from = set(localisation_data[
                            localisation_data["target_desc"] == localisation_GO_from]["source"].tolist())
                        genes_to = set(localisation_data[
                            localisation_data["target_desc"] == localisation_GO_to]["source"].tolist())

                        union_set = genes_from | genes_to
                        intersect_set = genes_from & genes_to

                        w_shared.write(f"{localisation_GO_from}\t{localisation_GO_to}\t{len(intersect_set)}\n")
                        w_iou.write(f"{localisation_GO_from}\t{localisation_GO_to}\t{len(intersect_set)/len(union_set)}\n")



