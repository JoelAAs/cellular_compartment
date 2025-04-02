library(EnrichDO)
library(org.Hs.eg.db)
library(tidyverse)
library(ggVennDiagram)

#### functions

get_top_rows <- function(gene_name_df, n_top=100){
  entrez_id_df <- map_entrez(gene_name_df)
  df_baits_entrez <- merge(gene_name_df, entrez_id_df, on="gene_name",how="inner")
  df_baits_entrez |> 
    group_by(gene_name,  entrez_id) |>
    summarise(n = n()) |>
    ungroup() -> df_bait_freq_df
  
  df_bait_freq_df$freq <- df_bait_freq_df$n/sum(df_bait_freq_df$n)
  df_bait_freq_df_top <- df_bait_freq_df |> top_n(n_top, freq)
  return(df_bait_freq_df_top)
}

enrichment_results <- function(df){
  enrichTest <- EnrichDO::doEnrich(df$entrez_id, method = "fdr")
  return(enrichTest@enrich)
}

#### args
args <- commandArgs(trailingOnly = TRUE)

bioid_bait_input              <- args[1]
ms_bait_input                <- args[2]
bioid_bait_enrichment_output <- args[3]
ms_bait_enrichment_output    <- args[4]
venn_plot_bait               <- args[5]
venn_plot_doid               <- args[6]
bait_count                   <- as.numeric(args[7])

### run
df_bioid_baits <- read.table(bioid_bait_input, header=T)
df_ms_baits <- read.table(ms_bait_input, header=T)

df_bioid_baits_top <- get_top_rows(df_bioid_baits, bait_count)
results_bioid <- enrichment_results(df_bioid_baits_top)
results_bioid |> write.table(bioid_bait_enrichment_output, sep ="\t")

df_ms_baits_top <- get_top_rows(df_ms_baits, bait_count)
results_ms <- enrichment_results(df_ms_baits_top)
results_ms |> write.table(ms_bait_enrichment_output, sep ="\t")

venn_doid <- list(Bioid = results_bioid$DOID, MS = results_ms$DOID)
venn_baits <- list(Bioid = df_bioid_baits_top$gene_name, MS = df_ms_baits_top$gene_name)

doid_vd <- ggVennDiagram(venn_doid) +
  ggtitle("Number of shared significantly enriched DOIDs")
ggsave(venn_plot_doid, doid_vd, dpi=300)
bait_vd <- ggVennDiagram(venn_baits) +
  ggtitle("Number of shared significantly enriched DOIDs")
ggsave(venn_plot_bait, bait_vd, dpi=300)

