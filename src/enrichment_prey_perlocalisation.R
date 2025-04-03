source("src/enrichment_support_functions.R")
library(tidyverse)
library(ggVennDiagram)


### Functions


enrichment_per_bait_group <- function(df){
  i = 1
  go_list <- list()
  do_list <- list()
  for (bait_go in unique(df$target_desc_bait)) {
    bait_go_df <- df |> 
      filter(target_desc_bait == bait_go)
    go_result <- enrichment_GO_results(bait_go_df)
    go_result$bait_go   <- bait_go
    go_result$bait_go_N <- nrow(bait_go_df) 
    do_result <- enrichment_DO_results(bait_go_df)
    do_result$bait_go   <- bait_go
    do_result$bait_go_N <- nrow(bait_go_df)
    go_list <- list()
    do_list <- list()
    
    go_list[[i]] = go_result
    do_list[[i]] = do_result
      i = i+1
  }
  return(c(go_list, do_list))
}

# Run
prey_probabilies <- read.table("work_folder/prey_probability/bait_all.csv", header=T)
prey_probabilies |>
  filter(in_bioid=="True") |>
  filter(in_permutation > 0) -> prey_probabilies

prey_probabilies <- map_entrez(prey_probabilies)
prey_probabilies$z <- (prey_probabilies$observed_value - prey_probabilies$probability_mean)/prey_probabilies$probability_std

prey_probabilies |>
  top_n(300, z) -> top_300

prey_probabilies |>
  top_n(3000, -z) -> bot_300

a = enrichment_per_bait_group(bot_300)


test <- enrichment_per_bait_group(top_300)

top_enrich <- enrichment_results(top_300)
bot_enrich <- enrichment_results(bot_300)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]


