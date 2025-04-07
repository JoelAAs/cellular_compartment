source("src/enrichment_support_functions.R")
library(tidyverse)
library(ggVennDiagram)
library(magrittr)


### Functions
enrichment_per_bait_group <- function(df){
  i = 1
  go_list <- list()
  do_list <- list()
  for (bait_go in unique(df$target_desc_bait)) {
    bait_go_df <- df |> 
      filter(target_desc_bait == bait_go)
    go_result <- enrichment_GO_results(bait_go_df$entrez_id)
    go_result$bait_go   <- bait_go
    go_result$bait_go_N <- nrow(bait_go_df) 
    do_result <- enrichment_DO_results(bait_go_df$entrez_id)
    do_result$bait_go   <- bait_go
    do_result$bait_go_N <- nrow(bait_go_df)
    go_list[[i]] = go_result
    do_list[[i]] = do_result
      i = i+1
  }
  return(list(
    go_results = go_list, 
    do_results = do_list)
    )
}

# Run
prey_probabilies <- read.table("work_folder/prey_probability/bait_all.csv", header=T)
prey_probabilies |>
  filter(in_bioid=="True") |>
  filter(in_permutation > 0) -> prey_probabilies

prey_probabilies <- map_entrez(prey_probabilies, "gene_name_prey")
prey_probabilies$z <- (prey_probabilies$observed_value - prey_probabilies$probability_mean)/prey_probabilies$probability_std

ggplot(prey_probabilies,
       aes(x=z)) +
  geom_histogram() +
  facet_wrap(target_desc_bait ~.)

prey_probabilies |>
  top_frac(0.05, z) -> top_5

prey_probabilies |>
  top_frac(0.05, -z) -> bot_5

top_enrichments <- enrichment_per_bait_group(top_5)
bot_enrichments = enrichment_per_bait_group(bot_5)

top_go_df <- bind_rows(top_enrichments$go_results)
top_do_df <- bind_rows(top_enrichments$do_results)

bot_go_df <- bind_rows(bot_enrichments$go_results)
bot_do_df <- bind_rows(bot_enrichments$do_results)

top_go_df %>% 
  filter(qvalue < 0.05) %>%  
  {table(.$ID)} %>% 
  stack() -> significant_enrichment_top
colnames(significant_enrichment_top) <- c(
  "top_count",
  "go_term"
)

bot_go_df %>% 
  filter(qvalue < 0.05) %>%  
  {table(.$ID)} %>% 
  stack() -> significant_enrichment_bot
colnames(significant_enrichment_bot) <- c(
  "bot_count",
  "go_term"
)

merge(significant_enrichment_bot, significant_enrichment_top, by="go_term", how=)

