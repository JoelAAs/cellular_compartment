library(tidyverse)
library(ggplot2)

permut_prey_probability = read.table(
  "combined.csv",
  sep = "\t",
  header = T
)

bioID_prey_probability = read.table(
  "work_folder/bioID_localisation.csv",
  sep="\t",
  header=T
)
bioID_prey_probability$target_desc_bait <- factor(
  bioID_prey_probability$target_desc_bait,
  levels = unique(c(bioID_prey_probability$target_desc_bait, bioID_prey_probability$target_desc_prey)))
bioID_prey_probability$target_desc_prey <- factor(
  bioID_prey_probability$target_desc_prey,
  levels = unique(c(bioID_prey_probability$target_desc_bait, bioID_prey_probability$target_desc_prey)))


permut_prey_probability |>
  filter(target_desc_bait == target_desc_prey) -> same_localisation_df

get_quant_value <- function(permut_prey_probability_df, observed_value_dfs, n_permuts){
  localisations = unique(permut_prey_probability_df$target_desc_bait)
  expected_computations = length(localisations)^2
  quant_value_df = data.frame(
    bait_localication = rep(character(0), expected_computations),
    prey_localication = rep(character(0), expected_computations),
    quant_value = rep(numeric(0), expected_computations),
    ms_data = rep(logical(0), expected_computations),
    biotiny = rep(logical(0), expected_computations)
  )
  i=0
  for (loc_bait in localisations) {
    for (loc_prey in localisations) {
      probs = rep(0, n_permuts)
      permuts_df = permut_prey_probability_df |>
        filter(target_desc_bait == loc_bait & target_desc_prey == loc_prey)
      permuts = permuts_df$likelihood_prey
      if (length(permuts) != 0) {
        probs[1:length(permuts)] = permuts
        ms = T
      } else {
        ms = F
      }
      ecdf_prob = ecdf(probs)
      observed_value = 0
      observed_value_row = observed_value_dfs |>
        filter(target_desc_bait == loc_bait & target_desc_prey == loc_prey)
      if (nrow(observed_value_row) != 0) {
        observed_value = observed_value_row$likelihood_prey
        biotin = T
      } else {
        biotin = F
      }
      quant_value = ecdf_prob(observed_value)
      quant_value_df[i,] = c(loc_bait, loc_prey, quant_value, ms, biotin)
      i = i +1
      print(paste0(i, " of ", expected_computations))
    }
  }
  return(quant_value_df)
}
            
            

ggplot(
  permut_prey_probability,
  aes(
    x=likelihood_prey,
    fill=target_desc_prey
  )
) + 
  geom_histogram() +
  facet_wrap(. ~ target_desc_bait) +
  theme(legend.position = "none")

ggplot(
  same_localisation_df,
  aes(
    x=log10(likelihood_prey),
    fill=target_desc_prey
  )
) + 
  geom_histogram(bins=30) +
  theme(legend.position = "none") + 
  facet_wrap(. ~ target_desc_bait) 





ggplot(
  bioID_prey_probability,
  aes(
    x=target_desc_bait,
    y=target_desc_prey,
    fill=likelihood_prey
  )
) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

permut_prey_probability_avg


ggplot(
  quant_value_df,
  aes(
    x=bait_localication,
    y=prey_localication,
    fill=quant_value
  )
) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + facet_wrap(ms_data ~ .)


function()
