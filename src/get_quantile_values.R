library(tidyverse)
library(ggplot2)

quant_value_df = read.table(
  "work_folder/bioid_quantile/bait_all.csv",
  sep = "\t",
  header = T
)

quant_value_df$same_local = quant_value_df$target_desc_bait == quant_value_df$target_desc_prey

quant_value_df |>
  filter(in_bioid=="True") |>
  filter(in_permutation == "True") -> all_seen

all_seen$z <- (all_seen$observed_value - all_seen$probability_mean)/all_seen$probability_std


ggplot(
  all_seen,
  aes(
    x=target_desc_bait,
    y=target_desc_prey,
    fill=z
  )
) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "green") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) + 
  theme(
    legend.position = "bottom"
  )
  
ggplot(
  quant_value_df,
  aes(
    fill=same_local,
    x=quantile_value
  )
) + 
  geom_histogram() 
  facet_wrap(target_desc_prey ~ .) +
  theme(legend.position = "none")
  

  
ggplot(
  all_seen,
  aes(
    x = log(probability_mean),
    y = log(observed_value),
#    color=target_desc_bait
  )
) +
  xlab("ln(P(p_loc|b_loc)) of other MS detection") +
  ylab("ln(P(p_loc|b_loc)) of BioID pulldown ") +
  geom_point() +
  geom_density2d() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme(
    legend.position="none"
  )
  facet_wrap(target_desc_prey ~.) 

  
sds = 50  
ggplot(  
  all_seen |> filter(between(z, -sds, sds)),
  aes(
    x=z
  )) +
    geom_histogram()
                   
  