##***************************
## Assignment 1 Code
##
## Group 8
##
## 2026-10-27
##
##***************************

rm(list = ls())

## Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

# Startup ends here

## Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)
#BA to mark suggestions made by assignment 2 partner

## Reading in data -----

df_trich_sp = read_csv("../data/trichinella_spiralis_te.csv")
df_trich_sp

df_trich_pap = read_csv("../data/GCA_001447755.1.csv")
df_trich_pap

## Filtering the data -----

ls_non_te = c("Simple_repeat", 
              "Low_complexity", 
              "Unknown", 
              "simple_repeat", 
              "Satellite", 
              "rRNA", 
              "snRNA")

#BA#1 Suggestion: Combine identical ls_non_te_sp and ls_non_te_pap objects into one ls_non_te list, replace this new object in the right steps

#The above matching classes were identified as non transposable elements based on the results of a Copilot prompt 

ls_retransposons_sp = c("LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-RTE", "SINE/tRNA-Deu-L2", "LTR/Gypsy", "LTR/Pao", "LTR/Ngaro", "LTR/DIRS")

#The above matching classes were identified as retrotransposons in this df based on the results of a Copilot prompt 

df_spiralis_retro = df_trich_sp |> 
  filter(matching_class %notin% ls_non_te) |>
  filter(matching_class %in% ls_retransposons_sp) 

df_spiralis_retro

#BA#2 Suggestion: I would suggest combining the two filtering steps and to create only one new filtered data frame for further analysis which could make the code easier to follow and keep the environment tidy. This would also make the code more aligned to your biological question as you are primarily analyzing retrotransposons between the species The same applies to the T. papuae filtering steps as well. 

#Filtering out all the non TEs and keeping only the retrotransposons in the Trichinella spiralis data frame

# ls_non_te_pap = c("Simple_repeat", "Low_complexity", "Unknown", "simple_repeat", "Satellite", "rRNA", "snRNA") -> this object was replaced by the shared ls_non_te list

ls_retransposons_pap = c("LINE", "LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-RTE", "LTR/Gypsy", "LTR/Pao", "LTR/Ngaro")

#The above matching classes were identified as retrotransposons in this df based on the results of a Copilot prompt 

df_papuae_retro = df_trich_pap |> 
  filter(matching_class %notin% ls_non_te) |>
  filter(matching_class %in% ls_retransposons_pap)

df_papuae_retro

#Filtering out all the non TEs and keeping only the retrotransposons in the Trichinella papuae data frame

## Combining the data frames -----

df_trich_sp_final = df_spiralis_retro |> 
  mutate(species = "T. spiralis") |> 
  mutate(qry_left = str_remove_all(string = qry_left, pattern = "-")) |> 
  mutate(no_bp_in_complement = 
           str_remove_all(string = no_bp_in_complement, pattern = "-")) |> 
  mutate(in_repeat_end = 
           str_remove_all(string = in_repeat_end, pattern = "-"))

df_trich_pap_final = df_papuae_retro |> 
  mutate(species = "T. papuae") |> 
  mutate(qry_left = str_remove_all(string = qry_left, pattern = "[()]")) |> 
  mutate(no_bp_in_complement = 
           str_remove_all(string = no_bp_in_complement, pattern = "[()]")) |> 
  mutate(in_repeat_end = 
           str_remove_all(string = in_repeat_end, pattern = "[()]"))

#Adding a column indicating which species before combining
#Formatting so that the bind_rows function can run without error.
#Formatting done was removing "()" and "-" from values in those 3 columns

df_trich_combined = bind_rows(df_trich_sp_final, df_trich_pap_final)

#Combined data frames for more analysis options

## Figures -----

df_trich_sp_retro |> 
  ggplot(aes(perc_insert)) +
  geom_histogram() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Spiralis")
ggsave("../output/fig1_perc insertion spiralis_histogram.png") 

df_trich_pap_retro |> 
  ggplot(aes(perc_insert)) +
  geom_histogram() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Papuae")
ggsave("../output/fig2_perc insertion papuae_histogram.png") 

df_trich_sp_retro |> 
  ggplot(aes(perc_insert)) +
  geom_freqpoly() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Spiralis")
ggsave("../output/fig3_perc insertion spiralis_poly.png") 

df_trich_pap_retro |> 
  ggplot(aes(perc_insert)) +
  geom_freqpoly() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Papuae")
ggsave("../output/fig4_perc insertion papuae_poly.png") 

df_trich_combined |> 
  ggplot(aes(x = perc_insert, colour = species)) +
  geom_freqpoly() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Papuae and Spiralis", colour = "Species") +
  theme_minimal()
ggsave("../output/fig5_perc insertion comparison_poly.png") 

#BA#3 Suggestion: For your comparative figure, the minor edits made to the code above to overlay the polygraphs over one another to make the visual comparison easier by removing the faceting and adding colors for the species. Adding ggsave commands to all figures  

## Potential novel code -----
  
summary(df_trich_sp_retro$perc_insert)
summary(df_trich_pap_retro$perc_insert)

df_trich_combined_lm = 
  lm(df_trich_combined$perc_insert~df_trich_combined$species)
summary(df_trich_combined_lm)

t.test(df_trich_pap_retro$perc_insert, df_trich_sp_retro$perc_insert)

#BA#4 Suggestion: Based on the figures generated, the perc insertion data is not a normal distribution so it does not align with the underlying assumptions of linear model generation and t tests, a Wilcoxon rank-sum test could be a better fit here to see if percentage insetion values are higher in one species vs. the other

wilcox.test(perc_insert ~ species, data = df_trich_combined)

#W = 130714, p-value = 0.4525 (not statistically significant, a=0.05)
#alternative hypothesis: true location shift is not equal to 0 