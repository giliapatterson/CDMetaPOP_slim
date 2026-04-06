library(tidyverse)
library(glue)

# Directory
directory = "./"
# Number of patches to plot
num_plot = 3
# Runvars
run_info <- read_csv(paste0(directory, "RunVars.csv")) |>
  mutate(run = row_number() - 1) |>
  uncount(mcruns, .id = 'mc') |>
  mutate(mc = mc - 1,
         rep = row_number() - 1) |>
  mutate(output_folder = glue("{directory}slim_output/run{run}batch0mc{mc}species0/"),
         subpop_file = paste0(output_folder, "summary_popAllTime.csv"),
         class_file = paste0(output_folder, "summary_classAllTime.csv"),
         qtl_overall_file = paste0(output_folder, "QTL_overall.csv"),
         qtl_subpops_file = paste0(output_folder, "QTL_subpops.csv")) |>
  select(!any_of(c('output_years', 'cdclimgentime', 'mutation_output_years', 'mutation_output_subpops')))

# Read in results
subpops_unsplit = read_csv(pull(run_info, subpop_file), id = "subpop_file") |>
  full_join(run_info, by = "subpop_file")
npatches = subpops_unsplit |> slice(1) |> pull(K) |> str_split_1(fixed("|")) |> length() - 1
plot_patches <- sample(0:(npatches - 1), min(num_plot, npatches))
subpops <- subpops_unsplit |>
  mutate(PatchID = paste(c("all", 0:(npatches - 1)), collapse = "|")) |>
  separate_longer_delim(everything(), "|") |>
  filter(PatchID %in% plot_patches)
classes = read_csv(pull(run_info, class_file), id = "class_file") |>
  full_join(run_info, by = "class_file") |>
  separate_longer_delim(everything(),"|")

qtl_overall_all_reps <- read_csv(pull(run_info, qtl_overall_file), id = "qtl_overall_file") |>
  full_join(run_info, by = "qtl_overall_file") |>
  mutate(PatchID = "all")
# Take mean for each year and rep
# First check if any reps didn't finish
not_finished <- qtl_overall_all_reps |>
  group_by(rep) |>
  reframe(year = max(year)) |>
  filter(year != max(year)) |>
  pull(rep)
if(length(not_finished) > 0){
  print(glue("Reps {not_finished} did not finish and may cause plots to not accurately represent results."))
}
qtl_overall <- qtl_overall_all_reps |>
  group_by(dispersal_prob, PatchID, year, sizecontrol) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) 

qtl_subpops_all_reps <- read_csv(pull(run_info, qtl_subpops_file), id = "qtl_subpops_file") |>
  full_join(run_info, by = "qtl_subpops_file") |>
  mutate(PatchID = factor(as.integer(PatchID)))
# First check if any reps didn't finish
not_finished <- qtl_subpops_all_reps |>
  group_by(rep) |>
  reframe(year = max(year)) |>
  filter(year != max(year)) |>
  pull(rep)
if(length(not_finished) > 0){
  print(glue("Reps {not_finished} did not finish and may cause plots to not accurately represent results."))
}

subpops_summary <- qtl_subpops_all_reps |>
  group_by(PatchID, sizecontrol) |>
  reframe(popsize_decline = popsize[year == 850] - popsize[year == 800],
          popsize_recovery = popsize[year == 1500] - popsize[year == 800],
          allele_decline = Topt_alleles[year == 850] - Topt_alleles[year == 800],
          allele_increase = Topt_alleles[year == 1500] - Topt_alleles[year == 850],
          phenotype_change = avg_Topt[year == 1500] - avg_Topt[year == 800],
          pre_change_phenotype = avg_Topt[year == 800],
          post_change_phenotype = avg_Topt[year == 1500],
          phenotype_match_pre = avg_Topt[year == 800] - temperature[year == 800],
          phenotype_match_post = avg_Topt[year == 1500] - temperature[year == 1500])
subpop_means <- subpops_summary |>
  group_by(PatchID, sizecontrol) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) 

# Take mean for each year rep and patch
qtl_subpops <- qtl_subpops_all_reps |>
  group_by(PatchID, sizecontrol, year)|>
  filter(PatchID %in% plot_patches) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) 

ggplot(qtl_subpops, aes(x = year, y = popsize, color = sizecontrol)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Population size")
ggsave("qtl_plots/patch_popsize.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_overall, aes(x = year, y = popsize, color = sizecontrol)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Population size")
ggsave("qtl_plots/overall_popsize.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_subpops, aes(x = year, y = avg_Topt, color = PatchID)) +
  geom_point(aes(y = avg_Topt), size = 0.5, alpha = 0.5) +
  facet_grid(~sizecontrol, labeller = label_both) +
  scale_color_viridis_d(name = "PatchID") +
  xlab("Year") +
  ylab("Thermal Optimum (Degrees Celsius)") +
  theme_dark(base_size = 14)
ggsave("qtl_plots/Topt.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_subpops, aes(x = year, y = avg_epsilon, color = PatchID)) +
  geom_point(aes(y = avg_Topt), size = 0.5, alpha = 0.5) +
  facet_grid(~sizecontrol, labeller = label_both) +
  scale_color_viridis_d(name = "PatchID") +
  xlab("Year") +
  ylab("epsilon (Degrees Celsius)") +
  theme_dark(base_size = 14)
ggsave("qtl_plots/epsilon.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_subpops, aes(x = year, y = Topt_alleles, color = PatchID)) +
  geom_point(alpha = 0.2, position = position_jitter()) +
  facet_grid(~sizecontrol, labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Number of alleles underlying thermal optimum") +
  scale_color_viridis_d(name = "PatchID")
ggsave("qtl_plots/patch_alleles.png", width = 7, height = 6, create.dir = TRUE)
ggplot(qtl_overall, aes(x = year, y = Topt_alleles)) +
  geom_point(alpha = 0.2) +
  facet_grid(~sizecontrol, labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Number of alleles underlying thermal optimum") +
  scale_color_viridis_d(name = "Dispersal probability")
ggsave("qtl_plots/overall_Topt_alleles.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_overall, aes(x = year, y = epsilon_alleles, color = sizecontrol)) +
  geom_point(alpha = 0.2) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Number of alleles underlying epsilon")
ggsave("qtl_plots/overall_epsilon_alleles.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_subpops, aes(x = year, y = Topt_heritability, color = PatchID)) +
  geom_point(alpha = 0.2) +
  facet_grid(~sizecontrol, labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Heritability") +
  scale_color_viridis_d(name = "PatchID")
ggsave("qtl_plots/patch_heritability.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_overall, aes(x = year, y = Topt_heritability, color = sizecontrol)) +
  geom_point(alpha = 0.2) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Topt Heritability")
ggsave("qtl_plots/overall_Topt_heritability.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_overall, aes(x = year, y = epsilon_heritability, color = sizecontrol)) +
  geom_point(alpha = 0.2) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("epsilon Heritability")
ggsave("qtl_plots/overall_epsilon_heritability.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_subpops, aes(x = year, y = neutral_pi, color = PatchID)) +
  geom_point(alpha = 0.2) +
  facet_grid(~sizecontrol, labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Neutral genetic diversity (pi)") +
  scale_color_viridis_d(name = "PatchID")
ggsave("qtl_plots/patch_neutral_genetic_diversity.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_overall, aes(x = year, y = neutral_pi, color =sizecontrol)) +
  geom_point(alpha = 0.2) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Neutral genetic diversity (pi)")
ggsave("qtl_plots/overall_neutral_genetic_diversity.png", width = 7, height = 6, create.dir = TRUE)


ggplot(qtl_subpops, aes(x = year, y = overall_pi, color = PatchID)) +
  geom_point(alpha = 0.2) +
  facet_grid(~sizecontrol, labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Overall genetic diversity (pi)") +
  scale_color_viridis_d(name = "PatchID")
ggsave("qtl_plots/patch_overall_genetic_diversity.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_overall, aes(x = year, y = overall_pi, color = sizecontrol)) +
  geom_point(alpha = 0.2) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Overall genetic diversity (pi)")
ggsave("qtl_plots/overall_genetic_diversity.png", width = 7, height = 6, create.dir = TRUE)

