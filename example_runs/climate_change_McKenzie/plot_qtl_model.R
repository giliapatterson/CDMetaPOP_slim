library(tidyverse)
library(glue)

# Directory
directory = "climate_change_McKenzie/"
# Runvars
run_info <- read_csv(paste0(directory, "RunVars.csv")) |>
  mutate(dispersal_prob = as.numeric(str_extract(Popvars, "\\d+\\.?\\d*"))) |>
  uncount(mcruns, .id = 'mc') |>
  mutate(mc = mc - 1,
         rep = row_number()) |>
  mutate(output_folder = glue("{directory}slim_output/run{row_number() -1}batch0mc{mc}species0/"),
         subpop_file = paste0(output_folder, "summary_popAllTime.csv"),
         class_file = paste0(output_folder, "summary_classAllTime.csv"),
         qtl_overall_file = paste0(output_folder, "QTL_overall.csv"),
         qtl_subpops_file = paste0(output_folder, "QTL_subpops.csv")) |>
  select(!any_of(c('output_years', 'cdclimgentime')))

# Read in results
subpops_unsplit = read_csv(pull(run_info, subpop_file), id = "subpop_file") |>
  full_join(run_info, by = "subpop_file")
npatches = subpops_unsplit |> slice(1) |> pull(K) |> str_split_1(fixed("|")) |> length() - 1
subpops <- subpops_unsplit |>
  mutate(PatchID = paste(c("all", 0:(npatches - 1)), collapse = "|")) |>
  separate_longer_delim(everything(), "|")
classes = read_csv(pull(run_info, class_file), id = "class_file") |>
  full_join(run_info, by = "class_file") |>
  separate_longer_delim(everything(),"|")

qtl_overall <- read_csv(pull(run_info, qtl_overall_file), id = "qtl_overall_file") |>
  full_join(run_info, by = "qtl_overall_file") |>
  mutate(PatchID = "all")
qtl_subpops <- read_csv(pull(run_info, qtl_subpops_file), id = "qtl_subpops_file") |>
  full_join(run_info, by = "qtl_subpops_file") |>
  mutate(PatchID = factor(as.integer(PatchID)))

ggplot(qtl_subpops, aes(x = year, y = popsize, color = factor(dispersal_prob))) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = qtl_overall, alpha = 0.5, size = 0.5) +
  theme_dark(base_size = 14) +
  facet_wrap(~PatchID, scales = "free", labeller = label_both) +
  xlab("Year") +
  ylab("Population size") +
  scale_color_viridis_d(name = "Dispersal probability")
ggsave("climate_change_McKenzie/qtl_plots/popsize.png", width = 7, height = 6, create.dir = TRUE)

ggplot(qtl_subpops, aes(x = year, y = avg_phenotype)) +
  geom_point(aes(color = 'Thermal Optimum'), size = 0.5, alpha = 0.5) +
  geom_line(aes(y = qtl_environment, color = 'Temperature')) +
  facet_grid(dispersal_prob~PatchID, labeller = label_both) +
  scale_color_viridis_d(name = "") +
  xlab("Year") +
  ylab("Degrees Celcius") +
  theme_dark(base_size = 14)
ggsave("climate_change_McKenzie/qtl_plots/phenotype.png", width = 7, height = 6, create.dir = TRUE)


ggplot(qtl_subpops, aes(x = year, y = qtl_alleles, color = factor(dispersal_prob))) +
  geom_point(alpha = 0.2, position = position_jitter()) +
  geom_point(data = qtl_overall, alpha = 0.2, position = position_jitter()) +
  facet_wrap(~PatchID, scales = "free", labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Number of alleles underlying thermal optimum") +
  scale_color_viridis_d(name = "Dispersal probability")
ggsave("climate_change_McKenzie/qtl_plots/alleles.png", width = 7, height = 6, create.dir = TRUE)


ggplot(qtl_subpops, aes(x = year, y = heritability, color = factor(dispersal_prob))) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = qtl_overall, alpha = 0.5, size = 0.5) +
  facet_wrap(~PatchID, scales = "free", labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Heritability") +
  scale_color_viridis_d(name = "Dispersal probability")
ggsave("climate_change_McKenzie/qtl_plots/heritability.png", width = 7, height = 6, create.dir = TRUE)


ggplot(qtl_subpops, aes(x = year, y = neutral_pi, color = factor(dispersal_prob))) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = qtl_overall, alpha = 0.5, size = 0.5) +
  facet_wrap(~PatchID, scales = "free", labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Neutral genetic variation (pi)") +
  scale_color_viridis_d(name = "Dispersal probability")
ggsave("climate_change_McKenzie/qtl_plots/neutral_pi.png", width = 7, height = 6, create.dir = TRUE)


ggplot(qtl_subpops, aes(x = year, y = overall_pi, color = factor(dispersal_prob))) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = qtl_overall, alpha = 0.5, size = 0.5) +
  facet_wrap(~PatchID, scales = "free", labeller = label_both) +
  theme_dark(base_size = 14) +
  xlab("Year") +
  ylab("Overall genetic variation (pi)") +
  scale_color_viridis_d(name = "Dispersal probability")
ggsave("climate_change_McKenzie/qtl_plots/genetic_diversity.png", width = 7, height = 6, create.dir = TRUE)

  
