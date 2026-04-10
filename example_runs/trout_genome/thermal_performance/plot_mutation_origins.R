library(tidyverse)
library(glue)

# Directory
directory = "./"
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
         qtl_subpops_file = paste0(output_folder, "QTL_subpops.csv"),
         mutation_origins_file = paste0(output_folder, "mutation_origins.csv")) |>
  select(!any_of(c('output_years', 'cdclimgentime', 'mutation_output_years', 'mutation_output_subpops')))

# Read in results
mutation_origins = read_csv(pull(run_info, mutation_origins_file), id = "mutation_origins_file") |>
  full_join(run_info, by = "mutation_origins_file") |>
  group_by(subpop, tick, mutation_origins_file) |>
  reframe(pick(everything()), Topt_prop = Topt_count/sum(Topt_count),
          neutral_prop = neutral_count/sum(neutral_count),
          epsilon_prop = epsilon_count/sum(epsilon_count))
# Sort focal subpops and origin subpops into cold, medium-cold, medium, medium-warm, and warm
# Patch info
patch_info <- read_csv(glue("{directory}slim_output/slim_parameters/run1/PatchVars_slim.csv")) |>
  filter(year == 1) |>
  mutate(temp_group = findInterval(GrowthTemperatureBack, quantile(GrowthTemperatureBack, probs = seq(0, 1, length.out = 5)))) |>
  select('PatchID', 'temp_group')
group_names = c("Cold patch", "Cool patch", "Medium patch", "Warm patch", "Hot patch")
# Add patch info to mutation origins
temp_groups <- mutation_origins |> left_join(patch_info, by = join_by('origin_subpop' == 'PatchID')) |>
  rename(origin_temp_group = temp_group) |>
  left_join(patch_info, by = join_by('subpop' == 'PatchID')) |>
  mutate(focal_temp_group = temp_group) |>
  group_by(sizecontrol, tick, origin_temp_group, focal_temp_group) |>
  reframe(Topt_count_group = sum(Topt_count),
          epsilon_count_group = sum(epsilon_count),
          neutral_count_group = sum(neutral_count),
          Topt_mean_contribution_group = mean(Topt_mean_contribution),
          epsilon_mean_contribution_group = mean(epsilon_mean_contribution),
          Topt_VA_group = sum(Topt_VA),
          epsilon_VA_group = sum(epsilon_VA),
          pick(everything())) |>
  group_by(sizecontrol, tick, focal_temp_group) |>
  mutate(Topt_prop_group = Topt_count_group/sum(Topt_count_group),
         epsilon_prop_group = epsilon_count_group/sum(epsilon_count_group),
         neutral_prop_group = neutral_count_group/sum(neutral_count_group),
         Topt_prop_VA = Topt_VA_group/sum(Topt_VA_group),
         epsilon_prop_VA = epsilon_VA_group/sum(epsilon_VA_group),
         pick(everything())) |>
  group_by(sizecontrol, tick, focal_temp_group, origin_temp_group) |>
  filter(tick %in% c(800, 850, 1500))
temp_groups_with_missing <- temp_groups |>
  group_by(sizecontrol, focal_temp_group, origin_temp_group) |>
  reframe(tick = unique(temp_groups$tick)[which(!(unique(temp_groups$tick) %in% unique(tick)))]) |>
  full_join(temp_groups)|> 
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  mutate(origin_group_name = fct_reorder(group_names[origin_temp_group], origin_temp_group),
         focal_group_name = fct_reorder(group_names[focal_temp_group], focal_temp_group),
         Sizecontrol = paste('Sizecontrol: ', sizecontrol))
ggplot(temp_groups_with_missing, aes(x = tick, y = Topt_prop_group, color = origin_group_name)) +
  geom_point() +
  geom_line() +
  facet_grid(sizecontrol~focal_group_name, scale = 'free_y') +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Proportion of variants underlying Topt')
ggsave("mutation_origins/Topt_variation.png", width = 7, height = 6, create.dir = TRUE)

ggplot(temp_groups_with_missing, aes(x = tick, y = Topt_mean_contribution_group, color = origin_group_name)) +
  geom_point() +
  geom_line() +
  facet_grid(sizecontrol~focal_group_name, scale = 'free_y') +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Mean contribution to Topt')
ggsave("mutation_origins/Topt_contribution.png", width = 7, height = 6, create.dir = TRUE)

ggplot(temp_groups_with_missing, aes(x = tick, y = epsilon_prop_group, color = origin_group_name)) +
  geom_point() +
  geom_line() +
  facet_grid(sizecontrol~focal_group_name, scale = 'free_y') +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Proportion of variants underlying epsilon')
ggsave("mutation_origins/epsilon_variation.png", width = 7, height = 6, create.dir = TRUE)

ggplot(temp_groups_with_missing, aes(x = tick, y = epsilon_mean_contribution_group, color = origin_group_name)) +
  geom_point() +
  geom_line() +
  facet_grid(sizecontrol~focal_group_name, scale = 'free_y') +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Mean contribution to epsilon')
ggsave("mutation_origins/epsilon_contribution.png", width = 7, height = 6, create.dir = TRUE)

ggplot(temp_groups_with_missing, aes(x = tick, y = neutral_prop_group, color = origin_group_name)) +
  geom_point() +
  geom_line() +
  facet_grid(sizecontrol~focal_group_name, scale = 'free_y') +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Proportion of neutral variants')
ggsave("mutation_origins/neutral_variation.png", width = 7, height = 6, create.dir = TRUE)

## Add QTL patch info
qtl_with_temp_group <- read_csv(pull(run_info, qtl_subpops_file), id = "qtl_subpops_file") |>
  full_join(run_info, by = "qtl_subpops_file") |>
  rename(tick = year) |>
  left_join(patch_info) |>
  mutate(patch_name = fct_reorder(group_names[temp_group], temp_group))
overall_group_names = c("Cold patches", "Cool patches", "Medium patches", "Warm patches", "Hot patches")
popsizes_summary <- qtl_with_temp_group |> group_by(temp_group, tick, sizecontrol) |>
  summarize(popsize = mean(popsize)) |>
  mutate(patch_name = fct_reorder(overall_group_names[temp_group], temp_group))
  
ggplot(popsizes_summary, aes(x = tick, y = popsize)) +
  geom_point() +
  geom_smooth() +
  facet_grid(sizecontrol~patch_name) +
  theme_bw() +
  xlab('Year') +
  ylab('Mean population size')
ggsave("mutation_origins/popsizes_by_temperature.png", width = 7, height = 6, create.dir = TRUE)

ggplot(filter(qtl_with_temp_group, PatchID %in% unique(temp_groups_with_missing$subpop)), aes(x = tick, y = avg_Topt)) +
  geom_point() +
  geom_smooth() +
  facet_grid(sizecontrol~patch_name) +
  theme_bw() +
  xlab('Year') +
  ylab('Average Topt')
ggsave("mutation_origins/focal_patches_Topt.png", width = 7, height = 6, create.dir = TRUE)

ggplot(filter(qtl_with_temp_group, PatchID %in% unique(temp_groups_with_missing$subpop)), aes(x = tick, y = avg_epsilon)) +
  geom_point() +
  geom_smooth() +
  facet_grid(sizecontrol~patch_name) +
  theme_bw() +
  xlab('Year') +
  ylab('Average epsilon')
ggsave("mutation_origins/focal_patches_epsilon.png", width = 7, height = 6, create.dir = TRUE)


