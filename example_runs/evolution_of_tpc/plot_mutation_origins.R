library(tidyverse)
library(glue)

# Directory
directory = "parameters/"
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

# Read in results and filter to six focal groups
focal_groups <- c(0, 3, 6, 9, 12, 15)
mutation_origins = read_csv(pull(run_info, mutation_origins_file), id = "mutation_origins_file") |>
  filter(subpop %in% focal_groups) |>
  full_join(run_info, by = "mutation_origins_file") |>
  group_by(subpop, tick, mutation_origins_file, environment, dispersal_distance) |>
  reframe(pick(everything()), Topt_prop = Topt_count/sum(Topt_count),
          neutral_prop = neutral_count/sum(neutral_count),
          epsilon_prop = epsilon_count/sum(epsilon_count))

## Group origin patches into temperature groups and summarize by each group
ngroups = 6
# Patch info
patch_info <- read_csv(glue("{directory}slim_output/slim_parameters/run1/PatchVars_slim.csv")) |>
  filter(year == 1) |>
  mutate(temp_group = findInterval(GrowthTemperatureBack,
                                   quantile(GrowthTemperatureBack,
                                            probs = seq(0, 1, length.out = ngroups + 1)) + c(rep(0, ngroups), 0.1))) |>
  select('PatchID', 'temp_group')
# Add patch info to mutation origins
vars_to_sum <- c('Topt_freq_origin_tick', 'epsilon_freq_origin_tick', 'neutral_freq_origin_tick')
temp_groups <- mutation_origins |>
  left_join(patch_info, by = join_by('origin_subpop' == 'PatchID')) |>
  rename(origin_temp_group = temp_group) |>
  group_by(environment, dispersal_distance, tick, origin_temp_group, subpop) |>
  reframe(across(all_of(vars_to_sum), \(x) sum(x))) |>
  filter(!is.na(Topt_freq_origin_tick))

ggplot(temp_groups, aes(x = tick, y = Topt_freq_origin_tick, color = factor(origin_temp_group))) +
  geom_point() +
  geom_smooth() +
  facet_grid(dispersal_distance+environment~subpop, scale = 'free_y') +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Proportion of variants underlying Topt')
ggsave("mutation_origins/Topt_variation.png", width = 7, height = 6, create.dir = TRUE)

ggplot(temp_groups, aes(x = tick, y = neutral_freq_origin_tick, color = factor(origin_temp_group))) +
  geom_point() +
  facet_grid(dispersal_distance+environment~subpop, scale = 'free_y') +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Proportion of neutral variants')
ggsave("mutation_origins/neutral_variation.png", width = 7, height = 6, create.dir = TRUE)

## Plot proportion of new mutations
vars_to_unique = c('Topt_prop_new', 'neutral_prop_new')
new_mutations <- mutation_origins |>
  group_by(subpop, environment, dispersal_distance, tick) |>
  reframe(across(all_of(vars_to_unique), \(x) x[1])) |>
  filter(!is.na(Topt_prop_new)) |>
  pivot_longer(cols = any_of(vars_to_unique),
               values_to = 'proportion_new')

ggplot(new_mutations, aes(x = tick, y = proportion_new)) +
  geom_point() +
  facet_grid(dispersal_distance+environment~subpop+name) +
  scale_color_brewer(name = "Origin", type = "seq", direction = 1, palette = 'Reds') +
  theme_dark() +
  xlab('Year') +
  ylab('Proportion')
ggsave("mutation_origins/new_mutations.png", width = 7, height = 6, create.dir = TRUE)


## Plot population sizes
qtl_with_temp_group <- read_csv(pull(run_info, qtl_subpops_file), id = "qtl_subpops_file") |>
  full_join(run_info, by = "qtl_subpops_file") |>
  rename(tick = year) |>
  left_join(patch_info)

focal_popsizes <- qtl_with_temp_group |>
  filter(PatchID %in% focal_groups)

## Plot population sizes for focal patches
ggplot(focal_popsizes, aes(x = tick, y = popsize)) +
  geom_point() +
  facet_grid(environment+dispersal_distance~PatchID) +
  theme_bw() +
  xlab('Year') +
  ylab('Total population size')
ggsave("mutation_origins/popsizes_by_patch.png", width = 7, height = 6, create.dir = TRUE)


## Plot population sizes for temperature groups
popsizes_summary <- qtl_with_temp_group |>
  group_by(temp_group, tick, environment, dispersal_distance) |>
  summarize(popsize = sum(popsize))
  
ggplot(popsizes_summary, aes(x = tick, y = popsize)) +
  geom_point() +
  geom_smooth() +
  facet_grid(environment+dispersal_distance~temp_group) +
  theme_bw() +
  xlab('Year') +
  ylab('Total population size')
ggsave("mutation_origins/popsizes_by_temperature.png", width = 7, height = 6, create.dir = TRUE)

## Plot overall population size
qtl_overall <- read_csv(pull(run_info, qtl_overall_file), id = "qtl_overall_file") |>
  full_join(run_info, by = "qtl_overall_file")

ggplot(qtl_overall, aes(x = year, y = popsize)) +
  geom_point() +
  geom_smooth() +
  facet_grid(environment~dispersal_distance) +
  theme_bw() +
  xlab('Year') +
  ylab('Total population size')
ggsave("mutation_origins/popsizes_by_temperature.png", width = 7, height = 6, create.dir = TRUE)