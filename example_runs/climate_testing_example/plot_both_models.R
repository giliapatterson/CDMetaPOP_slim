library(tidyverse)

split_variable <- function(variable_name, df, separated_by_string, includes_all){
  strings = df |> pull(variable_name)
  split = str_split_fixed(strings, "\\|", Inf)
  if(includes_all & ncol(split) > 1){colnames(split) = c("all", 0:(ncol(split) - 2))}
  if(!includes_all & ncol(split) > 1){colnames(split) = 0:(ncol(split)-1)}
  if(ncol(split) == 1){colnames(split) = "all"}
  split <- as_tibble(split) |> mutate(across(all_of(colnames(split)), as.numeric))|>
    mutate(Year = df$Year) |>
    pivot_longer(cols = !Year,
                 names_to = separated_by_string,
                 values_to = variable_name)
  return(split)
}

process_summary_pop <- function(results, separated_by_string = "Patch", includes_all = TRUE){
  variables = colnames(select(results, !Year))
  split_variables = map(variables, \(x) split_variable(x, results, separated_by_string, includes_all))
  out_df <- reduce(split_variables, full_join)
  out_df[out_df == -1] = NA
  return(out_df)
}
slim_results <- read_csv("slim_output/run0batch0mc0species0/summary_popAllTime.csv")
cdmp_results <- read_csv("cdmetapop_results/run0batch0mc0species0/summary_popAllTime.csv")
slim <- process_summary_pop(slim_results) |> mutate(Year = Year - 1)
cdmp <- process_summary_pop(cdmp_results)
all_patches <- bind_rows(list(slim = slim, cdmp = cdmp), .id = "method") |> filter(!is.na(K))

all <- filter(all_patches, Patch == "all")
patches <- filter(all_patches, Patch != 'all')

slim_class_results <- read_csv("slim_output/run0batch0mc0species0/summary_classAllTime.csv")
cdmp_class_results <- read_csv("cdmetapop_results/run0batch0mc0species0/summary_classAllTime.csv")
slim_class <- process_summary_pop(slim_class_results, separated_by_string = "Class", includes_all = FALSE)
cdmp_class <- process_summary_pop(cdmp_class_results, separated_by_string = "Class", includes_all = FALSE)
all_class <- bind_rows(list(slim = slim_class, cdmp = cdmp_class), .id = "method") 

ggplot(all, aes(x = Year, y = K, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = N_Initial, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(patches, aes(x = Year, y = N_Initial, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  facet_wrap("Patch", scale = "free_y")
ggplot(patches, aes(x = Year, y = N_beforePacking_AddAge0s, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  facet_wrap("Patch", scale = "free_y")
ggplot(patches, aes(x = Year, y = N_Emigration, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  facet_wrap("Patch", scale = "free_y")
ggplot(patches, aes(x = Year, y = N_Females, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  facet_wrap("Patch", scale = "free_y")
ggplot(patches, aes(x = Year, y = N_Males, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  facet_wrap("Patch", scale = "free_y")

ggplot(all, aes(x = Year, y = PackingDeaths_Immigration, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = PackingDeaths_Emigration, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all_class, aes(x = Year, y = N_Initial_Age, color = method)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Ages, scales = "free")
ggplot(all_class, aes(x = Year, y = N_Initial_Class, color = method)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Class, scales = "free")


ggplot(all_class, aes(x = Year, y = PackingDeaths_Emigration, color = method)) +
  geom_point(alpha = 0.4) +
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scales = "free")
ggplot(all_class, aes(x = Year, y = PackingDeaths_Immigration, color = method)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scales = "free")
ggplot(all_class, aes(x = Year, y = N_BeforePacking_AddAge0s, color = method)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Class, scales = "free")
ggplot(all_class, aes(x = Year, y = N_GrowthBack, color = method)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~Class, scales = "free")

ggplot(all, aes(x = Year, y = PopSizes_Mean, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)
ggplot(all, aes(x = Year, y = PopSizes_Std, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = GrowthRate, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  ylim(0, 1.3)


ggplot(patches, aes(x = Year, y = N_MatureFemales, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  facet_wrap(~Patch)
ggplot(patches, aes(x = Year, y = N_MatureMales, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5) +
  facet_wrap(~Patch)

ggplot(all, aes(x = Year, y = MatureCount/N_Initial, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)
ggplot(all, aes(x = Year, y = N_MatureFemales/N_Initial, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)
ggplot(all, aes(x = Year, y = N_MatureMales/N_Initial, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = EggLayEvents, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)
ggplot(all, aes(x = Year, y = EggDeaths, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = EggDeaths/N_MatureFemales, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = EggLayEvents/N_MatureFemales, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = EggLayEvents/N_MatureMales, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = N_Emigration, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = Births/EggLayEvents, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = Ho, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)
ggplot(all, aes(x = Year, y = He, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

ggplot(all, aes(x = Year, y = Alleles, color = method)) +
  geom_point(alpha = 0.5) + 
  geom_line(alpha = 0.5)

age_props <- all_class |> group_by(method, Year) |>
  summarize(N_Initial = sum(N_Initial_Age, na.rm = TRUE)) |>
  full_join(all_class) |>
  group_by(method, Year, Class) |>
  summarize(mean_N = mean(N_Initial_Age),
            prop = mean(N_Initial_Age/N_Initial),
            prop_class = mean(N_Initial_Class/N_Initial))

ggplot(age_props, aes(x = Year, y = prop, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scale = "free_y")
ggplot(age_props, aes(x = Year, y = prop_class, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scale = "free_y")

ggplot(all_class, aes(x = Year, y = AgeSize_Mean, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scale = "free_y")

ggplot(all_class, aes(x = Year, y = ClassSize_Mean, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scale = "free_y")

ggplot(all_class, aes(x = Year, y = PackingDeaths_Emigration, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scale = "free_y")

ggplot(all_class, aes(x = Year, y = PackingDeaths_Immigration, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5) +
  facet_wrap(~Class, scale = "free_y")

ggplot(all, aes(x = Year, y = Alleles, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5)


ggplot(all_patches, aes(x = Year, y = Ho, color = method)) +
  geom_point(alpha = 0.5)+ 
  geom_line(alpha = 0.5) +
  facet_wrap(~Patch)

all |> group_by(method) |> summarize(mean(N_MatureFemales))
patch_mature <- all_patches |> filter(Patch != "all") |> group_by(method, Patch) |> summarize(mature = mean(N_MatureFemales))

ggplot(patch_mature, aes(x = Patch, y = mature, color = method)) +
  geom_point()

class_sizes <- all_class |> group_by(method, Class) |>
  summarize(mean_size = mean(AgeSize_Mean, na.rm = TRUE),
            class_mean_size = mean(ClassSize_Mean, na.rm = TRUE),
            mean_size_std = mean(AgeSize_Std, na.rm = TRUE),
            class_size_std = mean(ClassSize_Std, na.rm = TRUE))

ggplot(class_sizes, aes(x = Class, y = mean_size, color = method)) +
  geom_point(alpha = 0.5)
ggplot(class_sizes, aes(x = Class, y = class_mean_size, color = method)) +
  geom_point(alpha = 0.5)
ggplot(class_sizes, aes(x = Class, y = mean_size_std, color = method)) +
  geom_point(alpha = 0.5)
ggplot(class_sizes, aes(x = Class, y = class_size_std, color = method)) +
  geom_point(alpha = 0.5)


