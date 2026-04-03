library(tidyverse)
library(glue)

## Popvars for size control = 'Y'
popvars <- read_csv("PopVars1.csv")

## thermal performance curves for a range of parameters including the ones used
tpc_df <- expand_grid(Topt = c(popvars$Topt_baseline - 5, popvars$Topt_baseline, popvars$Topt_baseline + 5),
                      epsilon = c(popvars$epsilon_baseline, popvars$epsilon_baseline + 5),
                      Topt_baseline = popvars$Topt_baseline,
                      epsilon_baseline = popvars$epsilon_baseline,
                      t = seq(-20, 28, length.out = 50),
                      cost_scale = c(0, 1)) |>
  mutate(Topt_cost = cost_scale*popvars$Topt_cost,
         epsilon_cost = cost_scale*popvars$epsilon_cost) |>
  mutate(penalty_scale = pmax(0.0, pmin(1 - abs(Topt - Topt_baseline)*Topt_cost - abs(epsilon - epsilon_baseline)*epsilon_cost, 1.0)),
         x = (t-Topt)/epsilon,
         P1 = exp(x)*(1-x)*penalty_scale,
         P = if_else(P1>0, P1, 0),
         Penalty = if_else(cost_scale == 0, "No", "Yes"))
ggplot(tpc_df, aes(x = t, y = P, color = Penalty)) +
  geom_line() +
  xlab("Temperature") +
  ylab("Performance") +
  facet_grid(epsilon~Topt, labeller = label_both, scales = "free_y")

## Von Bertalanffy growth
source('parameter_functions.R')

popvars <-  mutate(popvars, mature_eqn_slope_f = as.numeric(str_split_1(mature_eqn_slope, "~")[1]),
                   mature_eqn_slope_m = as.numeric(str_split_1(mature_eqn_slope, "~")[2]),
                   mature_eqn_int_f = as.numeric(str_split_1(mature_eqn_int, "~")[1]),
                   mature_eqn_int_m = as.numeric(str_split_1(mature_eqn_int, "~")[2]),
                   mature_age = gsub("age", "", mature_default),
                   mature_age_f = as.numeric(str_split_1(mature_age, "~")[1]),
                   mature_age_m = as.numeric(str_split_1(mature_age, "~")[2]),
                   mature_age = mature_age_f)

## Growth curves
ages <- c()
sizes <- c()
patch_temps <- c()
for(patch_temp in tpc_df$t){
  current_age <- 0
  current_size <- 49
  while(current_age <= 7){
    ages <- c(ages, current_age)
    sizes <- c(sizes, current_size)
    patch_temps <- c(patch_temps, patch_temp)
    current_age <- current_age + 1
    current_size <- grow_temperature(current_size, current_age, patch_temp, 365, popvars)
  }
}
growth_df <- tibble(ages, sizes, t = patch_temps)
ggplot(growth_df, aes(x = t, y = sizes, color = factor(ages))) +
  geom_line() +
  ylab("Body size (mm)") +
  xlab("Patch temperature\n(degrees Celsius)") +
  scale_colour_discrete(
    name = "Age",
  ) +
  theme_bw(base_size = 18)

# Add TPC to maturity as scaling factor
mat_tpc_df <- growth_df |>
  expand_grid(sex = c('F', 'M')) |>
  mutate(pmature = prob_mature_vector(sizes, sex, ages, popvars),
         mean_eggs = mean_num_eggs(sizes, popvars)*(sex == 'F')
  ) |>
  full_join(tpc_df, by = 't', relationship = 'many-to-many') |>
  mutate("Environmental effect only" = mean_eggs*pmature,
         "Environmental and genetic effects" = mean_eggs*pmature*P) |>
  group_by(ages, sex) |>
  reframe(pick(everything()), "Genetic effect only" = max(`Environmental effect only`)*P) |>
  pivot_longer(c("Environmental and genetic effects", "Genetic effect only"), names_to = "method", values_to = "fertility")

ggplot(mat_tpc_df |> filter(sex == 'F', ages == 7), aes(x = t, y = fertility, color = method)) +
  geom_line() +
  xlab("Patch temperature\n(degrees Celsius)") +
  ylab("Age 7 Fertility") +
  scale_colour_discrete(
    name = "Method",
  )  +
  facet_grid(epsilon*Penalty~Topt, labeller = label_both) +
  theme_bw(base_size = 14)

x = seq(0, 5, 0.1)
plot(x,  exp(x) -1)
plot(x, 1-dnorm(x, 0, 0.1))
distributions = tibble(x = seq(-4, 4, 0.1),
                       normal = dnorm(x, 0, 1),
                       exp10 = dexp(x, 10) + dexp(-x, 10)) |>
  pivot_longer(-x, names_to = "distribution", values_to = "density")

ggplot(distributions, aes(x = x, y = density, color = distribution)) +
  geom_line()

hist(rexp(1000, 10))