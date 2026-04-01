library(tidyverse)
library(glue)

## thermal performance curves
tpc_df <- expand_grid(Topt = c(10, 18), CTmax = c(12, 28), t = seq(-20, 28, 0.1),
                      range = 2,
                      cost = c(0.0001, 0.1)) |>
  mutate(epsilon = CTmax - Topt,
         penalty = epsilon*cost,
         x = (t-Topt)/epsilon,
         P1 = exp(x)*(1-x),
         P = if_else(P1>0, P1, 0),
         Fecundity = if_else(P-penalty > 0, P-penalty, 0),
         P1_integrated = exp(x+1)*(2-(x+1)) - exp(x-1)*(2-(x-1)),
         P_integrated = if_else(P1_integrated>0, P1_integrated, 0),
         Fecundity_int = if_else(P_integrated-penalty > 0, P_integrated-penalty, 0),
         Cost = glue("{cost}/degree"),
         Mortality = 1-if_else(P-penalty > 0, P-penalty, 0))
ggplot(tpc_df, aes(x = t, y = P, color = factor(CTmax))) +
  geom_line() +
  xlab("Temperature") +
  ylab("Performance") +
  facet_wrap(~Topt, labeller = label_both, scales = "free_y")
ggplot(tpc_df, aes(x = t, y = Fecundity, color = factor(CTmax))) +
  geom_line() +
  xlab("Temperature") +
  ylab("Fecundity") +
  facet_wrap(Cost~Topt, labeller = label_both, scales = "free_y")
ggplot(tpc_df, aes(x = t, y = Fecundity_int, color = factor(CTmax))) +
  geom_line() +
  xlab("Temperature") +
  ylab("Fecundity (integrated)") +
  facet_wrap(Cost~Topt, labeller = label_both, scales = "free_y")
ggplot(tpc_df, aes(x = t, y = Mortality, color = factor(CTmax))) +
  geom_line() +
  xlab("Temperature") +
  ylab("Mortality") +
  facet_wrap(Cost~Topt, labeller = label_both, scales = "free_y")

## Relationship between CTmax, cost, and 
max_f <- expand_grid(Topt = 10, CTmax = seq(10, 15, 0.1),
                      cost = seq(0, 0.5, 0.1), reference_CTmax = 10) |>
  mutate(epsilon = CTmax - Topt,
         reference_epsilon = reference_CTmax - Topt,
         penalty = (epsilon-reference_epsilon)*cost,
         max_F = if_else(1.0 - penalty > 0, 1.0-penalty, 0))
ggplot(max_f, aes(x = CTmax, y = max_F, color = factor(cost))) +
  geom_line()

## Von Bertalanffy growth
popvars <- read_csv("thermal_performance/PopVars1.csv")
source('thermal_performance/parameter_functions.R')

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
for(patch_temp in seq(0,20, length.out = 50)){
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
mat_tpc_df <- expand_grid(growth_df, Topt = c(5, 10, 15),
                          epsilon = c(5, 10),
                          cost = c(0, 0.01, 0.1)) |>
  mutate(CTmax = Topt + epsilon,
         penalty = epsilon*cost,
         x = (t-Topt)/epsilon,
         P1 = exp(x)*(1-x) - penalty,
         P = if_else(P1>0, P1, 0)) |>
  expand_grid(sex = c('F', 'M')) |>
  mutate(pmature = prob_mature_vector(sizes, sex, ages, popvars),
         mean_eggs = mean_num_eggs(sizes, popvars)*(sex == 'F'),
         "Environmental effect only" = mean_eggs*pmature,
         "Environmental and genetic effects" = mean_eggs*pmature*P
  ) |>
  group_by(ages, sex) |>
  reframe(pick(everything()), "Genetic effect only" = max(`Environmental effect only`)*P) |>
  pivot_longer(c("Environmental effect only", "Environmental and genetic effects", "Genetic effect only"), names_to = "method", values_to = "fertility")
ggplot(growth_tpc_df |> filter(ages == 7), aes(x = t, y = sizes)) +
  geom_line() +
  ylab("Age 7 body size (mm)") +
  xlab("Patch temperature\n(degrees Celsius)") +
  theme_bw(base_size = 18)
ggplot(mat_tpc_df |> filter(sex == 'F', ages == 7), aes(x = t, y = fertility, color = method)) +
  geom_line() +
  xlab("Patch temperature\n(degrees Celsius)") +
  ylab("Age 7 Fertility") +
  scale_colour_discrete(
    name = "Method",
  )  +
  facet_grid(epsilon*cost ~Topt, labeller = label_both) +
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