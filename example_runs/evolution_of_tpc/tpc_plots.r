library(tidyverse)
library(glue)

## Popvars for size control = 'Y'
popvars <- read_csv("parameters/PopVars1.csv")

## thermal performance curves for a range of parameters including the ones used
tpc_df <- expand_grid(Topt = c(popvars$Topt_baseline - 5, popvars$Topt_baseline, popvars$Topt_baseline + 5),
                      epsilon = popvars$epsilon_baseline,
                      Topt_baseline = popvars$Topt_baseline,
                      epsilon_baseline = popvars$epsilon_baseline,
                      t = seq(-20, 28, length.out = 50),
                      cost_scale = c(0)) |>
  mutate(Topt_cost = cost_scale*popvars$Topt_cost,
         epsilon_cost = cost_scale*popvars$epsilon_cost) |>
  mutate(penalty_scale = pmax(0.0, pmin(1 - abs(Topt - Topt_baseline)*Topt_cost - abs(epsilon - epsilon_baseline)*epsilon_cost, 1.0)),
         x = (t-Topt)/epsilon,
         P1 = exp(x)*(1-x)*penalty_scale,
         P = if_else(P1>0, P1, 0),
         Penalty = if_else(cost_scale == 0, "No", "Yes"),
         CTMax = Topt + epsilon)

## Environmental effect
environment_df = tibble(t = tpc_df$t,
                        best_environment = 10,
                        sd_environment = 5,
                        env_scaling_factor = dnorm(t, mean = best_environment, sd = sd_environment)/dnorm(best_environment, mean = best_environment, sd = sd_environment))
ggplot(environment_df, aes(x = t, y = env_scaling_factor)) +
  geom_line() +
  ylab("Scaling factor") +
  xlab("Patch temperature\n(degrees Celsius)") +
  scale_colour_discrete(
    name = "Age",
  ) +
  theme_bw(base_size = 18)

# Add environmental effect to TPC
tpc_and_environment <- tpc_df |>
  full_join(environment_df, by = 't', relationship = 'many-to-many') |>
  mutate("Environmental effect only" = env_scaling_factor,
         "Environmental and genetic effects" = P*env_scaling_factor,
         "Genetic effect only" = P) |>
  pivot_longer(c("Environmental effect only", "Environmental and genetic effects", "Genetic effect only"),
               names_to = "method", values_to = "performance")

ggplot(tpc_and_environment, aes(x = t, y = performance, color = method)) +
  geom_line() +
  xlab("Patch temperature\n(degrees Celsius)") +
  ylab("Performance") +
  facet_grid(~Topt, labeller = label_both) +
  theme_dark(base_size = 14) +
  scale_color_viridis_d(name = "Model", option = 'magma')
ggsave("tpc_with_environment.png", width = 10, height = 3)

area <- function(a, b, Topt, epsilon){
  exp((a - Topt)/epsilon)*(a - 2*epsilon - Topt) + exp((b - Topt)/epsilon)*(-b + 2*epsilon + Topt)
}

epsilon_cost_df <- expand_grid(epsilon_baseline = 3,
                               Topt_baseline = 13,
                               epsilon = c(2.0, 3.0, 4.0, 5.0, 6.0),
                               Topt = c(Topt_baseline),
                               t = seq(-20, 28, length.out = 100),
                               Topt_cost = 0,
                               epsilon_cost = 0,
                               tmin = -50) |>
  mutate(x = (t-Topt)/epsilon,
         CTMax = Topt + epsilon,
         area = area(tmin, CTMax, Topt, epsilon),
         penalty = pmin(area(tmin, Topt_baseline+epsilon_baseline, Topt_baseline, epsilon_baseline)/area(tmin, CTMax, Topt, epsilon), 1.0),
         P1 = exp(x)*(1-x)*penalty,
         P = if_else(P1>0, P1, 0))
ggplot(epsilon_cost_df, aes(x = t, y = P, color = factor(round(epsilon, 2)))) +
  geom_line(linewidth = 0.5) +
  xlab("Temperature") +
  ylab("Survival Probability") +
  theme_bw(base_size = 10) +
  scale_color_discrete(name = expression(epsilon)) +
  scale_x_continuous(breaks = c(-10, 0, 10, 20)) +
  labs(title = expression("TPCs with constant area under the curve and baseline of" ~ epsilon == 3))
ggsave("tpc_constant_area.png", width = 10, height = 3)
