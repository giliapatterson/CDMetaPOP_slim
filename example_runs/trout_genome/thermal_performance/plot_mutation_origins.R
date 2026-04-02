library(tidyverse)

mutation_origins <- read_csv("../../slimgui_testing/mutation_origins.csv") |>
  group_by(subpop, tick) |>
  reframe(pick(everything()), Topt_prop = Topt_count/sum(Topt_count),
          neutral_prop = neutral_count/sum(neutral_count))

ggplot(mutation_origins, aes(x = tick, y = Topt_prop, color = factor(origin_subpop))) +
  geom_line() +
  geom_vline(xintercept = 800) +
  facet_wrap(~subpop, scale = 'free_y')
ggplot(mutation_origins, aes(x = tick, y = neutral_prop, color = factor(origin_subpop))) +
  geom_line() +
  geom_vline(xintercept = 800) +
  facet_wrap(~subpop, scale = 'free_y')
ggplot(mutation_origins, aes(x = tick, y = Topt_prop-neutral_prop, color = factor(origin_subpop))) +
  geom_line() +
  geom_vline(xintercept = 800) +
  facet_wrap(~subpop, scale = 'free_y')

ggplot(mutation_origins, aes(x = tick, y = Topt_mean_contribution, color = factor(origin_subpop))) +
  geom_line() +
  geom_vline(xintercept = 800) +
  facet_wrap(~subpop, scale = 'free_y')
ggplot(mutation_origins, aes(x = tick, y = Topt_VA, color = factor(origin_subpop))) +
  geom_line() +
  geom_vline(xintercept = 800) +
  facet_wrap(~subpop, scale = 'free_y')
