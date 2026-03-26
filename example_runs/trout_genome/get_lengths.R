library(tidyverse)
library(readxl)


# Read in trout linkage groups from Flores A-M, Christensen KA, Godin T, Palti Y, Campbell MR, Waldbieser GC, et al. The genome assembly of the westslope cutthroat trout, Oncorhynchus lewisi , reveals interspecific chromosomal rearrangements with the rainbow trout, Oncorhynchus mykiss. G3: Genes, Genomes, Genetics. 2025;15:jkaf064. https://doi.org/10.1093/g3journal/jkaf064.
trout_contigs <- read_excel("File_S2_G3-2025-405746.xlsx", skip = 1)
# Compute size of each arm for each linkage group
trout_arms <- trout_contigs |> mutate(length = nchar(Consensus)) |>
  group_by(LinkageGroup, Arms, Centromere) |>
  reframe(Length = sum(length), Centricity = unique(Centricity)) |>
  rename(Chromosome = LinkageGroup)
# Compute size of each linkage group
trout_chroms <- trout_arms |> group_by(Chromosome) |>
  reframe(Length = sum(Length))
# Remove sex chromosome (29) and compute recombination rates for smaller chromosomes
scale = 0.01
trout_genome <- trout_chroms |> mutate(recombination_rate = 3) |>
  mutate(Length = Length*scale, recombination_rate = recombination_rate/scale,
         mutation_rate = 1e-8) |>
  filter(Chromosome != 29)
# Write lengths to file
write_csv(trout_genome,"climate_change_McKenzie/trout_genome.csv")

# # Compute average recombination rate for non-centromeres
# overall_rr <- 3 #cM/mB
# non_centromere_prop = trout_arms |> group_by(Centromere) |>
#   reframe(Length = sum(Length)) |>
#   mutate(prop_length = Length/sum(Length)) |>
#   filter(!Centromere) |>
#   pull(prop_length)
# non_centromere_rr <- overall_rr/non_centromere_prop
# # Compute centromere distances
# cent_dist <- trout_arms 
# 
# make_recomb_map <- function(trout_arms, chromosome, avg_recom){
#   # Make recombination map for each Chromosome
#   chrom <- trout_arms |> filter(Chromosome == 1)
#   # p arm first
#   p <- chrom |> filter(Arms == 'p') |> pull(Length)
#   start_positions_p = 0:(p-1)
#   # recombination rate is 0 at centromere, increases away
#   slope = -avg_recom/(median(start_positions_p))
#   intercept = -(median(start_positions_p)*2)*slope
#   recomb_rates = (start_positions_p)*slope + intercept
# }
# # Make recombination map for each Chromosome
# chrom <- trout_arms |> filter(Chromosome == 1)
# # p arm first
# p <- chrom |> filter(Arms == 'p') |> pull(Length)
# start_positions = 0:(p-1)
# # recombination rate is 0 at centromere, increases away
# avg_recom = 3
# slope = -avg_recom/(median(start_positions))
# intercept = -(median(start_positions)*2)*slope
# recomb_rates = (start_positions)*slope + intercept
# plot(start_positions, recomb_rates)
# # then centromere
# c <- chrom |> filter(Arms == 'c') |> pull(Length)
# start_positions = p:(c-1)
# # recombination rate is 0 at centromere, increases away
# avg_recom = 0
# recomb_rates = rep(0, length(start_positions))
# plot(start_positions, recomb_rates)
# # then q
# q <- chrom |> filter(Arms == 'q') |> pull(Length)
# start_positions = c:(q-1)
# # recombination rate is 0 at centromere, increases away
# avg_recom = 3
# # start_positions*slope + intercept
# # (min(start_positions) + (max(start_positions) - min(start_positions))/2)*slope + intercept = avg_recom
# # min(start_positions)*slope + intercept = 0
# # (min(start_positions) + (max(start_positions) - min(start_positions))/2)*slope - min(start_positions)*slope = avg_recom
# # (min(start_positions) + max(start_positions)/2 - min(start_positions)/2 - min(start_positions))*slope = avg_recom
# # (max(start_positions)/2 - min(start_positions)/2)*slope = avg_recom
# slope = avg_recom/(max(start_positions)/2 - min(start_positions)/2)
# intercept = -min(start_positions)*slope
# recomb_rates = (start_positions)*slope + intercept
# plot(start_positions, recomb_rates)
