library(tidyverse)
library(glue)

##------------------------------------------------------
## Make files for QTL simulations for CDMetaPOP_slim ###
## -----------------------------------------------------

# Directory for all files
datadir = "./parameters/"
dir.create(file.path(datadir), showWarnings = FALSE)

## Test low and high dispersal
# All movement is through dispersal
# If an individual disperses they can't stay in their natal patch
# Test independent environment effect
parameter_combinations = expand_grid(dispersal_distances = c(0, 1, 5),
                                     environment = c('Y', 'N'))
distance_traveled_vals = parameter_combinations$dispersal_distances
environment_vals = parameter_combinations$environment

mcruns = 1 # Number of runs for each parameter

Popvars_names = paste0("PopVars", 1:nrow(parameter_combinations), ".csv") # File with population level parameters
Patchvars_names = paste0("PatchVars", 1:nrow(parameter_combinations), ".csv")


# Years to change parameters
# Temperature is constant for 0 through Tdelta, then increases by temp_change_amount
# for the next change_years, then constant at new warmer temperature until the end
Tdelta = 799
change_years = 50
cdclimgentime =  paste(c(0, Tdelta:(Tdelta + change_years)), collapse = "|") 
temp_change_amount = 2/change_years
#cdclimgentime = "0"
# Total number of generations to simulate
runtime = 1500
#runtime = 200
output_years = paste0(c(1, Tdelta, Tdelta + change_years, runtime - 1),collapse = "|") # Years to write out individuals 

# Total population size to simulate
Ntotal = 1000

## Making RunVars file
gridformat = 'genalex' # Output for alleles

runvars_df <- tibble(Popvars = Popvars_names,sizecontrol = 'N',mcruns,runtime,output_years,gridformat,cdclimgentime,
                     dispersal_distance = distance_traveled_vals, mutation_output_years = NA,
                     mutation_output_subpops = NA, mutation_origin_year = NA, environment = environment_vals)
## ClassVars
classvars_file <- "ClassVars.csv"

# The following are sizecontrol parameters. Since sizecontrol is not active,
# they are used to calculate
# age dependent parameters

# Growth for size control is controlled by temperature
# Parameters of Von B's equation for Westslope cutthroat trout
growth_Loo = 72.2 
growth_R0 = 0.19 # Taken from example
growth_temp_max = 10.1 # 
growth_temp_CV = 10.0
growth_temp_t0 = -4.997

# Mean offspring number depends on age or sizes
Egg_Mean_ans = 'exp' 
Egg_Mean_par1 = 9.557633 
Egg_Mean_par2 = 0.018109 

# Describe probability of maturity based on size. Only used if size control is active.
mature_eqn_slope = '0.06~0.13'
mature_eqn_int = '-8.09~-20.28'

# Loop through dispersal values
for (i in 1:nrow(parameter_combinations)){
  ## Making PopVars and PatchVars ##
  avg_dist_traveled <- distance_traveled_vals[i]
  popvars_filename <- Popvars_names[i]
  xyfilename = Patchvars_names[i]
  
  popmodel = 'packing' # Density dependent model, competition between all classes
  popmodel_par1 = -0.6821 # Parameter that defines the distribution of carrying capacities among classes,
  # -0.6821 is for westslope cuthtroat trout
  
  mate_cdmat = "prob_matrix_mate.csv" # Cost distance matrix for mating
  migrateout_cdmat = "prob_matrix_disperse.csv" # Cost distance matrix for migration out of patches
  migrateback_cdmat = migrateout_cdmat # Cost distance matrix for migration back to natal patches
  stray_cdmat = migrateout_cdmat # Cost distance matrix for individuals straying to other patches when returning
  disperseLocal_cdmat = migrateout_cdmat # Cost distance matrix for dispersal between patches
  
  # Mating model
  matemoveno = '9' # Mates are chosen according to mate_matrix probabilities
  
  # Reproduction parameters
  mature_default = 'age6~age6' # Default age of maturity for each sex
  
  # Fertility parameters
  offno = '2' # Poisson draw from mean offspring number
  Egg_Mortality = 0 # Population level egg mortality
  
  ## QTL variables
  genome = 'genome.csv'
  qtl_prop_genome = 1.0
  Topt_baseline = 10.0
  epsilon_baseline = 7.0
  Topt_VE = 1.0
  epsilon_VE = 0.1
  Topt_cost = 0
  epsilon_cost = 0
  Topt_pheno_eff = 'sample(c(-1,1), 1)*rexp(1, mu = 0.5)'
  epsilon_pheno_eff = '0'
  qtl_mutations_initial = 1000
  
  ## Making patchvars file
  # X and Y coordinates of patches
  npatches = 15
  min_temp = 5
  max_temp = 15
  patches <- tibble(PatchID = 0:(npatches - 1),
                    id = PatchID,
                    X = PatchID, Y = 0,
                    mean_temp = min_temp + PatchID*(max_temp - min_temp)/(npatches -1))
  X = patches$X
  Y = patches$Y
  PatchID = patches$PatchID
  
  ## Make probability matrix so patch 1 connected to patch 2, patch 2 to patch 3, etc.
  distance_matrix = as.matrix(dist(0:(npatches - 1), diag = TRUE, upper = TRUE))
  set.seed(17889)
  prob_matrix <- dnorm(distance_matrix, mean = avg_dist_traveled, sd =  mean(abs(rnorm(1000, 0, avg_dist_traveled))))
  prob_matrix <- 1 - distance_matrix/max(distance_matrix)
  write(prob_matrix, paste0(datadir, disperseLocal_cdmat), sep = ",", ncolumns = npatches)
  prob_matrix_mate = matrix(0, nrow = npatches, ncol = npatches)
  diag(prob_matrix_mate) = 1
  write(prob_matrix_mate, paste0(datadir, mate_cdmat), sep = ",", ncolumns = npatches)
  

  # Carrying capacity
  K = rep(round(Ntotal/npatches), npatches)
  KStDev = rep(0, npatches)
  N0 = K # Initial pop size
  NatalGrounds = rep(1, npatches) # Individuals can spawn here
  MigrationGrounds = rep(1, npatches) # Individuals can occupy this patch during migration
  MigrationOutProb = 0 # Migration probability?
  MigrationBackProb = 0 # Migration probability?
  StrayingProb = '0'
  DispersalProb = 0.1
  ClassVars = classvars_file
  # Extra mortality
  MortalityOut = 0
  MortalityBack = 0
  MortalityEggs = 0.68
  # Temperature values and number of grow days for temperature-dependent growth (in degrees celsius)
  GrowthTemperatureOut = 0
  GrowthTemperatureOutStDev = 0
  GrowDaysOut = 0
  GrowDaysOutStDev = 0
  
  # Increase temperature over time
  change_times = as.numeric(strsplit(cdclimgentime, split = "|", fixed = TRUE)[[1]])
  nclimgen = length(change_times)
  if(nclimgen > 1){
    GrowthTemperatureBack_current = patches$mean_temp
    GrowthTemperatureBack_matrix = GrowthTemperatureBack_current
    for(time in change_times[2:nclimgen]){
      GrowthTemperatureBack_current = GrowthTemperatureBack_current + temp_change_amount
      GrowthTemperatureBack_matrix = cbind(GrowthTemperatureBack_matrix, GrowthTemperatureBack_current)
    }
    GrowthTemperatureBack = apply(GrowthTemperatureBack_matrix, 1, paste, collapse = "|")
  }
  if(nclimgen == 1){GrowthTemperatureBack = patches$mean_temp}
  GrowthTemperatureBackStDev = 0
  GrowDaysBack = 365
  GrowDaysBackStDev = 0
  patchvars_df <- tibble(PatchID,X,Y,K,'K StDev' = KStDev,N0,
                         'Natal Grounds' = NatalGrounds,
                         'Migration Grounds' = MigrationGrounds,
                         'Class Vars' = ClassVars,
                         'Mortality Out %' = MortalityOut,
                         'Mortality Back' = MortalityBack,
                         'Mortality Eggs' = MortalityEggs,
                         'Migration Out Prob' = MigrationOutProb,
                         'Migration Back Prob' = MigrationBackProb,
                         'Straying Prob' = StrayingProb,
                         'Dispersal Prob' = DispersalProb,
                         GrowthTemperatureOut,GrowthTemperatureOutStDev,
                         GrowDaysOut,GrowDaysOutStDev,
                         GrowthTemperatureBack,GrowthTemperatureBackStDev,
                         GrowDaysBack,GrowDaysBackStDev)
  write_csv(patchvars_df, paste0(datadir, xyfilename))
  
  # Years to write out origin of mutations 
  mutation_output_years = paste0(0:(runtime - 1), collapse = "|")
  # Focal subpops for mutation origins
  mutation_output_subpops = paste(PatchID, collapse = "|")
  runvars_df$mutation_output_years[i] = mutation_output_years
  runvars_df$mutation_output_subpops[i] = mutation_output_subpops
  # Year when climate starts changing
  runvars_df$mutation_origin_year[i] = Tdelta
  popvars_df <- tibble(xyfilename,
                       mate_cdmat,matemoveno,migrateout_cdmat, migrateback_cdmat,
                       stray_cdmat,disperseLocal_cdmat,
                       mature_default,offno,
                       Egg_Mortality,
                       popmodel,popmodel_par1,
                       genome, 
                       qtl_prop_genome,
                       Topt_baseline,
                       epsilon_baseline,
                       Topt_VE,
                       epsilon_VE,
                       Topt_cost,
                       epsilon_cost,
                       Topt_pheno_eff,
                       epsilon_pheno_eff,
                       qtl_mutations_initial)
  # Parameters for independent temperature effect
  if(environment_vals[i] == 'Y'){
    popvars_df <- mutate(popvars_df,
                         independent_environment_mean = mean(patches$mean_temp),
                         independent_environment_sd = 5)
  }
  write_csv(popvars_df, paste0(datadir, popvars_filename))
}
write_csv(runvars_df, paste0(datadir, "RunVars.csv"))


## Make ClassVars
# These parameters come from Westslope cutthroat trout
classvars = tibble("Age class" = c(0, 1, 2, 3, 4, 5, 6, 7),
                   "Body Size Mean (mm)" = c(49, 82, 122, 157, 178, 188, 198, 217),
                   "Body Size Std (mm)" = c(3, 11, 17, 18, 17, 19, 20, 20),
                   "Distribution" = c(0.57, 0.17, 0.11, 0.065, 0.04, 0.023, 0.014, 0.008),
                   "Sex Ratio" = c(".50~.50", ".50~.50", ".50~.50", ".50~.50", ".50~.50", ".50~.50", 
                                   ".50~.50", ".50~.50"),
                   "Age Mortality Out %" = c(0, 0, 0, 0, 0, 0, 0, 95),
                   "Age Mortality Out StDev" = c(0, 0, 0, 0, 0, 0, 0, 0),
                   "Age Mortality Back %" = c(0, 0, 0, 0, 0, 0, 0, 0),
                   "Age Mortality Back StDev" = c(0, 0, 0, 0, 0, 0, 0, 0),
                   "Migration Out Prob" = c(1, 1, 1, 1, 1, 1, 1, 1),
                   "Migration Back Prob" = c(1, 1, 1, 1, 1, 1, 1, 1),
                   "Straying Prob" = c(1, 1, 1, 1, 1, 1, 1, 1),
                   "Dispersal Prob" = c(1, 1, 1, 1, 1, 1, 1, 1))

# For age control we also need to specify "Maturation" and "Fecundity Ind"
# To do this, convert size control parameters to age control parameters
# by computing the maturation
# and fecundity probabilities at the mean size for each age class
source("parameter_functions.R")
# Size control parameters
popvars <- tibble(growth_temp_max, growth_temp_CV, growth_R0, growth_Loo, growth_temp_t0,
                  mature_default, mature_eqn_int, mature_eqn_slope,
                  Egg_Mean_par1, Egg_Mean_par2) |>
  mutate(mature_age = gsub("age", "", mature_default))

classvars_new <- classvars |>
  mutate(Maturation_F = round(prob_mature_vector(`Body Size Mean (mm)`, 'F', `Age class`, popvars), 3),
         Maturation_M = round(prob_mature_vector(`Body Size Mean (mm)`, 'M', `Age class`, popvars), 3),
         Maturation = map2_vec(Maturation_F, Maturation_M, \(x,y) paste0(c(x, y), collapse = "~")),
         `Fecundity Ind` = mean_num_eggs(`Body Size Mean (mm)`, popvars)
  ) |> select(-Maturation_F, -Maturation_M)
write_csv(classvars_new, paste0(datadir, classvars_file)) 

library(readxl)

# Trout genome
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
# Remove sex chromosome (30) and compute recombination rates for smaller chromosomes
scale = 0.01
trout_genome <- trout_chroms |> mutate(recombination_rate = 3) |>
  mutate(Length = Length*scale, recombination_rate = recombination_rate/scale,
         mutation_rate = 1e-7) |>
  filter(Chromosome != 30)
# Write lengths to file
write_csv(trout_genome, paste0(datadir, genome))
