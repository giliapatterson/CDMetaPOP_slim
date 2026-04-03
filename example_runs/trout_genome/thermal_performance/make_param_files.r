library(tidyverse)
library(glue)

##------------------------------------------------------
## Make files for QTL simulations for CDMetaPOP_slim ###
## -----------------------------------------------------

# Directory for all files
datadir = "./"
dir.create(file.path(datadir), showWarnings = FALSE)

## Test low and high dispersal
# All movement is through dispersal
# If an individual disperses they can't stay in their natal patch
dispersal_prob_vals = c(0.01, 0.01)
# Test age and size control
sizecontrol_vals = c('Y', 'N')

mcruns = 1 # Number of runs for each parameter

Popvars_names = paste0("PopVars", 1:length(dispersal_prob_vals), ".csv") # File with population level parameters
Patchvars_names = paste0("PatchVars", 1:length(dispersal_prob_vals), ".csv")


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
Ntotal = 10000
#Ntotal = 1000
# Shorten patches for testing
set.seed(100389)
npatches <- 125
#npatches = 5
patch_is <- sort(sample(1:125, npatches, replace = FALSE))

## Making RunVars file
gridformat = 'genalex' # Output for alleles

runvars_df <- tibble(Popvars = Popvars_names,sizecontrol = sizecontrol_vals,mcruns,runtime,output_years,gridformat,cdclimgentime,
                     dispersal_prob = dispersal_prob_vals, mutation_output_years = NA,
                     mutation_output_subpops = NA)
## ClassVars
classvars_file <- "ClassVars.csv"

# The following are sizecontrol parameters. If size control is active, they are
# specified in PopVars. If sizecontrol is not active, they are used to calculate
# age dependent parameters

# Growth for size control is controlled by temperature
# Parameters of Von B's equation for Westslope cutthroat trout
growth_Loo = 72.2 
growth_R0 = 0.19 # Taken from example
growth_temp_max = 10.1 # 
growth_temp_CV = 0.468 # These parameters are also
growth_temp_t0 = -4.997 # taken from example

# Mean offspring number depends on age or sizes
Egg_Mean_ans = 'exp' 
Egg_Mean_par1 = 9.557633 
Egg_Mean_par2 = 0.018109 

# Describe probability of maturity based on size. Only used if size control is active.
mature_eqn_slope = '0.06~0.13'
mature_eqn_int = '-8.09~-20.28'

# Loop through dispersal values
for (i in 1:length(dispersal_prob_vals)){
  sizecontrol = sizecontrol_vals[i]
  
  ## Making PopVars and PatchVars ##
  dispersal_prob_val <- dispersal_prob_vals[i]
  popvars_filename <- Popvars_names[i]
  xyfilename = Patchvars_names[i]
  
  popmodel = 'packing' # Density dependent model, competition between all classes
  popmodel_par1 = -0.6821 # Parameter that defines the distribution of carrying capacities among classes,
  # -0.6821 is for westslope cuthtroat trout
  
  mate_cdmat = "prob_matrix_trout.csv" # Cost distance matrix for mating
  migrateout_cdmat = "prob_matrix_0_diagonal.csv" # Cost distance matrix for migration out of patches
  migrateback_cdmat = migrateout_cdmat # Cost distance matrix for migration back to natal patches
  stray_cdmat = migrateout_cdmat # Cost distance matrix for individuals straying to other patches when returning
  disperseLocal_cdmat = migrateout_cdmat # Cost distance matrix for dispersal between patches
  
  ## Read in distance matrix
  distance_matrix <- read_csv("../stream_network/trout_patches_small_distance_matrix.csv", col_names = FALSE)
  distance_matrix <- distance_matrix[patch_is, patch_is]
  # Convert to probability matrix
  prob_matrix <- 1 - distance_matrix/max(distance_matrix)
  write_csv(prob_matrix, paste0(datadir, mate_cdmat), col_names = FALSE, quote = "none")
  # Convert diagonal to 0s
  diag0 <- as.matrix(prob_matrix)
  diag(diag0) <- 0
  diag0 <- data.frame(diag0)
  write_csv(diag0, paste0(datadir, migrateout_cdmat), col_names = FALSE, quote = "none")
  
  # Mating model
  matemoveno = '9' # Mates are chosen according to mate_matrix probabilities
  
  # Reproduction parameters
  mature_default = 'age6~age6' # Default age of maturity for each sex
  
  # Fertility parameters
  offno = '2' # Poisson draw from mean offspring number
  Egg_Mortality = 0 # Population level egg mortality
  
  ## QTL variables
  genome = 'trout_genome.csv'
  qtl_prop_genome = 0.5
  Topt_baseline = 10.0
  epsilon_baseline = 2.0
  Topt_VE = 1.0
  epsilon_VE = 0.1
  Topt_cost = 0.01
  epsilon_cost = 0.05
  Topt_pheno_eff = 'sample(c(-1,1), 1)*rexp(1, mu = 0.1)'
  epsilon_pheno_eff = 'sample(c(-1,1), 1)*rexp(1, mu = 0.1)'
  
  ## Making patchvars file
  # Read in big patches and convert to small
  # Read in patches
  patches <- read_csv("../stream_network/trout_patches_small.csv")
  patches <- patches[patch_is,]
  ## Populate PatchVars file
  # X and Y coordinates of patches
  num_patches <- nrow(patches)
  X = patches$X
  Y = patches$Y
  PatchID = 1:num_patches # Unique patch ids
  patches$id = PatchID
  
  # Carrying capacity
  K = rep(round(Ntotal/num_patches), num_patches)
  KStDev = rep(0, num_patches)
  N0 = K # Initial pop size
  NatalGrounds = rep(1, num_patches) # Individuals can spawn here
  MigrationGrounds = rep(1, num_patches) # Individuals can occupy this patch during migration
  MigrationOutProb = 0 # Migration probability?
  MigrationBackProb = 0 # Migration probability?
  StrayingProb = '0'
  DispersalProb = dispersal_prob_val
  GenesInitialize = rep('random', num_patches) # How to initialize allele frequencies
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
  mutation_output_years = paste0(c(1, Tdelta, Tdelta + change_years, runtime - 1), collapse = "|")
  # Focal subpops for mutation origins
  # Pick subpops with initial temperatures evenly spaced
  nsubpops = 5
  quantiles = quantile(patches$mean_temp, probs = seq(0, 1, length.out = nsubpops))
  mutation_output_subpops = patches |> expand_grid(q_temp = quantiles) |>
    mutate(d_from_q = abs(mean_temp - q_temp)) |>
    group_by(id) |>
    filter(d_from_q == min(d_from_q)) |>
    group_by(q_temp) |>
    filter(d_from_q == min(d_from_q)) |>
    group_by(q_temp) |>
    summarize(id = sample(id, 1)) |>
    pull(id) |>
    unique() |>
    paste(collapse = "|")
  runvars_df$mutation_output_years[i] = mutation_output_years
  runvars_df$mutation_output_subpops[i] = mutation_output_subpops
    
  if(sizecontrol == 'Y'){
    popvars_df <- tibble(xyfilename,
                         mate_cdmat,matemoveno,migrateout_cdmat, migrateback_cdmat,
                         stray_cdmat,disperseLocal_cdmat,
                         mature_default,mature_eqn_slope,mature_eqn_int,offno,
                         Egg_Mean_ans,Egg_Mean_par1,Egg_Mean_par2,Egg_Mortality,
                         growth_Loo,growth_R0,growth_temp_max,growth_temp_CV,growth_temp_t0,
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
                         epsilon_pheno_eff)
  }
  if(sizecontrol == 'N'){
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
                         epsilon_pheno_eff)
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
