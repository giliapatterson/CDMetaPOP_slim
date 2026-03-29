library(tidyverse)

##-----------------------------------------------
## Make files for QTL simulations for CDMetaPOP_slim ###
## ----------------------------------------------

## Test low and high dispersal
# All movement is through dispersal
# If an individual disperses they can't stay in their natal patch
dispersal_prob_vals = c(0.01, 0.1)

# Age or size control?
sizecontrol = 'N' # N means that age-based mortality controls population size

# Years to change parameters
# Temperature is constant for 0 through Tdelta, then increases by temp_change_amount
# for the next change_years, then constant at new warmer temperature until the end
Tdelta = 799
change_years = 50
cdclimgentime =  paste(c(0, Tdelta:(Tdelta + change_years)), collapse = "|") 
temp_change_amount = 0.2
#cdclimgentime = "0"
# Total number of generations to simulate
runtime = 1500
#runtime = 200
output_years = paste0(c(1, Tdelta, Tdelta + change_years, runtime - 1),collapse = "|") # Years to write out individuals 

# Directory for all files
datadir = "climate_change_McKenzie/"
# Total population size to simulate
#Ntotal = 10000
Ntotal = 1000
# Shorten patches for testing
set.seed(100389)
npatches <- 5
patch_is <- sort(sample(1:125, npatches, replace = FALSE))

## Making RunVars file
Popvars_names = paste0("PopVars", dispersal_prob_vals, ".csv") # File with population level parameters

gridformat = 'genalex' # Output for alleles

# Write to file
mcruns = 1 # Number of runs
runvars_df <- tibble(Popvars = Popvars_names,sizecontrol,mcruns,runtime,output_years,gridformat,cdclimgentime)
write_csv(runvars_df, paste0(datadir, "RunVars.csv"))

# Loop through dispersal values
for (i in 1:length(dispersal_prob_vals)){
  ## Making PopVars and PatchVars ##
  dispersal_prob_val <- dispersal_prob_vals[i]
  popvars_filename <- Popvars_names[i]
  
  ## ClassVars
  classvars_file <- "ClassVars.csv"

  xyfilename = paste0("PatchVars", dispersal_prob_val, ".csv") # Patches
  popmodel = 'packing' # Density dependent model, competition between all classes
  popmodel_par1 = -0.6821 # Parameter that defines the distribution of carrying capacities among classes,
  # -0.6821 is for westslope cuthtroat trout
  
  mate_cdmat = "prob_matrix_trout.csv" # Cost distance matrix for mating
  migrateout_cdmat = "prob_matrix_0_diagonal.csv" # Cost distance matrix for migration out of patches
  migrateback_cdmat = migrateout_cdmat # Cost distance matrix for migration back to natal patches
  stray_cdmat = migrateout_cdmat # Cost distance matrix for individuals straying to other patches when returning
  disperseLocal_cdmat = migrateout_cdmat # Cost distance matrix for dispersal between patches
  
  ## Read in big probability matrices and convert to small
  prob_matrix_big <- read_csv(paste0(datadir, "prob_matrix_trout_big.csv"), col_names = FALSE)
  prob_matrix <- prob_matrix_big[patch_is, patch_is]
  write_csv(prob_matrix, paste0(datadir, mate_cdmat), col_names = FALSE, quote = "none")
  prob_matrix_0_diagonal_big <- read_csv(paste0(datadir, "prob_matrix_0_diagonal_big.csv"), col_names = FALSE)
  prob_matrix_0_diagonal <- prob_matrix_0_diagonal_big[patch_is, patch_is]
  write_csv(prob_matrix_0_diagonal, paste0(datadir, migrateout_cdmat), col_names = FALSE, quote = "none")
  
  # Mating model
  matemoveno = '9' # Mates are chosen according to mate_matrix probabilities
  
  # Reproduction parameters
  mature_default = 'age6~age6' # Default age of maturity for each sex
  # Describe probability of maturity based on size. Only used if size control is active.
  if(sizecontrol == 'Y'){
    mature_eqn_slope = '0.06~0.13'
    mature_eqn_int = '-8.09~-20.28'
  }
  
  # Fertility parameters
  offno = '2' # Poisson draw from mean offspring number
  # Mean offspring number depends on age or sizes
  # Next three only used if size control is specified
  if(sizecontrol == 'Y'){
    Egg_Mean_ans = 'exp' 
    Egg_Mean_par1 = 9.557633 
    Egg_Mean_par2 = 0.018109 
  }
  Egg_Mortality = 0 # Population level egg mortality
  
  # Alleles from original CDMetaPOP
  startGenes = 0 # What time to start genetics
  loci = 50 # Number of loci
  alleles = 2 # Number of alleles per locus - 1
  muterate = 0.0001 # Mutation rate
  mutationtype = 'random' # Mutation model (random means all alleles can mutate to other alleles)
  
  # Growth for size control is controlled by temperature
  if(sizecontrol == 'Y'){
    # Parameters of Von B's equation for Westslope cutthroat trout
    growth_Loo = 72.2 
    growth_R0 = 0.19 # Taken from example
    growth_temp_max = 10.12 # These parameters 
    growth_temp_CV = 0.468 # are also
    growth_temp_t0 = -4.997 # taken from example
  }
  
  ## QTL variables
  genome_length = 100000
  qtl_prop_genome = 0.5
  qtl_pheno_eff = 'rnorm(1, 0, 5.0)'
  qtl_env_variable = 'GrowthTemperatureBack'
  # Mutation rate, recombination rate, Ve, and fitness function SD
  # can change over time
  qtl_muterate = 1e-8
  qtl_recrate = 1e-8
  qtl_ve = 1.0
  qtl_fit_sd = 1.0
  
  if(sizecontrol == 'Y'){
    popvars_df <- tibble(xyfilename,
                         mate_cdmat,matemoveno,migrateout_cdmat, migrateback_cdmat,
                        stray_cdmat,disperseLocal_cdmat,
                        mature_default,mature_eqn_slope,mature_eqn_int,offno,
                        Egg_Mean_ans,Egg_Mean_par1,Egg_Mean_par2,Egg_Mortality,
                        startGenes,loci,alleles,muterate,mutationtype,
                        growth_Loo,growth_R0,growth_temp_max,growth_temp_CV,growth_temp_t0,
                        popmodel,popmodel_par1,
                        genome_length, qtl_prop_genome, qtl_pheno_eff, qtl_env_variable,
                        qtl_muterate, qtl_recrate, qtl_ve, qtl_fit_sd)
  }
  if(sizecontrol == 'N'){
    popvars_df <- tibble(xyfilename,
                         mate_cdmat,matemoveno,migrateout_cdmat, migrateback_cdmat,
                         stray_cdmat,disperseLocal_cdmat,
                         mature_default,offno,
                         Egg_Mortality,
                         startGenes,loci,alleles,muterate,mutationtype,
                         popmodel,popmodel_par1,
                         genome_length, qtl_prop_genome, qtl_pheno_eff, qtl_env_variable,
                         qtl_muterate, qtl_recrate, qtl_ve, qtl_fit_sd)
  }
  write_csv(popvars_df, paste0(datadir, popvars_filename))
  
  ## Making patchvars file
  # Read in big patches and convert to small
  patches <- read_csv(paste0(datadir, "PatchVars_big.csv"))
  patches <- patches[patch_is,]
  ## Populate PatchVars file
  # X and Y coordinates of patches
  num_patches <- nrow(patches)
  X = patches$X
  Y = patches$Y
  PatchID = 1:num_patches # Unique patch ids
  
  # Carrying capacity
  K = rep(round(Ntotal/num_patches), num_patches)
  KStDev = rep(10, num_patches)
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
    GrowthTemperatureBack_current = patches$GrowthTemperatureBack
    GrowthTemperatureBack_matrix = GrowthTemperatureBack_current
    for(time in change_times[2:nclimgen]){
      GrowthTemperatureBack_current = GrowthTemperatureBack_current + temp_change_amount
      GrowthTemperatureBack_matrix = cbind(GrowthTemperatureBack_matrix, GrowthTemperatureBack_current)
    }
    GrowthTemperatureBack = apply(GrowthTemperatureBack_matrix, 1, paste, collapse = "|")
  }
  if(nclimgen == 1){GrowthTemperatureBack = patches$GrowthTemperatureBack}
  GrowthTemperatureBackStDev = 0
  GrowDaysBack = 365
  GrowDaysBackStDev = 0
  patchvars_df <- tibble(PatchID,X,Y,K,'K StDev' = KStDev,N0,
                         'Natal Grounds' = NatalGrounds,
                         'Migration Grounds' = MigrationGrounds,
                         'Genes Initialize' = GenesInitialize,
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
}
