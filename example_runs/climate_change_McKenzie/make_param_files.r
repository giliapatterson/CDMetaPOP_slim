library(tidyverse)
library(terra)
library(sf)
library(lwgeom)
library(sfnetworks)
library(tidygraph)
library(igraph)
library(units)
library(lubridate)

##-----------------------------------------------
## Make files for CDMetaPOP trout simulations ###
## ----------------------------------------------

# Directory for all files
datadir = ""
# Total population size to simulate
#Ntotal = 10000
Ntotal = 1000
# Shorten patches for testing
set.seed(100389)
npatches <- 5
patch_is <- sort(sample(1:125, npatches, replace = FALSE))

# All movement is through dispersal
# If an individual disperses they can't stay in their natal patch
dispersal_prob_val <- 0.1
# Age or size control?
sizecontrol = 'Y' # N means that age-based mortality controls population size

## Making RunVars file
# Helper variables
runvars_filename_base = "RunVars"
#burnin = 10 # Number of years of burn in
  
# Variables that need to be specified in the file
Popvars = paste0("PopVars", dispersal_prob_val, ".csv") # File with population level parameters
constMortans = 1 # 1 means that all events that lead to mortality are mutually exclusive (maybe? idk)
runtime = 200 # Number of generations to simulate
output_years = paste0(c(1, runtime - 1),collapse = "|") # Years to write out individuals 
#gridformat = 'genalex' # Output for alleles
gridformat = 'cdpop'
gridsampling = 'N' # Output all individuals
summaryOutput = 'N' # Whether to output summary files
# Years to change parameters
# Temperature is constant for 0 through Tdelta, then increases by 0.2 degrees
# for the next change_years, then constant at new warmer temperature until the end
Tdelta = 49
change_years = 10
cdclimgentime =  paste(c(0, Tdelta:(Tdelta + change_years)), collapse = "|") 
#cdclimgentime = "0"
startcomp = 0 # When to start competition (don't need yet, so leave at 0)
implementcomp = 'Back' # Don't need to implement competition yet so set to default

# Write to file
num_mcruns = 1 # Number of runs for each set of parameters
runvars_df <- tibble(Popvars,sizecontrol,constMortans,mcruns = 1,runtime,output_years,gridformat,gridsampling,summaryOutput,cdclimgentime,startcomp,implementcomp)
write_csv(runvars_df, paste0(datadir, runvars_filename_base, ".csv"))

## ClassVars
classvars_file <- "ClassVars.csv"

# Make PopVars and PatchVars
dispersal_prob <- dispersal_prob_val
popvars_filename <- Popvars

## Making PopVars file
xyfilename = paste0("PatchVars", dispersal_prob, ".csv") # Patches
popmodel = 'packing' # Density dependent model, competition between all classes
popmodel_par1 = -0.6821 # Parameter that defines the distribution of carrying capacities among classes,
# -0.6821 is for westslope cuthtroat trout
correlation_matrix = 'N' # Correlation of variables between patches, N means none, can specify file
subpopmort_file = 'N' # Mortality rates for individuals from one patch migrating into another
egg_delay = 0 # Number of years between mating and emergence
egg_add = 'mating' # Idk

# Variables that can change with climate (CDClimGen)
mate_cdmat = "prob_matrix_trout.csv" # Cost distance matrix for mating
migrateout_cdmat = "prob_matrix_0_diagonal.csv" # Cost distance matrix for migration out of patches
migrateback_cdmat = migrateout_cdmat # Cost distance matrix for migration back to natal patches
stray_cdmat = migrateout_cdmat # Cost distance matrix for individuals straying to other patches when returning
disperseLocal_cdmat = migrateout_cdmat # Cost distance matrix for dispersal between patches

## Read in big probability matrices and convert to small
prob_matrix_big <- read_csv("prob_matrix_trout_big.csv", col_names = FALSE)
prob_matrix <- prob_matrix_big[patch_is, patch_is]
write_csv(prob_matrix, paste0(datadir, mate_cdmat), col_names = FALSE, quote = "none")
prob_matrix_0_diagonal_big <- read_csv("prob_matrix_0_diagonal_big.csv", col_names = FALSE)
prob_matrix_0_diagonal <- prob_matrix_0_diagonal_big[patch_is, patch_is]
write_csv(prob_matrix_0_diagonal, paste0(datadir, migrateout_cdmat), col_names = FALSE, quote = "none")

# Parameters that control how costs are converted to probabilities
matemoveno = '9' # Probability of movement = 1 - (cd - min)/threshold (I think, should check)
matemovethresh = 'max' # Threshold is maximum cd value
migratemoveOutno = '9'
migratemoveOutthresh = 'max'
migratemoveBackno = '9'
migratemoveBackthresh = 'max'
StrayBackno = '6' # No straying
StrayBackthresh = 0 # Doesn't matter since there is no straying
HomeAttempt = 'mortality' # What happens if an individual can't make it home. Doesn't matter because there is no migration
disperseLocalno = '9' # Probability of movement = 1 - (cd - min)/threshold (I think, should check)
disperseLocalthresh = 'max' # Threshold is maximum cd value

# These params don't matter for numbers '6' or '8', so set to 0
matemoveparA = 0
matemoveparB = 0
matemoveparC = 0
migratemoveOutparA = matemoveparA
migratemoveOutparB = matemoveparB
migratemoveOutparC = matemoveparB
migratemoveBackparA = matemoveparA
migratemoveBackparB = matemoveparB
migratemoveBackparC = matemoveparB
StrayBackparA = matemoveparA
StrayBackparB = matemoveparB
StrayBackparC = matemoveparB
disperseLocalparA = matemoveparA
disperseLocalparB = matemoveparB
disperseLocalparC = matemoveparB

# Reproduction parameters
sex_chromo = 2 # Number of sex chromosomes
sexans = 'Y' # Sexual reproduction
selfans = 'N' # No selfing
Freplace = 'N' # Females can mate with multiple males
Mreplace = 'Y' # Males can mate with multiple females
AssortativeMate_Model = '1' # Random mating. I think this doesn't matter since there are no hybrids
AssortativeMate_Factor = '1' # Random mating. Not sure if this matters either
mature_default = 'age6~age6' # Size or age of maturity for each sex
# Describe probability of maturity based on size. Only used if size control is active.
# Idk what it means for size control to be active. In competition? Or in maturity?
mature_eqn_slope = '0.06~0.13'
mature_eqn_int = '-8.09~-20.28'

# Fertility parameters
offno = '2' # Poisson draw from mean offspring number
# Mean offspring number depends on age and is set in the classvars file
# Only matters if age control is active
offans_InheritClassVars = 'random' # I think only matters if multiple ClassVars files
# I only have one file so doesn't matter
equalClutchSize = 'N' # Does each mate pair have equal clutch size? Yes means that females will
# always have the same number of offspring regardless of number of mates
Egg_Freq_Mean = 1  # Mean and standard deviation of number of egg laying events
Egg_Freq_StDev = 0 # per year
# Next three only used if size control is specified
Egg_Mean_ans = 'exp' 
Egg_Mean_par1 = 9.557633 
Egg_Mean_par2 = 0.018109 
Egg_Mortality = 0 # Population level egg mortality. Taken from example file
Egg_Mortality_StDev = 0 # And standard deviation. Also from example file
Egg_FemaleProb = 0.5 # Probability an offspring will be female

# Genetic parameters
startGenes = 0 # What time to start genetics
loci = 50 # Number of loci
#loci = 5
alleles = 2 # Number of alleles per locus - 1
muterate = 0.0001 # Mutation rate
mutationtype = 'random' # Mutation model (random means all alleles can mutate to other alleles)
mtdna = 'N' # Don't model maternally inherited genes

# Selection
#cdevolveans = '1' # Selection on one locus
cdevolveans = 'N'
startSelection = 0 # When to start selection
implementSelection = 'Eggs' # Selection implemented only through egg mortality
betaFile_selection = 'N' # No polygenic selection

# Phenotypic plasticity
plasticgeneans = 'N' # No plasticity
# Rest of plasticity parameters don't matter
plasticSignalResponse = 0
plasticBehavioralResponse = 0
startPlasticgene = 0
implementPlasticgene = 'back'

# Disease
implement_disease = 'N' # No disease

# Growth
growth_option = 'temperature' # Von Bert;alskdfj; growth equation that also is influenced by temperature
# I don't understand this
growth_Loo = 72.2 # Parameters of Von B's equation
growth_R0 = 0.19 # Taken from example
growth_temp_max = 10.12 # These parameters 
growth_temp_CV = 0.468 # are also
growth_temp_t0 = -4.997 # taken from example

## QTL variables
genome_length = 10000
qtl_prop_genome = 0.5
qtl_pheno_eff = 'rnorm(1, 0, 5.0)'
qtl_env_variable = 'GrowthTemperatureBack'
# Mutation rate, recombination rate, Ve, and fitness function SD
# can change over time
change_times = as.numeric(strsplit(cdclimgentime, split = "|", fixed = TRUE)[[1]])
nclimgen = length(change_times)
qtl_muterate = paste(c(1e-4, rep(1e-8, nclimgen - 1)), collapse = "|")
qtl_recrate = paste(c(1e-4, rep(1e-8, nclimgen - 1)), collapse = "|")
qtl_ve = paste(c(0.5, rep(1.0, nclimgen - 1)), collapse = "|")
qtl_fit_sd = paste(c(10.0, rep(5.0, nclimgen - 1)), collapse = "|")

popvars_df <- tibble(xyfilename,mate_cdmat,matemoveno,matemoveparA,
                     matemoveparB,matemoveparC,matemovethresh,
                     migrateout_cdmat,migratemoveOutno,migratemoveOutparA,
                     migratemoveOutparB,migratemoveOutparC,migratemoveOutthresh,
                     migrateback_cdmat,migratemoveBackno,migratemoveBackparA,
                     migratemoveBackparB,migratemoveBackparC,migratemoveBackthresh,
                     stray_cdmat,StrayBackno,StrayBackparA,StrayBackparB,
                     StrayBackparC,StrayBackthresh,disperseLocal_cdmat,disperseLocalno,
                     disperseLocalparA,disperseLocalparB,disperseLocalparC,
                     disperseLocalthresh,HomeAttempt,sex_chromo,sexans,selfans,
                     Freplace,Mreplace,AssortativeMate_Model,AssortativeMate_Factor,
                     mature_default,mature_eqn_slope,mature_eqn_int,offno,
                     offans_InheritClassVars,equalClutchSize,Egg_Freq_Mean,Egg_Freq_StDev,
                     Egg_Mean_ans,Egg_Mean_par1,Egg_Mean_par2,Egg_Mortality,
                     Egg_Mortality_StDev,Egg_FemaleProb,startGenes,loci,alleles,
                     muterate,mutationtype,mtdna,cdevolveans,startSelection,
                     implementSelection,betaFile_selection,plasticgeneans,
                     plasticSignalResponse,plasticBehavioralResponse,startPlasticgene,
                     implementPlasticgene,growth_option,growth_Loo,growth_R0,
                     growth_temp_max,growth_temp_CV,growth_temp_t0,popmodel,
                     popmodel_par1,correlation_matrix,subpopmort_file,egg_delay,
                     egg_add,implement_disease,
                     genome_length, qtl_prop_genome, qtl_pheno_eff, qtl_env_variable,
                     qtl_muterate, qtl_recrate, qtl_ve, qtl_fit_sd)
write_csv(popvars_df, paste0(datadir, popvars_filename))

## Making patchvars file
# Read in big patches and convert to small
patches <- read_csv("PatchVars_big.csv")
patches <- patches[patch_is,]
## Populate PatchVars file
# X and Y coordinates of patches
num_patches <- nrow(patches)
X = patches$X
Y = patches$Y
PatchID = 1:num_patches # Unique patch ids
SubpatchNO = rep(1, num_patches) # I have no idea what this is

# Carrying capacity
K = rep(round(Ntotal/num_patches), num_patches)
KStDev = rep(10, num_patches)
N0 = K # Initial pop size
NatalGrounds = rep(1, num_patches) # Individuals can spawn here
MigrationGrounds = rep(1, num_patches) # Individuals can occupy this patch during migration
MigrationOutProb = 0 # Migration probability?
MigrationBackProb = 0 # Migration probability?
SetMigration = 'N' # Whether individuals that become migrants stay migrants
StrayingProb = '0'
DispersalProb = dispersal_prob
GenesInitialize = rep('random', num_patches) # How to initialize allele frequencies
ClassVars = classvars_file
# Extra mortality
MortalityOut = 0
MortalityBack = 0
MortalityEggs = 0.68
MortalityOutStDev = 0
MortalityBackStDev = 0
MortalityEggsStDev = 0
# Temperature values and number of grow days for temperature-dependent growth (in degrees celsius)
GrowthTemperatureOut = 0
GrowthTemperatureOutStDev = 0
GrowDaysOut = 0
GrowDaysOutStDev = 0

# Increase temperature over time
change_times = as.numeric(strsplit(cdclimgentime, split = "|", fixed = TRUE)[[1]])
nclimgen = length(change_times)
if(nclimgen > 1){
  change_amount = 0.5
  GrowthTemperatureBack_current = patches$GrowthTemperatureBack
  GrowthTemperatureBack_matrix = GrowthTemperatureBack_current
  for(time in change_times[2:nclimgen]){
    GrowthTemperatureBack_current = GrowthTemperatureBack_current + change_amount
    GrowthTemperatureBack_matrix = cbind(GrowthTemperatureBack_matrix, GrowthTemperatureBack_current)
  }
  GrowthTemperatureBack = apply(GrowthTemperatureBack_matrix, 1, paste, collapse = "|")
}
if(nclimgen == 1){GrowthTemperatureBack = patches$GrowthTemperatureBack}
GrowthTemperatureBackStDev = 0
GrowDaysBack = 365
GrowDaysBackStDev = 0
# Capture probs
CaptureProbabilityOut = 'N'
CaptureProbabilityBack = 'N'
# Habitat quality, which doesn't matter because there is no plasticity
HabitatOut = 1
HabitatBack = 1
# No idea
comp_coef = 0
# Fitness values
Fitness_AA = 0
Fitness_Aa = 0
Fitness_aa = 0
Fitness_BB = 0
Fitness_Bb = 0
Fitness_bb = 0
Fitness_AABB = 0
Fitness_AaBB = 0
Fitness_aaBB = 0
Fitness_AABb = 0
Fitness_AaBb = 0
Fitness_aaBb = 0
Fitness_AAbb = 0
Fitness_Aabb = 0
Fitness_aabb = 0
patchvars_df <- tibble(PatchID,X,Y,SubpatchNO,K,'K StDev' = KStDev,N0,
                       'Natal Grounds' = NatalGrounds,
                       'Migration Grounds' = MigrationGrounds,
                       'Genes Initialize' = GenesInitialize,
                       'Class Vars' = ClassVars,
                       'Mortality Out %' = MortalityOut,
                       'Mortality Out StDev' = MortalityOutStDev,
                       'Mortality Back' = MortalityBack,
                       'Mortality Back StDev' = MortalityBackStDev,
                       'Mortality Eggs' = MortalityEggs,
                       'Mortality Eggs StDev' = MortalityEggsStDev,
                       'Migration Out Prob' = MigrationOutProb,
                       'Set Migration' = SetMigration,
                       'Migration Back Prob' = MigrationBackProb,
                       'Straying Prob' = StrayingProb,
                       'Dispersal Prob' = DispersalProb,
                       GrowthTemperatureOut,GrowthTemperatureOutStDev,
                       GrowDaysOut,GrowDaysOutStDev,
                       GrowthTemperatureBack,GrowthTemperatureBackStDev,
                       GrowDaysBack,GrowDaysBackStDev,
                       'Capture Probability Out' = CaptureProbabilityOut,
                       'Capture Probability Back' = CaptureProbabilityBack,
                       HabitatOut, HabitatBack,
                       Fitness_AA,Fitness_Aa,Fitness_aa,
                       Fitness_BB,Fitness_Bb,Fitness_bb,
                       Fitness_AABB,Fitness_AaBB,Fitness_aaBB,
                       Fitness_AABb,Fitness_AaBb,Fitness_aaBb,
                       Fitness_AAbb,Fitness_Aabb,Fitness_aabb,
                       comp_coef)
write_csv(patchvars_df, paste0(datadir, xyfilename))
