suppressPackageStartupMessages(library(tidyverse))
library(argparse)
library(glue)

thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}
script_directory <- dirname(thisFile())

source(paste0(script_directory,"/file_processing_functions.r"))
# source("file_processing_functions.r")
parser <- ArgumentParser()

parser$add_argument(
    "-d",
    "--parameter_directory",
    default = '../../../trout-and-frogs/trout_sim/trout_sim_parameters/',
    help = "Directory containing CDMetaPOP input files"
)

parser$add_argument(
    "-r",
    "--runvars_file_name",
    default = '../../../trout-and-frogs/trout_sim/trout_sim_parameters/RunVars.csv',
    help = "Name of RunVars file for CDMetaPOP"
)

parser$add_argument(
    "-o",
    "--output_directory",
    default = '../../../trout-and-frogs/trout_sim/trout_sim_parameters/slim_test_params/',
    help = "Directory for SLiM input files"
)

args <- parser$parse_args()
param_directory = args$parameter_directory # Directory containing CDMetaPOP parameter files
output_overall = args$output_directory
dir.create(file.path(output_overall), showWarnings = FALSE)

### Process RunVars ###
runvars_all <- read_csv(args$runvars_file_name, show_col_types = FALSE) 
nruns <- nrow(runvars_all)

for(run in 1:nruns){
  print(glue('Processing run {run} of {nruns}'))
  runvars <- runvars_all |> slice(run)
  output_directory = paste0(output_overall, "run", run, "/") # Directory to output SLiM parameter files
  dir.create(file.path(output_directory), showWarnings = FALSE) # Create output directory if it doesn't already exist

  # Does this run change the climate?
  climchangeyears <- as.numeric(str_split_1(as.character(pull(runvars, cdclimgentime)), fixed("|"))) + 1
  climate_change <- length(climchangeyears) > 1
  
  # SLiM code doesn't currently support all options of CDMetaPOP
  # Extract the options used by SLiM and add 1 to variables containing years to be consistent with SLiM
  add_one <- function(year_string){
    return(paste(as.numeric(str_split_1(year_string, fixed("|"))) + 1, collapse = "|"))
  }
  runvars <- mutate(runvars, output_years = add_one(as.character(output_years)),
                    cdclimgentime = add_one(as.character(cdclimgentime)))
  if(runvars$gridformat != 'genalex'){
    print(glue("gridformat {runvars$gridformat} is not supported. Using 'genalex' grid format."))
  }

  if(climate_change){runvars_used <- c("Popvars", "sizecontrol","runtime", "output_years", "cdclimgentime")}
  if(!climate_change){runvars_used <- c("Popvars", "sizecontrol", "runtime", "output_years")}
  runvars <- select(runvars, all_of(runvars_used))
  
  # Change entry for Popvars to match the Popvars file used for SLiM
  popvars_file_out <- paste0(output_directory, "PopVars_slim.csv")
  runvars_out <- runvars |>
    mutate(Popvars = popvars_file_out)
  
  # Output runvars
  write_csv(runvars_out, paste0(output_directory, "RunVars_slim.csv"))
  
  
  ### Process PopVars ###
  # Read in CDMetaPOP popvars
  popvars <- read_csv(paste0(param_directory, pull(runvars, Popvars)), show_col_types = FALSE)
  
  ## 1. Process patchvars ##
  # Read in CDMetaPOP patchvars
  patchvars <- read_csv(paste0(param_directory, pull(popvars, xyfilename)), show_col_types = FALSE)
  # Output file for patchvars for slim
  patchvars_file_out <- paste0(output_directory, "PatchVars_slim.csv")
  # Convert patchID to be 0 indexed. This is necessary for SLiM
  patchvars <- mutate(patchvars, PatchID = row_number() - 1)
  
  ## 1. (a) Process Genes ##
  
  # First check if `Genes Initialize` is random or random_var
  # random_var is unsupported
  genes_initialize = pull(patchvars, `Genes Initialize`)
  if (length(unique(genes_initialize)) == 1 & genes_initialize[1] == "random"){
    print(glue("Using random gene initialization."))
    genes_method = "random"
  }
  if (length(unique(genes_initialize)) == 1 & genes_initialize[1] == "random_var"){
    print(glue("Option random_var not supported for initializing genes, using random instead."))
    genes_method = "random"
  }
  
  if (length(unique(genes_initialize)) > 1 | !(genes_initialize[1] %in% c("random", "random_var"))){
    print(glue("Using gene initialization from file."))
    genes_method = "file"
  }
  
  ## 1. (a) (i)
  ## For method "file"
  ## Merge all gene files and assign a position to each locus.##
  ## This step is necessary to ensure that loci with the same name get assigned to 
  ## the same position in the genome in SLiM
  if(genes_method == "file"){
    # Read in all gene files and merge
    genes_list <- map2(paste0(param_directory, pull(patchvars, `Genes Initialize`)),
                       patchvars$PatchID,
                       split_genes)
    all_genes <- reduce(genes_list, merge_genes) |>
      rename("Frequency_0" = Frequency, "PatchID_0" = PatchID)
    all_genes <- mutate(all_genes, position = match(locus, unique(all_genes$locus)) - 1)
    
    ## 1. (a) (ii) Split gene files back into patches and write to file##
    # Folder for storing gene files for SLiM
    genes_folder = paste0(output_directory, "genes/")
    dir.create(file.path(genes_folder), showWarnings = FALSE)
    # Create a new column in patchvars with the names of the gene files for slim
    patchvars <- mutate(patchvars, genes_file_slim = paste0(genes_folder, "genes_", PatchID, ".csv"))
    # Split all genes into patches and write to appropriate files
    map2(patchvars$genes_file_slim, patchvars$PatchID, patch_genes, all_genes = all_genes)
  }
  
  ## 1. (b) Process classvars
  # CDMetaPOP classvars
  classvars = read_csv(paste0(param_directory, pull(patchvars, `Class Vars`)[1]), show_col_types = FALSE)
  # Remove unused columns
  if(runvars$sizecontrol == 'Y'){
    classvars_used = c("Age class", "Body Size Mean (mm)", "Body Size Std (mm)", "Distribution",
                       "Age Mortality Out %", "Age Mortality Back %", "Migration Out Prob",
                       "Migration Back Prob", "Straying Prob", "Dispersal Prob")
  }
  else{
    classvars <- mutate(classvars, Maturation_F = as.numeric(str_split_i(Maturation, "~", 1)),
                        Maturation_M = as.numeric(str_split_i(Maturation, "~", 2)))
    classvars_used = c("Age class", "Distribution",
                       "Age Mortality Out %", "Age Mortality Back %", "Migration Out Prob",
                       "Migration Back Prob", "Straying Prob", "Dispersal Prob",
                       "Maturation_F", "Maturation_M", "Fecundity Ind")
  }
  classvars_out <- classvars |> select(all_of(classvars_used))
  
  # Write to file
  classvars_out_file = paste0(output_directory, "classvars.csv")
  write_csv(classvars_out, classvars_out_file)
  
  # Update file in patchvars
  patchvars <- mutate(patchvars, classvars = classvars_out_file)
  
  ## 1. (c) Remove variables that aren't used by SLiM and write to file ##
  patchvars_used <- c("PatchID", "X", "Y", "K", "K StDev", "N0", "Mortality Eggs",
                      "Migration Out Prob", "Migration Back Prob", "Straying Prob",
                      "Dispersal Prob", "GrowthTemperatureOut", "GrowthTemperatureOutStDev",
                      "GrowDaysOut", "GrowDaysOutStDev", "GrowthTemperatureBack", "GrowthTemperatureBackStDev",
                      "GrowDaysBack", "GrowDaysBackStDev", "classvars", "Natal Grounds", "Migration Grounds")
  if(genes_method == "file"){
     patchvars_used = c(patchvars_used, "genes_file_slim")
  }
  
  patchvars_out <- select(patchvars, all_of(patchvars_used))                    
  write_csv(patchvars_out, patchvars_file_out)
  
  # Update patchvars file name in popvars
  popvars <- popvars |> mutate(xyfilename = patchvars_file_out)
  
  ## 2. Process matrices ##
  # If the climate changes over time, popvars has multiple rows, one for each
  # time step
  if(climate_change){
    # Make popvars multiple rows
    popvars_new <- bind_rows(rep(list(popvars), length(climchangeyears))) |>
      mutate(year = climchangeyears)
    # Split file names
    popvars_new <- popvars_new |> mutate(
                       mate_cdmat_old = str_split_1(popvars$mate_cdmat[1], fixed("|")),
                       migrateout_cdmat_old = str_split_1(popvars$migrateout_cdmat[1], fixed("|")),
                       migrateback_cdmat_old = str_split_1(popvars$migrateback_cdmat[1], fixed("|")),
                       stray_cdmat_old = str_split_1(popvars$stray_cdmat[1], fixed("|")),
                       disperse_cdmat_old = str_split_1(popvars$disperseLocal_cdmat[1], fixed("|")))
  }
  if(!climate_change){
    popvars_new <- popvars |> mutate(mate_cdmat_old = mate_cdmat,
                                     migrateout_cdmat_old = migrateout_cdmat,
                                     migrateback_cdmat_old = migrateback_cdmat,
                                     stray_cdmat_old = stray_cdmat,
                                     disperse_cdmat_old = disperseLocal_cdmat)
  }
  
  # Make new file names
  cdmat_dir = paste0(output_directory,"cdmats")
  dir.create(file.path(cdmat_dir), showWarnings = FALSE) # Create output directory if it doesn't already exist
  
  popvars_new <- popvars_new |> mutate(mate_cdmat = new_file_name(cdmat_dir, mate_cdmat_old),
                                       migrateout_cdmat = new_file_name(cdmat_dir, migrateout_cdmat_old),
                                       migrateback_cdmat = new_file_name(cdmat_dir, migrateback_cdmat_old),
                                       stray_cdmat = new_file_name(cdmat_dir, stray_cdmat_old),
                                       disperse_cdmat = new_file_name(cdmat_dir, disperse_cdmat_old))
  # Copy matrices over
  for(i in 1:length(climchangeyears)){
    write_csv(read_csv(paste0(param_directory, popvars_new$mate_cdmat_old[i]), col_names = FALSE, show_col_types = FALSE),
              popvars_new$mate_cdmat[i],col_names = FALSE, quote = "none")
    write_csv(read_csv(paste0(param_directory, popvars_new$migrateout_cdmat_old[i]), col_names = FALSE, show_col_types = FALSE),
              popvars_new$migrateout_cdmat[i],col_names = FALSE, quote = "none")
    write_csv(read_csv(paste0(param_directory, popvars_new$migrateback_cdmat_old[i]), col_names = FALSE, show_col_types = FALSE),
              popvars_new$migrateback_cdmat[i],col_names = FALSE, quote = "none")
    write_csv(read_csv(paste0(param_directory, popvars_new$stray_cdmat_old[i]), col_names = FALSE, show_col_types = FALSE),
              popvars_new$stray_cdmat[i],col_names = FALSE, quote = "none")
    write_csv(read_csv(paste0(param_directory, popvars_new$disperse_cdmat_old[i]), col_names = FALSE, show_col_types = FALSE),
              popvars_new$disperse_cdmat[i],col_names = FALSE, quote = "none")
  }
  
  ## 3. Process remaining variables ##
  popvars_new  = mutate(popvars_new, mature_eqn_slope_f = str_split_1(popvars_new$mature_eqn_slope[1], "~")[1],
                    mature_eqn_slope_m = str_split_1(popvars$mature_eqn_slope[1], "~")[2],
                    mature_eqn_int_f = str_split_1(popvars$mature_eqn_int[1], "~")[1],
                    mature_eqn_int_m = str_split_1(popvars$mature_eqn_int[1], "~")[2],
                    mature_age = gsub("age", "", mature_default))
  ## 4. Remove unused variables and write to file
  if(popvars$mutationtype != 'random'){
    print(glue("Mutation type {popvars$mutationtype} not supported, using 'random'"))
  }
  if(runvars$sizecontrol == 'Y'){
    popvars_used <- c("xyfilename", "mate_cdmat", "matemoveno", "migrateout_cdmat",
                      "migrateback_cdmat", "stray_cdmat", "disperse_cdmat",
                      "mature_eqn_slope", "mature_eqn_int", "Egg_Mean_ans", "Egg_Mean_par1", "Egg_Mean_par2",
                      "Egg_Mortality", "offno", "loci", "growth_Loo", "growth_R0", "growth_temp_max",
                      "growth_temp_CV", "growth_temp_t0", "popmodel_par1", "mature_eqn_slope_f",
                      "mature_eqn_slope_m", "mature_eqn_int_f", "mature_eqn_int_m", "mature_age",
                      "popmodel", "startGenes", "muterate")
  }
  else{
    popvars_used <- c("xyfilename", "mate_cdmat", "matemoveno", "migrateout_cdmat",
                      "migrateback_cdmat", "stray_cdmat", "disperse_cdmat",
                      "Egg_Mortality", "offno", "loci",
                      "popmodel_par1",
                      "popmodel", "startGenes", "muterate")
  }
  # For gene initialization method "random", check that number of alleles is a single number
  if(genes_method == "random"){
    
    if(!is.numeric(popvars_new$alleles)){
      popvars_new$alleles = as.numeric(str_split_1(popvars_new$alleles, fixed(":")))[1]
      print(glue("Specifying different numbers of alleles for each locus is unsupported.",
            "Using {popvars_new$alleles} alleles per locus."))
    }
    popvars_used = c(popvars_used, "alleles")
  }
  if(climate_change){popvars_used = c("year", popvars_used)}
  popvars_out <- select(popvars_new, all_of(popvars_used)) |> mutate(startGenes = startGenes + 1)
  write_csv(popvars_out, popvars_file_out)
  print(glue('Processing run {run} of {nruns} finished'))
}




