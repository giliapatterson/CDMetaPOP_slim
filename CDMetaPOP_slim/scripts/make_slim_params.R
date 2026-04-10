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

#source("file_processing_functions.r")
source(paste0(script_directory,"/file_processing_functions.r"))
parser <- ArgumentParser()

parser$add_argument(
    "-d",
    "--parameter_directory",
    default = '../../example_runs/trout_genome/thermal_performance/',
    help = "Directory containing CDMetaPOP input files"
)

parser$add_argument(
    "-r",
    "--runvars_file_name",
    default = '../../example_runs/trout_genome/thermal_performance/RunVars.csv',
    help = "Name of RunVars file for CDMetaPOP"
)

parser$add_argument(
    "-o",
    "--output_directory",
    default = '../../example_runs/trout_genome/thermal_performance/slim_parameters_test/',
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
  # Slice run 1
  runvars <- runvars_all |> slice(run)
  
  # Directory to output SLiM parameter files
  output_directory = paste0(output_overall, "run", run, "/")
  # Create output directory if it doesn't already exist
  dir.create(file.path(output_directory), showWarnings = FALSE)
  
  #Add 1 to variables containing years to be consistent with SLiM
  runvars <- mutate(runvars,
                    cdclimgentime = add_one(as.character(cdclimgentime)))
  if(grepl("|", runvars$output_years, fixed = TRUE)){
    runvars <- mutate(runvars, output_years = add_one(as.character(output_years)))
  }
  if(has_name(runvars, 'mutation_output_years')){
    runvars <- mutate(runvars, mutation_output_years = add_one(as.character(mutation_output_years)),
      mutation_origin_year = mutation_origin_year + 1)
  }
  
  # Does this run change the climate?
  climchangeyears <- as.numeric(str_split_1(as.character(pull(runvars, cdclimgentime)), fixed("|"))) 
  climate_change <- length(climchangeyears) > 1
  
  # Check output format
  if(runvars$gridformat != 'genalex'){
    print(glue("gridformat {runvars$gridformat} is not supported. Using 'genalex' grid format."))
  }
  
  # SLiM code doesn't currently support all options of CDMetaPOP
  # Extract the options used by SLiM
  if(climate_change){runvars_used <- c("Popvars", "sizecontrol","runtime", "output_years", "cdclimgentime")}
  if(!climate_change){runvars_used <- c("Popvars", "sizecontrol", "runtime", "output_years")}
  if(has_name(runvars, 'mutation_output_years')){runvars_used <- c(runvars_used, 'mutation_output_years', 'mutation_output_subpops', 'mutation_origin_year')}
  runvars <- select(runvars, all_of(runvars_used))
  
  #### POPVARS ######
  # Change entry for Popvars to match the Popvars file used for SLiM
  popvars_file_out <- paste0(output_directory, "PopVars_slim.csv")
  runvars_out <- runvars |>
    mutate(Popvars = popvars_file_out)
  
  ### Process PopVars ###
  # Read in CDMetaPOP popvars
  popvars <- read_csv(paste0(param_directory, pull(runvars, Popvars)), show_col_types = FALSE) 
  
  ## 1. Process patchvars ##
  # Read in CDMetaPOP patchvars
  patchvars <- read_csv(paste0(param_directory, pull(popvars, xyfilename)), show_col_types = FALSE)
  # Output file for patchvars for slim
  patchvars_file_out <- paste0(output_directory, "PatchVars_slim.csv")
  # Convert patchID to be 0 indexed. This is necessary for SLiM
  patchvars <- mutate(patchvars, PatchID_old = PatchID, PatchID = row_number() - 1)
  if(has_name(runvars, 'mutation_output_subpops')){
    # Convert subpops to output mutations for in RunVars to SLiM PatchID format
    output_subpops_old = runvars_out$mutation_output_subpops |> str_split_1(fixed("|")) |> as.numeric()
    output_subpops_new = patchvars |> filter(PatchID_old %in% output_subpops_old) |> pull(PatchID)
    runvars_out$mutation_output_subpops = paste(output_subpops_new, collapse = "|")
  }
  
  ## 1. (a) Process Genes (if needed) ##
  if(!has_name(patchvars, "Genes Initialize")){
    genes_method = 'none'
  }
  if(has_name(patchvars, "Genes Initialize")){
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
  if(runvars$sizecontrol != 'Y'){
    if(grepl("~", classvars$Maturation)[1]){
      classvars <- mutate(classvars, Maturation_F = as.numeric(str_split_i(Maturation, "~", 1)),
                          Maturation_M = as.numeric(str_split_i(Maturation, "~", 2)))
    }
    else{
      classvars <- mutate(classvars, Maturation_F = Maturation,
                          Maturation_M = Maturation)
    }
    classvars_used = c("Age class", "Body Size Mean (mm)", "Body Size Std (mm)", "Distribution",
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
  
  ## Split patchvars by climate
  patchvars <- patchvars |> mutate(year = paste(climchangeyears, collapse = "|")) |>
                                     separate_longer_delim(everything(), "|")
  
  ## 1. (c) Remove variables that aren't used by SLiM and write to file ##
  patchvars_used <- c("year", "PatchID", "X", "Y", "K", "K StDev", "N0", "Mortality Eggs",
                      "Migration Out Prob", "Migration Back Prob", "Straying Prob",
                      "Dispersal Prob", "GrowthTemperatureOut", "GrowthTemperatureOutStDev",
                      "GrowDaysOut", "GrowDaysOutStDev", "GrowthTemperatureBack", "GrowthTemperatureBackStDev",
                      "GrowDaysBack", "GrowDaysBackStDev", "classvars", "Natal Grounds", "Migration Grounds")
  if(genes_method == "file"){
     patchvars_used = c(patchvars_used, "genes_file_slim")
  }
  if(has_name(patchvars, 'qtl_loci_initial') & !has_name(popvars, 'genome')){
    patchvars_used = c(patchvars_used, 'qtl_dpe_initial', 'qtl_loci_initial')
  }
  
  patchvars_out <- select(patchvars, all_of(patchvars_used))                    
  write_csv(patchvars_out, patchvars_file_out)
  
  # Update patchvars file name in popvars
  popvars <- popvars |> mutate(xyfilename = patchvars_file_out)
  
  ## 2. Process matrices ##
  # If the climate changes over time, popvars can have multiple rows, one for each
  # time step
  # Separate popvars by year
  popvars <- popvars |> separate_longer_delim(everything(), delim = "|")
  # If none of the variables in popvars change over time, replicate rows
  if(climate_change & nrow(popvars) == 1){popvars <- popvars |> slice(rep(1, length(climchangeyears)))}
  if(nrow(popvars) > 1 & nrow(popvars) != length(climchangeyears)){
    print(glue("Error: not enough values specified for CDClimGen"))
    print(glue("{length(climchangeyears)} values needed for at least one variable in PopVars"))
    print(glue("Only {nrow(popvars)} specified"))
  }
  if(climate_change){popvars <- mutate(popvars, year = climchangeyears)}
  
  ## Process matrices
  cdmat_dir = paste0(output_directory,"cdmats")
  dir.create(file.path(cdmat_dir), showWarnings = FALSE) # Create output directory if it doesn't already exist
  popvars_new <- mutate(popvars, mate_cdmat_old = mate_cdmat,
                        migrateout_cdmat_old = migrateout_cdmat,
                        migrateback_cdmat_old = migrateback_cdmat,
                        stray_cdmat_old = stray_cdmat,
                        disperse_cdmat_old = disperseLocal_cdmat)
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
  popvars_new = mutate(popvars_new, mature_age = gsub("age", "", mature_default))
  if(runvars$sizecontrol == 'Y'){
    popvars_new  = mutate(popvars_new, mature_eqn_slope_f = str_split_i(mature_eqn_slope, "~", 1),
                      mature_eqn_slope_m = str_split_i(mature_eqn_slope, "~", 2),
                      mature_eqn_int_f = str_split_i(mature_eqn_int, "~", 1),
                      mature_eqn_int_m = str_split_i(mature_eqn_int, "~", 2))
  }
  if(grepl("~", popvars_new$mature_age[1])){
    popvars_new <- mutate(popvars_new,
                          mature_age_f = str_split_i(mature_age, "~", 1),
                          mature_age_m = str_split_i(mature_age, "~", 1))
  }
  if(!grepl("~", popvars_new$mature_age[1])){
    popvars_new <- mutate(popvars_new,
                          mature_age_f = mature_age,
                          mature_age_m = mature_age)
  }
  ## 4. Remove unused variables and write to file
  # First, set default genetic parameters if loci is not specified
  if(!has_name(popvars_new, 'loci')){
    popvars_new$loci = 0
    popvars_new$startGenes = 0
    popvars_new$muterate = 0
  }
  if(has_name(popvars, 'mutationtype')){
    if(popvars$mutationtype[1] != 'random'){
      print(glue("Mutation type {popvars$mutationtype} not supported, using 'random'"))
    }
  }
  if(runvars$sizecontrol == 'Y'){
    popvars_used <- c("xyfilename", "mate_cdmat", "matemoveno", "migrateout_cdmat",
                      "migrateback_cdmat", "stray_cdmat", "disperse_cdmat",
                      "mature_eqn_slope", "mature_eqn_int", "Egg_Mean_ans", "Egg_Mean_par1", "Egg_Mean_par2",
                      "Egg_Mortality", "offno", "loci", "growth_Loo", "growth_R0", "growth_temp_max",
                      "growth_temp_CV", "growth_temp_t0", "popmodel_par1", "mature_eqn_slope_f",
                      "mature_eqn_slope_m", "mature_eqn_int_f", "mature_eqn_int_m", "mature_age_f", "mature_age_m",
                      "popmodel", "startGenes", "muterate")
  }
  if(runvars$sizecontrol != 'Y'){
    popvars_used <- c("xyfilename", "mate_cdmat", "matemoveno", "migrateout_cdmat",
                      "migrateback_cdmat", "stray_cdmat", "disperse_cdmat",
                      "Egg_Mortality", "offno", "loci", "mature_age_f", "mature_age_m",
                      "popmodel_par1",
                      "popmodel", "startGenes", "muterate")
  }
  # Is there a QTL?
  if(has_name(popvars, "genome_length")){
    popvars_used <- c(popvars_used, "genome_length", "qtl_prop_genome", "qtl_pheno_eff", "qtl_env_variable",
                      "qtl_muterate", "qtl_recrate", "qtl_ve", "qtl_fit_sd")
    if(has_name(popvars, "qtl_mutations_initial")){
      popvars_used <- c(popvars_used, "qtl_mutations_initial")
    }
  }
  # Genome file
  if(has_name(popvars, "genome")){
    genome <- read_csv(paste0(param_directory, popvars$genome[1]), show_col_types = FALSE)
    genome_outfile <- paste0(output_directory, "genome.csv")
    write_csv(genome, genome_outfile)
    popvars_new$genome <- genome_outfile
    # Is this a QTL model with one trait?
    if(has_name(popvars, "qtl_pheno_eff")){
      popvars_used <- c(popvars_used, "genome", "qtl_prop_genome", "qtl_pheno_eff", "qtl_env_variable",
                        "qtl_ve", "qtl_fit_sd")
      if(has_name(popvars, "qtl_mutations_initial")){
        popvars_used <- c(popvars_used, "qtl_mutations_initial")
      }
    }
    # Or is this a QTL model for thermal performance curve?
    if(has_name(popvars, "Topt_pheno_eff")){
      popvars_used <- c(popvars_used,
                        "genome", 
                        "qtl_prop_genome",
                        "Topt_baseline",
                        "epsilon_baseline",
                        "Topt_VE",
                        "epsilon_VE",
                        "Topt_cost",
                        "epsilon_cost",
                        "Topt_pheno_eff",
                        "epsilon_pheno_eff")
      if(has_name(popvars, "qtl_mutations_initial")){
        popvars_used <- c(popvars_used, "qtl_mutations_initial")
      }
    }
  }
  # For gene initialization method "random", check that number of alleles is a single number
  if(genes_method == "random"){
    if(grepl(":",popvars_new$alleles[1])){
      popvars_new$alleles = as.numeric(str_split_1(popvars_new$alleles[1], fixed(":")))[1]
      print(glue("Specifying different numbers of alleles for each locus is unsupported.",
            "Using {popvars_new$alleles} alleles per locus."))
    }
    popvars_used = c(popvars_used, "alleles")
  }
  if(climate_change){popvars_used = c("year", popvars_used)}
  popvars_out <- select(popvars_new, all_of(popvars_used)) |> mutate(startGenes = as.numeric(startGenes) + 1)
  write_csv(popvars_out, popvars_file_out)
  # Output runvars
  write_csv(runvars_out, paste0(output_directory, "RunVars_slim.csv"))
  print(glue('Processing run {run} of {nruns} finished'))
}




