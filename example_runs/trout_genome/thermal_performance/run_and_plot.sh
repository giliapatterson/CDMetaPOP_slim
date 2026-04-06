#!/bin/bash
# Make parameter files
Rscript make_param_files.r
# Plot thermal performance curve examples
Rscript tpc_plots.r
# Run simulations
python ../../../CDMetaPOP_slim/CDmetaPOP_slim.py -c8 -d . -i RunVars.csv -o slim_output --no-filetime -s 300
# Plot results
Rscript plot_qtl_model.R
Rscript plot_mutation_origins.R
