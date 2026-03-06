#!/bin/bash
Rscript climate_change_McKenzie/make_param_files.R
python ../CDMetaPOP_slim/CDmetaPOP_slim.py -c2 -d climate_change_McKenzie -i RunVars.csv -o slim_output --no-filetime -s 300
Rscript climate_change_McKenzie/plot_qtl_model.R
