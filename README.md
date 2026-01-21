# CDMetaPOP Version 3S
Single species CDMetaPOP 3 with source code switched to SLiM for increased efficiency when implementing genetic models. Multi-species models are in progress. Models implemented and input file formats are identical to CDMetaPOP 3. For more info on these models, see the CDMetaPOP 3 documentation https://github.com/ComputationalEcologyLab/CDMetaPOP.

## How to install

Clone or download the CDMetaPOP_slim repository

```
git clone git@github.com:giliapatterson/CDMetaPOP_slim.git
```

Navigate to the folder containing the source code

```
cd CDMetaPOP_slim/CDMetaPOP_slim
```

Install the necessary packages using conda or mamba and activate the environment

```
conda env create -f environment.yml
conda activate cdmetapop_slim
```

or 

```
mamba env create -f environment.yml
mamba activate cdmetapop_slim
```

## How to run

- For the same behavior as CDMetaPOP3:

    ```
    python CDMetaPOP_slim.py -d ../example_runs/small_WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results_small
    ```

- To prevent CDMetaPOP from appending a timestamp to the end of the output folder:

    ```
    python CDMetaPOP_slim.py -d ../example_runs/small_WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results_small --no-filetime
    ```

- To set a random seed for reproducibility:

    ```
    python CDMetaPOP_slim.py -d ../example_runs/small_WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results_small --no-filetime -s 20329
    ```

    This option is very helpful for running on the cluster, because CDMetaPOP won't rerun replicates with the same seed. If some of the replicates don't finish due to time limits or out of memory errors, you can simply rerun the command with the same seed and CDMetaPOP will only rerun replicates that didn't finish the first time.

- To paralellize across replicates:

    ```
    python CDMetaPOP_slim.py -d ../example_runs/small_WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results_small --no-filetime -s 20329 --cores 20
    ```
    
    CDMetaPOP_slim will use as many cores as possible up to the number specified by `--cores`.

## Example runs

### Coastal cutthroat trout in the McKenzie

```
example_runs/McKenzie_CCT
```

Random initialization of 100 loci.

```
python CDMetaPOP_slim.py -d ../example_runs/McKenzie_CCT/ -i RunVars.csv -o cdmetapop_slim_results --no-filetime -s 233
```

Runtime: ~24 seconds in CDMetaPOP Version 3S, ~12 minutes in CDMetaPOP Version 3

### Westslope cutthroat trout - small version

```
example_runs/small_WCT_parameters
```

5 patches, genes initialized from file, climate changes over time, 5 runs of each set of parameters.

```
python CDMetaPOP_slim.py -d ../example_runs/small_WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results_small --no-filetime -s 20329 --cores 1
```

Runtime: With no parallelization, ~10 minutes in CDMetaPOP Version 3S, ~124 minutes in CDMetaPOP Version 3.

### Westslope cutthroat trout - big version

```
example_runs/WCT_parameters
```

415 patches, genes initialized from file, climate changes over time, 5 runs of each set of parameters.

```
python CDMetaPOP_slim.py -d ../example_runs/WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results --no-filetime -s 20329 --cores 8
```

Runtime: Parallelized over 8 cores: ~216 minutes in CDMetaPOP Version 3S



