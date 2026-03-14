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

- To force rerunning simulations:

    ```
    python CDMetaPOP_slim.py -d ../example_runs/small_WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results_small --no-filetime -s 20329 --rerun
    ```

- To paralellize across replicates:

    ```
    python CDMetaPOP_slim.py -d ../example_runs/small_WCT_parameters/ -i RunVars_WCT.csv -o cdmetapop_slim_results_small --no-filetime -s 20329 --cores 20
    ```
    
    CDMetaPOP_slim will use as many cores as possible up to the number specified by `--cores`.

## Quantitative traits

### Model

CDMetaPOP SLiM allows users to model a quantitative trait linked to an environmental variable. In this model, an individuals phenotype is equal to
$$P_j = \sum_{h=1}^{2}\sum_{i=1}^{n}a_{ijh} + B + E_j$$
where $a_ijh$ is the phenotyic effect of the allele at locus $i$, on the $h^{th}$ copy of the genome of individual $j$, $B$ is the baseline phenotype, and $E_j$ is the additive effect of the environment for individual $j$. $E_j$ is drawn from a normal distribution with mean 0 and variance $V_E$. The baseline phenotype, $B$, is the mean value of the environmental variable specified in `qtl_env_variable` across all patches at the beginning of the simulation and does not change over time. The additive genetic variance of the trait is then
$$V_A = var\left(\sum_{h=1}^{2}\sum_{i=1}^{n}a_{ijh}\right)$$
and the heritability of the trait is $V_A/(V_E + V_A)$. This model does not include dominance or interactions between genotype and environment.

The simulation is initialized with no genetic variation and with a genome of length `genome_length`. At this point all individuals have the same phenotype, $B$. The user specifies the proportion of the genome where mutations that arise will affect the phenotype, `qtl_prop_genome`, as well as a mutation rate, `qtl_muterate`. When a new mutation arises in the quantitative trait portion of the genome, its phenotypic effect is drawn from the distribution of phenotypic effects specified by the user, `qtl_pheno_eff`. For example, if `qtl_pheno_eff = rnorm(1, 0, 5.0)`, then the phenotypic effect is drawn from a normal distribution with mean 0 and standard deviation 5.0. Mutations that arise in the remaining portion of the genome are neutral. The user also specifies a recombination rate for the genome, `qtl_recrate`. `qtl_muterate` and `qtl_recrate` can change over time, but `genome_length`, `qtl_prop_genome`, and `qtl_pheno_eff` cannot.

Selection on the trait takes place immediately after dispersal, before the "Packing back" step. Selection takes place through excess mortality. The probability an individual dies in the selection step is determined by the deviation of its phenotype from the environment it experiences. The user specifies which environmental variable the phenotype is linked to in `qtl_env_variable`. `qtl_env_variable` must be the name of a column in the PatchVars file. The probability an individual dies is:
$$f((P_j-\text{environment experienced}))/f(0)$$
where $f$ is the PDF of a normal distribution with mean 0 and standard deviation $\sigma = $ `qtl_fit_sd`,

$$f(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left( -\frac{1}{2}\left(\frac{x}{\sigma}\right)^{2}\right)$$

### Parameters

The parameters for the quantitative trait model are specified in PopVars and are summarized below:

`genome_length`: The total number of base pairs in the simulated genome. Cannot change over time.

`qtl_prop_genome`: The proportion of the genome where new mutations influence phenotype. In the remainder of the genome, mutations are neutral. Cannot change over time.

`qtl_pheno_eff`: A piece of code for drawing the phenotypic effect of a new mutation. This code should return a single value. Any of the distribution functions from R will work, e.g. "rnorm(1, 0, 5.0)" or "rexp(1, 1)". So will sampling from a fixed set of phenotypic effects, e.g. "sample(c(-1, 0, 1), 1)". Cannot change over time.

`qtl_env_variable`: A column name in PatchVars. It can be either an existing column used for other parts of the model, such as GrowthTemperatureBack, or a new column used only for the quantitative trait model. The column name specified cannot change over time, but the value of the environmental variable can. Use PatchVars and cdclimgen to do this. In addition, if GrowthTemperatureBack or GrowthTemperatureOut are used and the standard deviation of GrowthTemperatureBack or GrowthTemperatureOut is not 0, then the value of the environmental variable will vary randomly year to year. 

`qtl_muterate`: Mutation rate. Can change over time, e.g. 1e-8|0.0|0.0.

`qtl_recrate`: Recombination rate. Can change over time.

`qtl_ve`: V_E. Can change over time.

`qtl_fit_sd`: Standard deviation of the normal distribution used to determine fitness. Can change over time.

`qtl_mutations_initial`: Number of initial mutations.

In PatchVars:

Setting the initial genotypes of every individual in the patch (optional). All individuals in a patch will have the same initial genotype.

`qtl_loci_initial`: Number of loci underlying trait (e.g. 100)

`qtl_dpe_initial`: Distribution of effect sizes for those loci. (e.g. 'rep(1/200, 100)').

### Output

If the quantitative trait model is used, CDMetaPOP_slim will output two additional files, QTL_overall.csv and QTL subpops.csv. These files contain:

`phenotypic_variance`: Total variance in phenotype.

`Va`: Additive genetic variance.

`Ve`: Environmental variance.

`heritability`: Equal to $V_A/(V_A + V_E)$.

`qtl_alleles`: Number of alleles that influence the quantitative trait.

`avg_effect_size`: Average phenotypic effect.

`neutral_pi`: Genetic diversity for only the neutral portion of the genome.

`overall_pi`: genetic diversity for the entire genome.

`qtl_environment`: Value of the environmental variable in a patch.

`avg_phenotype`: Average phenotype across all individuals in a patch.

## Example runs

### Modeling thermal optimum as a quantitative trait

```
example_runs/climate_change_McKenzie
```

This example models the evolution of thermal tolerance as the climate warms. The example folder contains a bash script, `run_and_plot.sh` that generates the input files for CDMetaPOP_slim, runs the simulations, and plots population sizes, phenotypes, heritability, and genetic diversity over time. The resulting plots will be in the folder `example_runs/climate_change_McKenzie/qtl_plots`.

To run the script, navigate to the `example_runs` directory (`cd example_runs`), then run the bash script:

```
bash climate_change_McKenzie/run_and_plot.sh
```

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



