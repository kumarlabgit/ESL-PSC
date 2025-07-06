# Supplementary Code #

## Table of Contents ##

1. [Description and dependencies](#Description_and_dependencies)
2. [asr_convergence_checker.py](#asr_convergence_checker.py)
3. [busted_parallel_run.py](#busted_parallel_run.py)
4. [ccs_method.py](#ccs_method.py)
5. [csubst_sim_parallel_run.py](#csubst_sim_parallel_run.py)
6. [esl_multi_trait_parallel_run.py](#esl_multi_trait_parallel_run.py)
7. [esl_parallel_run.py](#esl_parallel_run.py)
8. [random_convergent_pairs.py](#random_convergent_pairs.py)
9. [run_busted_on_many_alignments.py](#run_busted_on_many_alignments.py)
10. [run_ccs_parallel.py](#run_ccs_parallel.py)
11. [run_conv_sim_grid.py](#run_conv_sim_grid.py)
12. [run_esl_psc_on_cluster.py](#run_esl_psc_on_cluster.py)
13. [run_many_csubst_simulations_one_input.py](#run_many_csubst_simulations_one_input.py)
14. [run_multi_trait_ccs_parallel.py](#run_multi_trait_ccs_parallel.py)
15. [run_multi_trait_esl_on_cluster.py](#run_multi_trait_esl_on_cluster.py)
16. [run_parallel_csubst_analyze.py](#run_parallel_csubst_analyze.py)
17. [run_random_species_sims.py](#run_random_species_sims.py)
18. [sim_preprocess_for_esl_benchmarks.py](#sim_preprocess_for_esl_benchmarks.py)
19. [summarize_conv_benchmark_data.py](#summarize_conv_benchmark_data.py)
20. [summarize_conv_multi_trait_benchmark_data.py](#summarize_conv_multi_trait_benchmark_data.py)
21. [TAAS_method.py](#TAAS_method.py)

## Description and dependencies ##

This directory includes scritps used to implement the benchmarking and simulation analyses in (Allard et al. 2024). Many of these are designed to be run on an HPC cluster using the Torque resource manager. Several of these scripts automatically generate and submit PBS job scripts in order to run many thousands of analyses in parallel. A `requirements.txt` file is include in this directory that lists dependencies for these scripts on an HPC cluster. In order to use these scripts, it is necessary to create and use a virtual environment (env) that conatins these dependencies. The path to the env, along with the locations of local scratch memory etc. will need to be adjusted in each script as needed for the cluster on which they are running.

#### Note on software version ####

These supplementary scripts were developed and used in conjunction with the version of ESL-PSC as of commit `be9533e6e656ddc47df2e4f18043b9219cb2ade0`, which corresponds to the code version used in (Allard et al. 2025). Subsequent changes to ESL-PSC (e.g., improvements to default output paths and command-line options) may break the hard-coded assumptions in the scripts below. We recommend checking out that commit before using any of the scripts:

```bash
git checkout be9533e6e656ddc47df2e4f18043b9219cb2ade0
```

Brief overviews of each included script are given below. Information on inputs and parameters can be viewed by running `python [script name] --help`

## asr_convergence_checker.py ##

This script identifies convergent amino acid residues across different species by comparing ancestral sequence alignments. The script takes as input directories containing ancestral sequence alignments, a file listing species in foreground clades for comparison, and generates a CSV file with convergence counts for each protein.

## busted_parallel_run.py ##

This is not designed for stand alone use. It manages an individual run of the BUSTED method on an alignment on an HPC cluster as part of a large number of runs orchestrated by the `run_busted_on_many_alignments.py` script. 

## ccs_method.py ##

This script implements the Convergence at Conservative Sites (CCS) method. It can be run as a stand alone script and is also called by the scripts that run CCS analyses on an HPC cluster.

## csubst_sim_parallel_run.py ##

This manages a single paralell run of a simulation using the CSUBST simulate function.  It is designed to be called as part of a large number of simulations run by the `run_conv_sim_grid.py` script on an HPC cluster.

## esl_multi_trait_parallel_run.py ##

This manages a single paralell run of ESL-PSC as part of a large number of runs on an HPC cluster orchestrated by the `run_multi_trait_esl_on_cluster.py` script on an HPC cluster.  This version is designed for implementation of runs across many combinations of randomly selected species designated as convergent in order to test ESL-PSC on simulated data across different numbers of species. "Multi trait" referrs to simulated traits, i.e. the randomly selected species used for simulations and testing.

## esl_parallel_run.py ##

This manages a single paralell run of ESL-PSC as part of a large number of runs on an HPC cluster orchestrated by the `run_esl_psc_on_cluster.py` script. This version is designed to run ESL-PSC on the same set of combinations of species from the same two clades (i.e. the null echolocation simulation analysis).

## random_convergent_pairs.py ##

This script generates random combinations of species in pairs with one convergent-designated and one control-designated species per pair.  Constraints are enforced in order to ensure that each combination meets the topological requiremnts to run ESL-PSC and CCS.  A file can be submitted to exclude a list of species from the selections in order to restrict the choices to a certian clade and to exclude outgroup speices to be used for CCS.

## run_busted_on_many_alignments.py ##

This script orchestrates runs of the BUSTED method on many alignments on an HPC cluster.  It requires that the Hyphy package be installed (https://github.com/veg/hyphy) as well as the hyphy-analyses repository (https://github.com/veg/hyphy-analyses). 

## run_ccs_parallel.py ##

This script orchestrates many runs of the CCS method on an HPC cluster. Outputs consist of csv files listing convergence counts for each simulated alignment in a given set and are stored in the parent directory of the directory containing the input alignments. 

## run_conv_sim_grid.py ##

This script orchestrates many runs of partition replacement simulations using the CSUBST simulate function.  It is designed to run simulations for many replicates involving a given number of randomly selected alignments per replicate, for which n simulated sites affected by simulated convergent selection will be substituted in for n sites in the original alignments.  For each replicate, a grid of the number of sites affected and foreground scaling factors can be defined such that simulation-containing alignments over all combinations of parameters will be generated.

## run_esl_psc_on_cluster.py ##

This orchestrates many parallel runs of ESL-PSC on an HPC cluster.  This version runs multi-combination ESL runs using the same species combinations for each replicate, as we did using the null echolocation combinations. This script takes a directory of subdirectories with simulation-containing alignments as generated by the `run_conv_sim_grid.py` script (but these must be translated to generate AA sequences first) and will set up ESL-PSC runs for each set of simulation-containing alignments by transferring compressed alignment files to local node scratch memory and substituting in the simulation-containing versions of alignments for their respective seed alignments, and then will run ESL-PSC on the full proteome-scale datasets and move the output gene ranks files to their appropriate directories in the original input directory. 

## run_many_csubst_simulations_one_input.py ##

This script can be run on a desktop computer.  It will generate a series of fully simulated alignments with a certain number each containing a defined percentage of sites affected by convergence, using the CSUBST simulate function.

## run_multi_trait_ccs_parallel.py ##

This script orchestrates a run on an HPC cluster of the CCS method on many sets of simulated alignments generated by the `run_random_species_sims.py` scripts. The output csv files are stored in the appropriate directories of each 

## run_multi_trait_esl_on_cluster.py ##

This script runs many whole proteome runs of ESL-PSC, like the `run_esl_psc_on_cluster.py` script, but in this case it is designed to run different species combinations (i.e. sets of randomly selected convergent-designated and control pairs) for which only a single combination is analyzed per ESL-PSC run. 

## run_parallel_csubst_analyze.py ##

This script orchestrates many parallel runs of the CSUBST analyze function on an HPC cluster. The output is one csv file generated by CSUBST for each individual alignment that is analyzed. These are moved to the same directory where the input alignment is located. CSUBST must be installed for the user on the HPC cluster for this to work.

## run_random_species_sims.py ##

This script is designed to run many simulations and generate partition-replaced alignments with n simulated sites substituted for n original sites in a set of seed alignments. This takes a directory containing text files specifying species combinations (as generated by the `random_convergent_pairs.py` script) and it generates a structure of nested subdirectories with one for each set of simulated alignments for each set of parameters for each replicate for each species combination.  One PBS job script is created for each replicate. 

## sim_preprocess_for_esl_benchmarks.py ##

This is not designed to be run in a stand alone way but is called by the `run_conv_sim_grid.py` script in order to generate a list of bash commands to run the CSUBST simulate function for each alignment that needs a simulated partition of sites.

## summarize_conv_benchmark_data.py ##

This script iterates through the directories containing ESL-PSC, CCS, and CSUBST analysis results from runs on many replicates simulation-containing alignments. For ESL-PSC, each analysis output file is the result of a full proteome scale run so the script simply counts the number of simulations in the top proteins. For CCS and CSUBST, results of full analyses by each method on the full background set of original simulation-free alignments must be included, so the scores from each set of simulated alignments can be compared with the empirical background to determine the rankings of the simulations according to the respective method. The output is a CSV file containing a summary of the number of simulation-containing alignments in the top ranking 100 proteins from each method. In addition, it also records the number of convergent sites on average detected by the CCS method for each simulated alignment. 

## summarize_conv_multi_trait_benchmark_data.py ##

This script perfroms the same summarization as the `summarize_conv_benchmark_data.py` script, but this is designed to summarize output from analyses of the random species combinations simulations for which only ESL-PSC were run due to computational constriants.

## TAAS_method.py ##

This script implements the Target species-specific Amino Acid Substitution method as described in (Zhang et al. 2014).
