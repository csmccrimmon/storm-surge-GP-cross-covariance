The folders simulation and storm_surge contain code to obtain the results in sections three and four of the manuscript. Each folder contains separate scripts for the SEP, NONSEP and PP-ARCokrig and PP models. Each emulator requires functions from the scripts located in the functions folder. The data for both the simulation and application are located in the data folder. The data for the simulation was generated using the script “data/code_data_generation_simulation.R.” 

To obtain the results for each model in Table 1 of the manuscript, use the following scripts:
- SEP env_results_SEP.R
- NONSEP env_results_nonsep.R
- PP-ARCokrig env_results_parallel_partial_ARCokrig.R
- PP env_results_parallel_partial.R.

To obtain the results for each model in Table 2 of the manuscript, use the following scripts:
- SEP storm_surge_results_SEP.R
- NONSEP storm_surge_results_nonsep_cokriging.R
- PP-ARCokrig storm_surge_results_parallel_partial_ARCokrig.R
- PP storm_surge_results_parallel_partial.R. 

The columns of Table 3 of the manuscript can be recreated using the following script:
- storm_surge_results_nonsep_cokriging.R

The run time for each emulator in the numerical example and application is provided in Tables 1 and 2 of the manuscript. We included the run time in the final summary table for each script, but the numbers may vary slightly from what is displayed in the manuscript. It only takes around a minute to recreate the results in Table 3 of the manuscript. For the storm surge application, computing the energy score for each emulator requires at least 8 GB of local memory. Because of the memory requirement, we provided a separate table with the non energy score performance metrics and included code to calculate the energy score at the end of the script. The memory requirements for all other scripts are well within the limits of a standard local workstation.