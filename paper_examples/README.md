## Organiziation

Code to reproduce each figure of our manuscript is provided in the directory `figure_code`, which contains a separate sub-directory for each figure.
The directory `simulation_study1` contains code to run a single replication of our first simulation study, which compared VCBART to several competitor models in a setting with $p = 5$ covariates, $R = 20$ modifiers, and $N = 1000$ total observations.
Our second simulation study investigated how VCBART scaled to larger datasets; code to run all 25 replications for a given value of $N$ is provided in the directory `simuulation_study2.`
Code to pre-process and analyze the HRS and Philadelphia crime data are available in the sub-directories `hrs_analysis` and `philly_analysis.` 
Although we do not provide the raw HRS data, we provide instructions for downloading the necessary datasets from [the HRS website](https://hrs.isr.umich.edu/data-products). 
