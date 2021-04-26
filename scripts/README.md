## Reproducing results in paper

The scripts in this folder can be used to reproduce all of the examples and results from the main paper.

### Helper functions

While performing our simulations, we found it helpful to write wrapper functions for each of the methods considered.
These wrapper functions are impleneted in the scripts `*_wrapper.R`.

### Simulation with p = 5 and R = 20

To reproduce the main simulation, first run the script `generate_p5R20_data.R`.
Assuming you are working in the main directory (i.e. `..` relative to this folder), this script will save the synthetic dataset in the file `data/p5R20_data.RData`.

The main script to carry out the simulation study is `simulation_p5R20.R`.
It was written to be run as a batch job, where the fold number is provided as a command line argument.
For instance, to run the 10th fold and dump all of the printed output to a file, you could run the following from the command line
```
R --no-save --args 10 < scripts/simulation_p5R20.R > output.txt
```

If, instead, you wish to run a single fold interactive (e.g. within RStudio), you should comment out lines 17-19 of `simulation_p5R20.R` and assign the variable `sim_number` to the specific fold value you wish to run.


Once all folds are run, you can run the script `tabulate_p5R20.R` to combine and tabulate the simulation results.
Upon tabulating the results, you should be able to generate Figures 2 & 3

### Preparing the HRS dataset

We used publically available data from the Health and Retirement Study.
In order to replicate our main analysis, you will first need to follow the instructions [here](https://hrs.isr.umich.edu/data-products/access-to-public-data?_ga=2.12924862.827996587.1597077533-6044607.1597077533) to register and gain access to the public data.

### Predictive performance for HRS data

The predictive performance comparisons for HRS data can be run using the scripts `hrs_*.R`.





 
