# VCBART: Bayesian trees for varying coefficients

An R package for fitting a linear varying coefficient model using Bayesian Additive Regression Trees.
For more details about the VCBART procedure, see [our paper](https://arxiv.org/abs/2003.06416).



## Installation

The package source files are contained in the sub-directory `VCBART`. 
To install, you can download those files and then build and install the package from the command line.
Alternatively, you can install the package using the following command in R:

```
devtools::install_github(repo = "skdeshpande91/VCBART", subdir = "VCBART").
```

If you are using macOS, you will need to have previously set up the macOS toolchain for R.
This is typically something you have to do only once (and if you have used packages with **Rcpp** dependence, you probably already have done this).
But in case you haven't and are unable to install **VCBART**, you might find the following links helpful:

+ [R for macOS](https://cran.r-project.org/bin/macosx/tools/) for more information.
+ [Instructions from The Coatless Professor](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/). Note that these instructions may not be relevant for later versions of macOS (e.g. ones shipped with an Apple silicon processor)


## Basic usage and reproducibility

The directory `paper_examples` contains code to reproduce all figures from our paper.
It further contains code for a single replication of our synthetic simulation study; to reproduce the entire simulation study, one would need to run that code repeatedly.
Finally, it contains code used to pre-process and analyze the HRS dataset and the Philadelphia crime data.

