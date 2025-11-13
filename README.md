# VCBART: Bayesian trees for varying coefficients
<!-- badges: start -->
  [![R-CMD-check](https://github.com/skdeshpande91/VCBART/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/skdeshpande91/VCBART/actions/workflows/R-CMD-check.yaml)
  [![CRAN status](https://www.r-pkg.org/badges/version/VCBART)](https://CRAN.R-project.org/package=VCBART)
<!-- badges: end -->



An R package for fitting a linear varying coefficient model using Bayesian Additive Regression Trees.
For more details about the VCBART procedure, see [our paper](https://arxiv.org/abs/2003.06416).



## Installation

You can install the package using the following command in R:

```
devtools::install_github(repo = "skdeshpande91/VCBART").
```
Alternatively, you can fork or clone this repository and then build & install the package directly from the command line (be sure run this from the directory itself).
```
R CMD BUILD .
R CMD INSTALL flexBART_1.2.0.tar.gz
```

If you are using macOS, you will need to have previously set up the macOS toolchain for R.
This is typically something you have to do only once (and if you have used packages with **Rcpp** dependence, you probably already have done this).
But in case you haven't and are unable to install **VCBART**, you might find the following links helpful:

  + [R for macOS](https://cran.r-project.org/bin/macosx/tools/) for more information.
  + [Instructions from The Coatless Professor](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/). Note that these instructions may not be relevant for later versions of macOS (e.g. ones shipped with an Apple silicon processor)
  + [This StackOverflow post](https://stackoverflow.com/questions/69639782/installing-gfortran-on-macbook-with-apple-m1-chip-for-use-in-r/72997915#72997915), which outlines how to install the necessary C++ compiler, gfortran, and set the necessary paths. 


## Basic usage and reproducibility

The directory `paper_examples` contains code to reproduce all figures from our paper.
It further contains code for a single replication of our synthetic simulation study; to reproduce the entire simulation study, one would need to run that code repeatedly.
Finally, it contains code used to pre-process and analyze the HRS dataset and the Philadelphia crime data.

