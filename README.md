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
R CMD INSTALL VCBART_1.2.5.tar.gz
```

It is highly recommended that you install R version 4.0.0 or later before installing **VCBART**. Before installing **VCBART**, please ensure that you have set up an appropriate C++ toolchain for your system:

  + For macOS, we recommend using the [**macrtools**](https://github.com/coatless-mac/macrtools) to set up your system
  + For Windows, we recommend using Rtools, which can be downloaded [here](https://cran.r-project.org/bin/windows/Rtools/). Please make sure you download the version of Rtools corresponding to your R version (e.g., RTools45 for R version 4.5)
  + For Linux, we recommend following [these instructions](https://github.com/stan-dev/rstan/wiki/Configuring-C-Toolchain-for-Linux) from the Stan Development Team


## Basic usage and reproducibility

The directory `paper_examples` contains code to reproduce all figures from our paper.
It further contains code for a single replication of our synthetic simulation study; to reproduce the entire simulation study, one would need to run that code repeatedly.
Finally, it contains code used to pre-process and analyze the HRS dataset and the Philadelphia crime data.

