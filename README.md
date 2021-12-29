# VCBART: Bayesian trees for varying coefficients

An R package for fitting a linear varying coefficient model using Bayesian Additive Regression Trees.
For more details about the VCBART procedure, see [our paper](https://arxiv.org/abs/2003.06416).

## Refactoring plan

1. No longer make copies of X,Y, and Z (just reference pointers to these objects directly)
2. In tree updates, avoid copying trees. 
3. Make a map with (bottom node pointer, vector of vector of observation indices at node)
4. Update Dirichlet concentration parameter with a MH step instead of discretizing
5. Better default tuning.
6. Handle the re-centering and re-scaling outside of the main functions for readability
 

