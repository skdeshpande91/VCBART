# VCBART: Bayesian trees for varying coefficients

An R package for fitting a linear varying coefficient model using Bayesian Additive Regression Trees.
For more details about the VCBART procedure, see [our paper](https://arxiv.org/abs/2003.06416).

This branch is dedicated to a major refactoring (started Dec 29, 2021) that *hopefully* speeds up VCBART substantially, mainly by avoiding repeated tree traversals in each iteration of the MCMC procedure

## Refactoring plan

1. No longer make copies of X,Y, and Z (just reference pointers to these objects directly)

*Update* we can use .begin() with Rcpp objects but they're stored in column-major order. In the code it's useful to be able to do something like *(x + j + i*p) to get the j-th variable for i-th observation. So we should pass the transposes of X,Z as inputs to the main routine (do the transpose in a wrapper function in R)

2. New sufficient stat object: it's a map whose key is the node id and whose value is a vector of vector of integers. In outer vector, elements correspond to subjects. These inner elements record which observations from each subject get assigned to which bottom node. In main tree updating loop we populate the sufficient stat object with `tree_traversal`. This is the **only** time we do a full tree traversal (previous implementations had at least 3+ additional traversals). In `grow_tree` we copy the sufficient statistic map, add elements correspond to the newly proposed leaves, and delete the element corresponding to node that we proposed splitting. In `prue_tree` we combine the two elements of sufficient statistic map corresponding to the two leafs we are pruning and then erase those elements in the map. 

3. We need to write a log-integrated-likelihood function that accepts a sufficient statistic map as input. There is hopefully a lot of cancellation / opportunities to avoid repeated computation

4. Update Dirichlet concentration parameter with an MH step instead of discretizing

5. Do all the re-centering & re-scaling outside of the main functions (improves readability)  


