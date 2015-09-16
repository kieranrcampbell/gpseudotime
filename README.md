# gpseudotime

This repo accompanies our preprint [Bayesian Gaussian Process Latent Variable Models for pseudotime inference in single-cell RNA-seq data](http://biorxiv.org/content/early/2015/09/15/026872). To reproduce the figures in the paper:
* The Jupyter notebooks `gpseudotime_all.ipynb` and `monocle_analysis_inversegamma.ipynb` in julia_notebooks will reproduce the synthetic and 'monocle' workflows respectively
* The R scripts `synthetic_plots.R` and `monocle_plots.R` in plotting will then reproduce the figures

Note to construct the monocle representation you need to run `R_notebooks/vignette.Rmd` to get the Laplacian Eigenmaps representation. This all relies heavily on HDF5 (through the rhdf5 and HDF5.jl libraries).

## Using bgplvm.jl

The main MH algorithm is in `bgplvm.jl`. Briefly, it is invoked via
```julia
B_GPLVM_MH(X, n_iter, burn, thin, 
    t, tvar, lambda, lvar, sigma, svar, 
    r = 1, return_burn = false, cell_swap_probability = 0,
    gamma = 1.0)
```
where
* `X` - cell-by-feature matrix
* `t`, `lambda`, `sigma` - initial values for the markov chain (note all are vectors)
* `tvar`, `lvar`, `svar` - variances for the proposal distributions of t, lambda and sigma respectively
* `r` - repulsion parameter for Corp prior
* `return_burn` - should the burn period of the traces be returned?
* `cell_swap_probability` - leave as 0
* `gamma` - rate parameter for exponential prior on lambda


## Structure of repo

* `bgplvm.jl` contains the MH algorithms written in Julia for inference
* `R_notebooks` has R markdown notebooks to convert the embeddings to HDF5
* `julia_notebooks` contains Jupyter notebooks for all the analysis (synthetic + moncole)
* `plotting` contains R scripts to create the plots for the paper
* `data` contains MCMC trace data in the form of HDF5 & CSV