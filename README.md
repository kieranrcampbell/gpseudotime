# gpseudotime
Bayesian Gaussian Process Latent Variable Models (B-GP-LVM) for single-cell pseudotime ordering

## Structure of repo

* `bgplvm.jl` contains the MH algorithms written in Julia for inference
* `R_notebooks` has R markdown notebooks to convert the embeddings to HDF5
* `julia_notebooks` contains Jupyter notebooks for all the analysis (synthetic + moncole)
* `plotting` contains R scripts to create the plots for the paper
* `data` contains MCMC trace data in the form of HDF5 & CSV