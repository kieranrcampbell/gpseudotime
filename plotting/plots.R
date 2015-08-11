
# Generate plots for Bayesian Gaussian process latent variable models

library(ggplot2)
library(rhdf5)
#library(cowplot)
library(dplyr)
library(reshape2)

setwd("/net/isi-scratch/kieran/GP/gpseudotime/plotting/")

h5file <- "../data/traces.h5"
h5ls(h5file)

# Synthetic parameters ----------------------------------------------------
t_gt <- h5read(h5file, "synthetic/t_gt")
X <- h5read(h5file, "synthetic/X")
sigma <- h5read(h5file, "synthetic/sigma")
n <- h5read(h5file, "synthetic/n")
lambda <- h5read(h5file, "synthetic/lambda")

pred <- h5read(h5file, "ppm/X")

df_synth <- data.frame(x = X[,1], y = X[,2], p1 = pred[,1], p2 = pred[,2], Pseudotime = t_gt)

synth_plt <- ggplot(df_synth) + geom_point(aes(x=x, y=y, fill = Pseudotime), 
                              shape = 21, size=3.5, colour = "gray70", alpha = .65) +
  scale_fill_gradient(low = "gold", high = "darkred") + xlab("") + ylab("") +
  geom_line(aes(x = p1, y = p2), linetype = 1, size = 1.5, colour = "gray10", alpha=.5)
synth_plt


# Pseudotime Traces ------------------------------------------------------------------

tchain = h5read(h5file, "mh/tchain")

## start by thinning the tchain a little
set.seed(123)
ttchain = tchain[c(TRUE, rep(FALSE, 9)),]
dfttchain <- data.frame(ttchain)
names(dfttchain) <- paste0('t', 1:ncol(dfttchain))

to_sample <- sample(ncol(dfttchain), 30, replace = FALSE)
dfttchain <- dfttchain[,to_sample]
dfttchain$Iteration <- 10*(1:nrow(dfttchain))

dfm <- melt(dfttchain, id.vars = "Iteration", variable.name = "Cell", value.name = "Pseudotime")

pst_trace_plt <- ggplot(dfm) + geom_line(aes(x=Iteration, y=Pseudotime, colour = Cell)) +
  ggtitle("MCMC traces") + guides(colour = FALSE) + cowplot::theme_cowplot()
pst_trace_plt


# Kernel parameter traces -------------------------------------------------
theta_chain = h5read(h5file, "mh/theta_chain")
dft <- data.frame(theta_chain)
names(dft) <- c("Lambda", "Sigma")
dft <- dft[c(TRUE, rep(FALSE, 9)),]
dft$Iteration <- 10*(1:nrow(dft))
dftm <- melt(dft, id.vars = "Iteration", variable.name = "Parameter", value.name = "Value")

ggplot(dftm) + geom_line(aes(x=Iteration, y=Value)) + 
  facet_wrap(~Parameter, scales = "free_y", nrow = 2) +
  cowplot::theme_cowplot() + ggtitle("Kernel parameter MCMC traces")




# Save all plots ----------------------------------------------------------



