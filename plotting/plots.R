
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
                              shape = 21, size=3.5, colour = "gray20", alpha = .65) +
  scale_fill_gradient(low = "gold", high = "darkred") + xlab("") + ylab("") +
  geom_line(aes(x = p1, y = p2), linetype = 1, size = 1.2, colour = "black", alpha=.9) +
  cowplot::theme_cowplot()
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
theta_chain <- h5read(h5file, "mh/theta_chain")
dft <- data.frame(theta_chain)
names(dft) <- c("Lambda", "Sigma")
dft <- dft[c(TRUE, rep(FALSE, 9)),]
dft$Iteration <- 10*(1:nrow(dft))
dftm <- melt(dft, id.vars = "Iteration", variable.name = "Parameter", value.name = "Value")

kern_plt <- ggplot(dftm) + geom_line(aes(x=Iteration, y=Value)) + 
  facet_wrap(~Parameter, scales = "free_y", nrow = 2) +
  cowplot::theme_cowplot() + ggtitle("Kernel parameter MCMC traces")




# Multiple chain vs ground truth ------------------------------------------

nchains <- h5read(h5file, "multi/nchains")
tmeans <- vector("list", nchains)
theta_means <- vector("list", nchains)

for(i in 1:nchains) {
  tm <- paste0("multi/chain", i)
  tchaini <- h5read(h5file, paste0(tm, "/tchain"))
  theta_chaini <- h5read(h5file, paste0(tm, "/theta_chain"))
  
  ## thin
  tchaini <- tchaini[c(TRUE, rep(FALSE, 9)), ]
  theta_chaini <- theta_chaini[c(TRUE, rep(FALSE, 9)),]
  
  tmeans[[i]] <- colMeans(tchaini)
  theta_means[[i]] <- colMeans(theta_chaini)
}

df_tmeans <- data.frame(do.call(cbind, tmeans))
df_thmeans <- data.frame(do.call(cbind, theta_means))

names(df_tmeans) <- as.character(1:nchains) # paste0("Chain", 1:nchains)
df_tmeans$Pseudotime <- t_gt

df_melted <- melt(df_tmeans, id.vars="Pseudotime", variable.name="Chain", value.name="Posterior")

pm_plt <- ggplot(df_melted) + 
  geom_point(aes(x=Pseudotime, y=Posterior, fill=Chain), size = 3, alpha = .8, shape = 21, colour = "gray20") +
  cowplot::theme_cowplot() + xlab("True pseudotime value") + ylab("Posterior mean")
pm_plt

# Save all plots ----------------------------------------------------------

ggsave("synth.png", synth_plt, width=8, height = 5)
ggsave("pst_traces.png", pst_trace_plt, width=8, height =5)
ggsave("kern_plt.png", kern_plt, width=8, height =5)
ggsave("pos_mean.png", pm_plt, width=8, height=5)

