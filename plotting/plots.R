
# Generate plots for Bayesian Gaussian process latent variable models
# kieran.campbell@sjc.ox.ac.uk

library(ggplot2)
library(rhdf5)
library(dplyr)
library(reshape2)
library(magrittr)

# Functions for covariance matrix calcualtion ---------

cov_mat <- function(t, lambda, sigma) {
  tt <- as.matrix(dist(t))^2
  S <- exp(-lambda * tt) + sigma * diag(length(t))
  S
}

cross_cov_mat <- function(t1, t2, lambda, sigma) {
  n1 <- length(t1)
  n2 <- length(t2)
  tt <- matrix(0, nrow = n1, ncol = n2)
  for(i in 1:n1) {
    for(j in 1:n2) {
      tt[i, j] <- (t1[i] - t2[j])^2
    }
  }
  # tt <- tt + t(tt)
  S <- exp(-lambda * tt)
  S
}

# Analysis starts here --------------

setwd("/net/isi-scratch/kieran/GP/gpseudotime/plotting/")

h5file <- "../data/traces.h5"
h5ls(h5file)

# Synthetic parameters ----------------------------------------------------
t_gt <- h5read(h5file, "synthetic/t_gt")
X <- h5read(h5file, "synthetic/X")
sigma <- h5read(h5file, "synthetic/sigma")
n <- h5read(h5file, "synthetic/n")
lambda <- h5read(h5file, "synthetic/lambda")

# Read posterior traces from file ----------
tchain <- h5read(h5file, "mh/tchain")
lchain <- h5read(h5file, "mh/lambda_chain")
schain <- h5read(h5file, "mh/sigma_chain")
ll <- h5read(h5file, "mh/log_lik")
prior <- h5read(h5file, "mh/prior")

# Compute posterior mean estimates -----
burn <- as.integer(nrow(tchain) / 2)
tmap <- colMeans(tchain[burn:nrow(tchain),])
lmap <- colMeans(lchain[burn:nrow(lchain),])
smap <- colMeans(schain[burn:nrow(schain),])


tp <- runif(200, min = 0, max = 1)
mu <- matrix(NA, ncol = 2, nrow = length(tp))

for(i in 1:2) {
  K_map <- cov_mat(tmap, lmap[i], smap[i])
  K_star_transpose <- cross_cov_mat(tp, tmap, lmap[i], smap[i])
  matrix_prefactor <- K_star_transpose %*% solve(K_map)
  mu[,i] <- matrix_prefactor %*% X[,i]
}

# Main plot -------

df_x <- data.frame(cbind(X, tmap))
names(df_x) <- c("x1", "x2", "tmap")
df_mu <- data.frame(cbind(mu, tp))
names(df_mu) <- c("mu1", "mu2", "t")
df_mu <- arrange(df_mu, t)

embedded_plt <- ggplot() + geom_point(data = df_x, aes(x = x1, y = x2, fill = tmap), color = "gray20", 
                                      alpha = 0.65, size = 3.5, shape = 21) +
  geom_path(data = df_mu, aes(x = mu1, y = mu2, colour = t), size = 3, alpha = 0.7) +
  scale_fill_continuous_tableau() + scale_color_continuous_tableau()
embedded_plt + ggtitle("Posterior mean curve")

# Pseudotime Traces ------------------------------------------------------------------

## start by thinning the tchain a little
set.seed(123)
dftchain <- data.frame(tchain)
names(dftchain) <- paste0('t', 1:ncol(dftchain))

to_sample_bb <- paste0('t', sample(ncol(dftchain), 30, replace = FALSE))
to_sample_trace <- paste0('t', sample(ncol(dftchain), 10, replace = FALSE))

dftchain$Iteration <- 10*(1:nrow(dftchain))

dfm <- melt(dftchain, id.vars = "Iteration", variable.name = "Cell", value.name = "Pseudotime")

# data frame for converged trace (after burn in)
df_converge <- dfm %>% 
  filter(Cell %in% to_sample_trace) %>%
  filter(Iteration > (as.integer(max(Iteration) / 2)))


df_bigbang <- dfm %>%
  filter(Cell %in% to_sample_bb) %>%
  filter(Iteration < as.integer(max(Iteration) / 50))

pst_trace_plt <- ggplot(df_converge) + geom_line(aes(x=Iteration, y=Pseudotime, colour = Cell)) +
  ggtitle("MCMC traces") + guides(colour = FALSE) + cowplot::theme_cowplot()
pst_trace_plt

bb_plt <- ggplot(df_bigbang) + geom_line(aes(x=Iteration, y=Pseudotime, colour = Cell)) +
  ggtitle("MCMC traces") + guides(colour = FALSE) + cowplot::theme_cowplot()
bb_plt


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

