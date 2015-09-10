
# Generate plots for Bayesian Gaussian process latent variable models
# kieran.campbell@sjc.ox.ac.uk

library(ggplot2)
library(rhdf5)
library(dplyr)
library(reshape2)
library(magrittr)
library(modeest)
library(coda)

thin <- 100

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
get_map <- function(x) mlv(x,method = "HSM")$M
burn <- as.integer(nrow(tchain) / 2)
tmap <- apply(tchain[burn:nrow(tchain),], 2, get_map )
lmap <- apply(lchain[burn:nrow(lchain),], 2, get_map )
smap <- apply(schain[burn:nrow(schain),], 2, get_map )


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
  scale_fill_continuous_tableau() + scale_color_continuous_tableau() + guides(fill = FALSE, colour = FALSE)
embedded_plt 

# Pseudotime Traces ------------------------------------------------------------------

## start by thinning the tchain a little
set.seed(123)
dftchain <- data.frame(tchain)
names(dftchain) <- paste0('t', 1:ncol(dftchain))

to_sample_bb <- paste0('t', sample(ncol(dftchain), 30, replace = FALSE))
to_sample_trace <- paste0('t', sample(ncol(dftchain), 10, replace = FALSE))

dftchain$Iteration <- thin*(1:nrow(dftchain))

dfm <- melt(dftchain, id.vars = "Iteration", variable.name = "Cell", value.name = "Pseudotime")

# data frame for converged trace (after burn in)
df_converge <- dfm %>% 
  filter(Cell %in% to_sample_trace) %>%
  filter(Iteration > (as.integer(max(Iteration) / 2)))

# data for big-bang initialisation
df_bigbang <- dfm %>%
  filter(Cell %in% to_sample_bb) %>%
  filter(Iteration < as.integer(max(Iteration) / 50))

pst_trace_plt <- ggplot(df_converge) + geom_line(aes(x=Iteration, y=Pseudotime, colour = Cell)) +
   guides(fill = FALSE, colour = FALSE) + cowplot::theme_cowplot()
# pst_trace_plt

bb_plt <- ggplot(df_bigbang) + geom_line(aes(x=Iteration, y=Pseudotime, colour = Cell)) +
    guides(colour = FALSE) + cowplot::theme_cowplot()
#bb_plt

# Multiple chain vs ground truth ------------------------------------------

df <- data.frame(t = t_gt, tmap = tmap)
tburn <- tchain[burn:nrow(tchain),]
thpd <- apply(tburn, 2, function(x) HPDinterval(mcmc(x))[1,])
df$lower <- thpd[1,]
df$higher <- thpd[2,]

comp_plt <- ggplot(df) + geom_point(aes(x = t, y = tmap), size = 3, shape = 2) +
  geom_errorbar(aes(ymax = higher, ymin = lower, x = t), width = 0.01, alpha = 0.5, colour = "darkred") +
  stat_function(fun = function(x) x ) + ylab("Pseudotime MAP estimate") +
  xlab("'True' pseudotime") 

# Save all plots ----------------------------------------------------------

ggsave("embedding.png", embedded_plt, width=8, height = 5)
ggsave("pst_traces.png", pst_trace_plt, width=8, height =5)
ggsave("bb.png", bb_plt, width=8, height =5)
ggsave("pos_mean.png", comp_plt, width=8, height=5)

grid <- plot_grid(embedded_plt, comp_plt, pst_trace_plt, bb_plt, 
                  nrow = 2, labels = c("A", "B", "C", "D"))

cowplot::ggsave(filename = "all.png", plot = grid, width = 3, height = 2, scale = 3.5)
