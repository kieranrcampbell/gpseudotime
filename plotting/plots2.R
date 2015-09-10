library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggthemes)



setwd("/net/isi-scratch/kieran/GP/gpseudotime")

plot_all <- function(fnames) { 
  set.seed(123)
  tchain <- read.csv(fnames[1], header =F)
  tchain <- data.frame(tchain)
  names(tchain) <- paste0("t", 1:ncol(tchain))
  
  
  # Pseudotime chain plots --------------------------------------------------
  
  tchain <- tchain[,sample(1:ncol(tchain), 8)]
  tchain$Iteration <- 1:nrow(tchain)
  
  tm <- melt(tchain, id.vars = "Iteration", variable.name = "Cell", value.name = "Value")
  
  pst_trace_plt <- ggplot(tm) + geom_line(aes(x=Iteration, y=Value, colour = Cell)) +
    ggtitle("Pseudotime") + guides(colour = FALSE) + cowplot::theme_cowplot()
  pst_trace_plt
  
  
  
  # Sigma chain plots -------------------------------------------------------
  schain <- read.csv(fnames[2], header =F)
  schain <- data.frame(schain)
  
  names(schain) <- paste0("sigma", 1:2)
  schain$Iteration <- 1:nrow(schain)
  sc <- melt(schain, id.vars = "Iteration", variable.name = "Sigma", value.name = "Value")
  
  sc_trace_plt <- ggplot(sc) + geom_line(aes(x=Iteration, y=Value, colour = Sigma)) +
    ggtitle("Sigma") + cowplot::theme_cowplot()
  sc_trace_plt + scale_y_log10()
  
  
  # Lambda chain plots ------------------------------------------------------
  
  lchain <- read.csv(fnames[3], header =F)
  lchain <- data.frame(lchain)
  
  names(lchain) <- paste0("lambda", 1:2)
  lchain$Iteration <- 1:nrow(lchain)
  lc <- melt(lchain, id.vars = "Iteration", variable.name = "Lambda", value.name = "Value")
  
  lc_trace_plt <- ggplot(lc) + geom_line(aes(x=Iteration, y=Value, colour = Lambda)) +
    ggtitle("Lambda") +  cowplot::theme_cowplot()
  lc_trace_plt
  
  e_plt <- plot_mean_curve(fnames)
  
  plot_grid(pst_trace_plt, sc_trace_plt, lc_trace_plt, e_plt)
}



plot_mean_curve <- function(fname) {
  ## predicted expression
  burn <- as.integer(nrow(tchain) / 2)
  
  tchain <- read.csv(fname[1], header = F)
  schain <- read.csv(fname[2], header = F)
  lchain <- read.csv(fname[3], header = F)
  
  tmap <- colMeans(tchain[burn:nrow(tchain),])
  lmap <- colMeans(lchain[burn:nrow(lchain),])
  smap <- colMeans(schain[burn:nrow(schain),])
  
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
  
  tp <- runif(200, min = 0, max = 1)
  X <- read.csv("data/X.csv", header=F)
  X <- X[X[,1] > 0,]
  X[,1] <- (X[,1] - mean(X[,1])) / sd(X[,1])
  X[,2] <- (X[,2] - mean(X[,2])) / sd(X[,2])
  
  mu <- matrix(NA, ncol = 2, nrow = length(tp))
  
  for(i in 1:2) {
    K_map <- cov_mat(tmap, lmap[i], smap[i])
    K_star_transpose <- cross_cov_mat(tp, tmap, lmap[i], smap[i])
    matrix_prefactor <- K_star_transpose %*% solve(K_map)
    mu[,i] <- matrix_prefactor %*% X[,i]
  }
  
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
}


gnames <- list()
gnames[[1]] <- c("data/gamma1.0_tchain.csv", "data/gamma1.0_sigma_chain.csv", "data/gamma1.0_lambda_chain.csv")
gnames[[2]] <- c("data/gamma10.0_tchain.csv", "data/gamma10.0_sigma_chain.csv", "data/gamma10.0_lambda_chain.csv")
gnames[[3]] <- c("data/gamma50.0_tchain.csv", "data/gamma50.0_sigma_chain.csv", "data/gamma50.0_lambda_chain.csv")

plots <- lapply(gnames, plot_all)

pdf(file = "plotting/vary_gamma.pdf")
plots[[1]]
plots[[2]]
plots[[3]]
dev.off()

rnames <- list()
rnames[[1]] <- c("data/r0.1_tchain.csv", "data/r0.1_sigma_chain.csv", "data/r0.1_lambda_chain.csv")
rnames[[2]] <- c("data/r0.01_tchain.csv", "data/r0.01_sigma_chain.csv", "data/r0.01_lambda_chain.csv")
rnames[[3]] <- c("data/r0.001_tchain.csv", "data/r0.001_sigma_chain.csv", "data/r0.001_lambda_chain.csv")

rplots <- lapply(rnames, plot_all)
pdf(file = "plotting/vary_r.pdf")
rplots[[1]]
rplots[[2]]
rplots[[3]]
dev.off()

tnames <- c("data/template_large_iter_tchain.csv", 
            "data/template_large_iter_sigma_chain.csv", 
            "data/template_large_iter_lambda_chain.csv")
pdf("plotting/longiter.pdf", width=9, height=10)
plot_all(tnames)
dev.off()
