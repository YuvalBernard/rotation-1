
# Preliminaries -----------------------------------------------------------

# load cmdstanr and build cmdstan with optimizations
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# library(cmdstanr)
# check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# cpp_options = list(STAN_CPP_OPTIMS = TRUE, "CXXFLAGS += -O3 -march=native -mtune=native")
# cmdstan_make_local(cpp_options = cpp_options, append = TRUE)
# install_cmdstan(cores = 4) or rebuild_cmdstan(cores = 4) if already installed without optimizations


rm(list = ls())
gc()
# Load libraries
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
color_scheme_set("red")
bayesplot_theme_set(new = theme_default())

modelName <- "CEST"

setwd("~/Weizmann/rotation-1/CEST-fitting/script")
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- file.path(projectDir, "model")
outDir <- file.path(modelDir, modelName)


# Model Configuration -----------------------------------------------------

stan_file <- file.path(modelDir, paste(modelName, ".stan", sep = ""))
mod <- cmdstan_model(stan_file, compile = FALSE)

# Before compiling a stan model, check syntax for improvement with:
mod$check_syntax(pedantic = TRUE)

# Compile stan model
mod$compile(stanc_options = list("O1"))


# Data Configuration ------------------------------------------------------

# Read data
Xdata <- read.csv(file.path(dataDir, "LP30_323_500.csv"),header = FALSE)
xZ = with(Xdata,V1)
Z = with(Xdata,V2)

# Specify current experiment name. e.g. 'LP30_323_500_{sim/exp}'
# for Material: LP30, temperature: 323K, w1: 500Hz
expName <- "LP30_323_500_exp"

# Create data set
data_list <- list(
  N = length(xZ),
  R1a = 8,
  R2a = 393,
  w1 = 500,
  dw = -260,
  tp = 0.2,
  xZ = xZ,
  Z = Z
)

# Specify initial estimates
init_estimates = function() list(
  R1b_unif = runif(n=1, min=0, max=1),
  R2b_std = runif(n=1, min=0, max=1),
  f_unif = runif(n=1, min=0, max=1),
  k_unif = runif(n=1, min=0, max=1),
  sigma_unif = runif(n=1, min=0, max=1)
)


# Run MCMC ----------------------------------------------------------------

fit <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  init = init_estimates,
)

# Save fit object
rds_file <- file.path(outDir, paste(modelName, "_fit.RDS", sep = ""))
fit$save_object(file = rds_file)

# Load fit object
rds_file <- file.path(outDir, paste(modelName, "_fit.RDS", sep = ""))
fit <- readRDS(rds_file)

# Posterior Summary Statistics --------------------------------------------

# Specify parameters for posterior diagnostics
pars_to_fit <- c("R1b", "R2b", "f", "k")

posterior <- as.array(fit$draws())

# Print summary statistics
fit$summary(pars_to_fit)

# Plot MCMC diagnostics 
p <- mcmc_combo(
  posterior,
  combo = c("dens_overlay", "trace"),
  pars = pars_to_fit)
print(p)
ggsave(file.path(figDir, paste(expName,"_mcmc_combo.pdf", sep = "")), plot = p, dpi = 300)
mcmc_rhat(rhat(fit,pars = pars_to_fit)) + yaxis_text(hjust = 1)
ggsave(file.path(figDir, paste(expName,"_mcmc_rhat.pdf", sep = "")),dpi = 300)
mcmc_neff(neff_ratio(fit,pars = pars_to_fit)) + yaxis_text(hjust = 1)
ggsave(file.path(figDir, paste(expName,"_mcmc_neff.pdf", sep = "")),dpi = 300)
mcmc_pairs(posterior, pars = pars_to_fit)
ggsave(file.path(figDir, paste(expName,"_mcmc_pairs.pdf", sep = "")),dpi = 300)


# Penalized Maximum Likelihood -------------------------------------------

# Find the mode (per parameter) of the joint posterior distribution 
fit_pml <- mod$optimize(
  data = data_list,
  init = init_estimates
)

# Summarize results
optimized_pars <- fit_pml$summary(pars_to_fit)
optimized_pars

p1 <- mcmc_hist(fit$draws("f")) + 
  vline_at(fit_pml$mle("f"), size = 1.5)

p2 <- mcmc_hist(fit$draws("k")) + 
  vline_at(fit_pml$mle("k"), size = 1.5)

p3 <- mcmc_hist(fit$draws("R1b")) + 
  vline_at(fit_pml$mle("R1b"), size = 1.5)

p4 <- mcmc_hist(fit$draws("R2b")) + 
  vline_at(fit_pml$mle("R2b"), size = 1.5)

# Plot draws from posterior and add MLE estimate (solid line) 
cowplot::plot_grid(p1,p2,p3,p4, nrow = 2, ncol = 2)


# Posterior Predictive Checking -------------------------------------------

Z_rep <- as_draws_matrix(fit$draws("Z_rep"), .nchains = 1)
Z_rep_mean = colMeans(Z_rep)

ppc_dens_overlay(Z,Z_rep[1:25, ])
ggsave(file.path(figDir, paste(expName,"_ppc_dens_overlay.pdf", sep = "")),dpi = 300)

plot(xZ,Z)
lines(xZ,Z_rep_mean,col = "blue")

# Check if Z_tilde (generated from model) fits the data
Z_pred <- as_draws_matrix(fit$draws("Z_tilde"), .nchains = 1)
Z_pred_mean = colMeans(Z_pred)

ppc_dens_overlay(Z,Z_pred[1:25, ])

plot(xZ,Z)
lines(xZ,Z_pred_mean, col = "red")

# Profiling ---------------------------------------------------------------

fit$profiles()
