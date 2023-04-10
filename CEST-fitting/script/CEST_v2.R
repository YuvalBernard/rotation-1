
# Preliminaries -----------------------------------------------------------

rm(list = ls())
gc()
# Load libraries
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
library(ggplot2)
color_scheme_set("darkgray")
bayesplot_theme_set(new = theme_gray())
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
Xdata <- read.csv(file.path(dataDir, "LP30_323.csv"),header = FALSE)
xZ = with(Xdata,V1)
Z = with(Xdata,V2)

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
  R1b = runif(n=1, min=1, max=150),
  R2b = runif(n=1, min=15000, max=35000),
  f = runif(n=1, min=0.0001, max=0.1),
  k = runif(n=1, min=200, max=400),
  sigma = runif(n=1, min=0, max=0.1)
)


# Run MCMC ----------------------------------------------------------------

fit <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  init = init_estimates,
)

# Save fit object
temp_rds_file <- tempfile(fileext = ".RDS")
fit$save_object(file = temp_rds_file)

# Load fit object
fit <- readRDS(temp_rds_file)

# Posterior Summary Statistics --------------------------------------------

# Specify parameters for posterior diagnostics
pars_to_fit <- c("R1b", "R2b", "f", "k")

posterior <- as.array(fit$draws())

# Print summary statistics
fit$summary(pars_to_fit)

# Plot MCMC diagnostics 
mcmc_combo(
  posterior,
  combo = c("dens_overlay", "trace"),
  pars = pars_to_fit
)
mcmc_rhat(rhat(fit,pars = pars_to_fit)) + yaxis_text(hjust = 1)
mcmc_neff(neff_ratio(fit,pars = pars_to_fit)) + yaxis_text(hjust = 1)
mcmc_pairs(posterior, pars = pars_to_fit)


# Penalized Maximum Likelihood -------------------------------------------

# Find the mode (per parameter) of the joint posterior distribution 
fit_pml <- mod$optimize(
  data = data_list,
  init = init_estimates
)

# Summarize results
optimized_pars <- fit_pml$summary(pars_to_fit)

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

plot(xZ,Z)
lines(xZ,Z_rep_mean,col = "blue")