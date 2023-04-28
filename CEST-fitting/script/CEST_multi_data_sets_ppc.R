
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
library(ggthemes)
color_scheme_set("red")
bayesplot_theme_set(new = theme_gray())

modelName <- "CEST_multi_data_sets_ppc"

setwd("~/Weizmann/rotation-1/CEST-fitting/script")
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- file.path(projectDir, "model")
outDir <- file.path(modelDir, modelName)

if (dir.exists(outDir) == FALSE)
{
  dir.create(outDir)
  dir.create(figDir)
}

# Model Configuration -----------------------------------------------------

stan_file <- file.path(modelDir, paste(modelName, ".stan", sep = ""))
mod <- cmdstan_model(stan_file, compile = FALSE)

# Before compiling a stan model, check syntax for improvement with:
mod$check_syntax(pedantic = TRUE)

# Compile stan model
mod$compile(stanc_options = list("O1"))


# Data Configuration ------------------------------------------------------

# Read data
Xdata <- read.csv(file.path(dataDir, "LP30_dendrotes_CEST_exp_fit.csv"),header = FALSE, skip = 3)
xZ = with(Xdata,V1) * 16.546 * 9.4; # Convert from ppm to Hz and Lithium.
Z1 = with(Xdata,V20)
Z2 = with(Xdata,V22)
Z3 = with(Xdata,V24)
Z4 = with(Xdata,V26)

Z = c(Z1,Z2,Z3,Z4)

# Specify current experiment name. e.g. 'LP30_323_500_{sim/exp}'
# for Material: LP30, temperature: 323K, w1: 500Hz
expName <- "LP30_323_500_1000_1500_2000_exp"

# Create data set
data_list <- list(
  K = 4,
  N = c(length(Z1),length(Z2),length(Z3),length(Z4)),
  R1a = 8,
  R2a = 393,
  w1 = c(500,1000,1500,2000),
  dw = -260 * 9.4 * 16.546,
  tp = 0.2,
  xZ = rep(xZ,times = 4),
  Z = c(Z1,Z2,Z3,Z4)
)

# Run MCMC ----------------------------------------------------------------

fit <- mod$sample(
  data = data_list,
  chains = 1,
  parallel_chains = 1,
  iter_sampling = 3000,
  fixed_param = TRUE
)

# Save fit object
rds_file <- file.path(outDir, paste(expName, "_fit.RDS", sep = ""))
fit$save_object(file = rds_file)

# Load fit object
rds_file <- file.path(outDir, paste(expName, "_fit.RDS", sep = ""))
fit <- readRDS(rds_file)

# Posterior Summary Statistics --------------------------------------------

# Specify parameters for posterior diagnostics
pars_to_fit <- c("R1b", "R2b", "f", "k", "sigma")

posterior <- as.array(fit$draws())

# Posterior Predictive Checking -------------------------------------------

Z_rep <- as_draws_matrix(fit$draws("Z_sim"), .nchains = 1)
Z_rep_mean = colMeans(Z_rep)

ppc_dens_overlay(Z,Z_rep[1:50, ])
ggsave(file.path(figDir, paste(expName,"_ppc_dens_overlay.pdf", sep = "")),dpi = 300)

plot(rep(xZ,times=4),Z)
lines(rep(xZ,times=4),Z_rep_mean,col = "blue")

# Check if Z_tilde (generated from model) fits the data
Z_pred <- as_draws_matrix(fit$draws("Z_tilde"), .nchains = 1)
Z_pred_mean = colMeans(Z_pred)

ppc_dens_overlay(Z,Z_pred[1:25, ])
ggsave(file.path(figDir, paste(expName,"_ppc_dens_overlay_tilde.pdf", sep = "")),dpi = 300)

plot(xZ,Z)
lines(xZ,Z_pred_mean, col = "red")
