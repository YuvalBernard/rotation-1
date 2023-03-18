# Fit 2 pool CEST model.

rm(list = ls())
gc()

library(rstan)
library(bayesplot)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

modelName <- "CEST_2_pools"

## Relative paths assuming the working directory is the script directory
## containing this script
setwd("~/Weizmann/rotation-1/CEST-fitting/script")
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- file.path(projectDir, "model")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")

## Read data
Xdata <- read.csv(file.path(dataDir, "LP30_323.csv"),header = FALSE)

## Create data set
data <- with(Xdata,
             list(
               N = nrow(Xdata),
               w1 = V1[1],
               tp = V2[1],
               xZ = V3,
               Z = V4
             ))
## Specify the variables for which you want history and density plots
parametersToPlot <- c("R1a", "R1b", "R2a", "R2b", "f", "k", "dw", "sigma")

## Additional variables to monitor
otherRVs <- c("ZPred")

parameters <- c(parametersToPlot, otherRVs)

##################################################################
# Run Stan

nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinnning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

dir.create(outDir)

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin,
            chains = nChains)

save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

##################################################################
## posterior distributions of parameters

dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

options(bayesplot.base_size = 12,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "brightblue")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

posterior <- as.array(fit)
mcmc_acf(posterior, pars = parametersToPlot) + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 4, myTheme = myTheme)

mcmc_dens_overlay(posterior, parametersToPlot)
mcmc_dens(posterior, parametersToPlot) + myTheme

pairs(fit, pars = parametersToPlot)

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive distributions

pred <- as.data.frame(fit, pars = "fxa24Pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(xdata %>% dplyr:::select(cavg, fxa.inh.avg))

p1 <- ggplot(pred, aes(x = cavg, y = fxa.inh.avg))
p1 <- p1 + geom_point() +
  labs(x = "average ME-2 plasma concentration (ng/mL)",
       y = "average factor Xa inhibition (%)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = cavg, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()