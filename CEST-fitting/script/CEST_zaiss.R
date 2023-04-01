# Fit a 2 pool Z-spectrum

rm(list = ls())
gc()

modelName <- "CEST_zaiss"

setwd("~/Weizmann/rotation-1/CEST-fitting/script")
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- file.path(projectDir, "model")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")

require(rstan)
require(parallel)
require(bayesplot)
require(dplyr)
require(tidyr)
require(fdrtool)
require(ggplot2)
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Read data
Xdata <- read.csv(file.path(dataDir, "LP30_323.csv"),header = FALSE)
xZ = with(Xdata,V1)
Z = with(Xdata,V2)

## Create data set
data <- list(
  N = nrow(Xdata),
  R1a = 8,
  R2a = 393,
  w1 = 500,
  dw = -260,
  tp = 0.2,
  xZ = xZ,
  Z = Z
)

## Specify initial estimates
init = function() list(
  R2b = rnorm(1,27000,2000),
  f = rnorm(1,0.02,0.01),
  k = rnorm(1,285,10)
)

## Specify the variables for which you want history and density plots
parametersToPlot <- c("R2b", "f", "k")

## Additional variables to monitor
otherRVs <- c("Z_pred")

parameters <- c(parametersToPlot, otherRVs)

##################################################################
# Run Stan

nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
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
            #control = list(adapt_delta = 0.9),
            chains = nChains)

save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

##################################################################
## posterior distributions of parameters

dir.create(figDir)
dir.create(tabDir)

## open graphics device
# pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
# width = 6, height = 6, onefile = F)

options(bayesplot.base_size = 12,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "blue")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 5, myTheme = myTheme)
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, byChain = TRUE, 
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, 
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))

pairs(fit, pars = parametersToPlot[!grepl("rho", parametersToPlot)])

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
pred <- as.data.frame(fit, pars = "Z_pred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(Xdata %>% dplyr:::select(V1, V2))

p1 <- ggplot(pred, aes(x = xZ, y = Z))
p1 <- p1 + geom_point() +
  labs(x = "offsets (Hz)",
       y = "Z spectrum") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = xZ, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

# dev.off()