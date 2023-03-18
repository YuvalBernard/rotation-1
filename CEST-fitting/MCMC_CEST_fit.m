function theta = MCMC_CEST_fit(data,n,varargin)

% This function utilizes Baye's rule to estimate the most probable set
% of parameters that fit the experimental Z spectrum.

% According to Baye's rule, for each parameter theta, its posterior
% pdf is proportional to its prior pdf times the likelihood of observing
% that parameter, given emphirical evidence.
% The prior should incorporate all knowledge of the parameter, and in the
% most extreme case of uncertainty should be taken as a uniform pdf.
% The likelihood, assuming additive white noise in experimental data,
% should be taken as a product of gaussians centered at the
% proposed simulated Z-spectrum.

% After constructing the priors and likelihood function, the
% Metropolis-Hastings algorithm is used to sample from the posteriors.

