##########################
## SCRIPT to obtain posteriors for Wang and Loreau(2014)'s estimates of alpha, beta, gamma metapopulation temporal variability (http://doi.wiley.com/10.1111/ele.12292)
## from a generalized linear mixed model describing variations in abundance at each patch
## here version for a 9-patch experiment ran by Stefano Masier and Frederik Mortier
## v0 June 2019
## initial script writer and maintained: Maxime Dahirel
#########################


### IMPORTANT!!! the script assumes the associated raw dataset is pre-loaded in R under the name data


### metadata on raw dataset content (to document later)

#
#
#
#
#


####### STEP 0: loading needed packages
### Note:  this script uses the Bayesian language Stan for model fitting
### The RStan package containing Stan has a non standard installation process
### check online: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

library(tidyverse)

library(rstan) ## Stan

library(coda) ## post-processing packages
library(bayesplot)
library(tidybayes)
library(matrixStats)

library(brms) ## the interface we are using

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) ## options to use all parallel cores during model fitting (1 chain per core)
## if this slows down too much, change
Nchains <- 4
Niter <- 2000 ## some default seetings, to be adjusted if not enough for model convergence


library(QGglmm) ### important package to obtain variance-covariance matrices on the observed scale (http://www.genetics.org/content/204/3/1281)

####### END OF STEP 0


####### STEP 1: DATA CLEANING, WRANGLING AND PREPARATION


data$REPLICATE_fullID <- data$TREAT ### this is just a renaming to make code clearer

data$PATCH2 <- paste("P", data$PATCH, sep = "") ##there is bug here that sometimes paste two Ps instead of one when executing the whole script
## this creates errors when trying to fit model, I need to check this
## if this happens, simply restart from the txt file import and then proceed line by line through this part
data$PATCH2 <- str_remove(data$PATCH2, "[.]")
### the two lines above make patches names (a) usable as column names (letter as 1st character)
### and consistent with brms standards on variable names (best to avoid dots and underscore,
### as some functions will remove them for output names and then matching input and output becomes harder

Npatches <- length(unique(data$PATCH2)) ### counting how many patches there are

tab <- spread(data = data, key = PATCH2, value = AFEMA)
## reshape the dataset so that each of the 9 patch of the metapopulation has its own column
## this allows us to estimate the patch by patch variance covariance matrix, which is the basis of Wang and Loreau's method


### TO DO: there are some NA; need to manage them
### for the moment brms just ignore their row, but this may mean throwing some info that can be salvaged using the subset() acat


####### END OF STEP 1


####### STEP 2: FIT MODEL, SIMPLEST CASE

## We fit a multivariate generalized linear model to the abundance data, with one submodel for each of the 9 patches
## these models include "random" effects for metapopulation (maybe some replicate have on average higher/lower pop sizes than other)
## and for observation week, to account for tempral variability
## this is this 2nd effect we are most interested in

## the |P| and |Q| are "covariance" indices. random effects sharing the same index across models should be considered to belong to the same variance covariance matrix
## which is then estimated


## NB: you can start with a Poisson model instead but it's easy to see that it's does not fit well
## pp_check(mod, rsp="P11") or pp_check(mod,"stat_2d",resp="P11") so a negative binomial is the best bet
## posterior predictive checks also show that there is no evidence of zero-inflation

## this is the simplest model to sdemonstrate the method: no effect of treatment either on mean patch population (ie. no fixed effect)
## and the variance covariance matrix are common to all treatments

## 9 patchs, 9 submodels
bf11 <- bf(P11 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf12 <- bf(P12 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf13 <- bf(P13 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf21 <- bf(P21 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf22 <- bf(P22 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf23 <- bf(P23 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf31 <- bf(P31 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf32 <- bf(P32 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)
bf33 <- bf(P33 ~ 1 + (1 | P | REPLICATE_fullID) + (1 | Q | WEEK), family = negbinomial)

## a weakly informative model (see McElreath course and book for rationale)
prior <- c(
  set_prior("normal(0,5)", class = "Intercept", resp = c("P11", "P12", "P13", "P21", "P22", "P23", "P31", "P32", "P33")),
  set_prior("exponential(1)", class = "sd", resp = c("P11", "P12", "P13", "P21", "P22", "P23", "P31", "P32", "P33")),
  set_prior("lkj(2)", class = "cor") # midly skeptical of very high correlations (positive or negative)
)

## the model
 mod <- brm(mvbf(bf11 + bf12 + bf13 +
                bf21 + bf22 + bf23 +
                bf31 + bf32 + bf33),
 data = tab, iter = Niter, chains=Nchains, prior = prior
 )
# takes max 15 min to fit on my laptop

summary(mod) ## check if model has converged and how it looks
pp_check(mod,resp="P11") ## can the model generate the observed data? (there are many ppcheck possibles, see the help)
##in multivariate model, you need to specify the response (so here which patch)
####### END OF STEP 2




####### STEP 3: EXTRACTION AND CALCULATION OF RELEVANT METRICS FROM MODEL

## functions to calculate wang and loreau metrics given a variance covariance matrix and a vector of patch means
## NB need to add some error messages to defend against common errors (non numeric and non matrix inputs, VCV matrix and means not the same dimension...)
## but they work

alpha_wang_loreau <- function(varcorr, means) {
  return((sum(sqrt(diag(varcorr))) / sum(means))^2)
}

phi_wang_loreau <- function(varcorr) {
  return(sum(varcorr) / (sum(sqrt(diag(varcorr))))^2)
}

gamma_wang_loreau <- function(varcorr, means) {
  sqrt_w_bar <- sum(sqrt(diag(varcorr))) / 9 ## because 9 patches ### to modify so the function detects it from input
  return((sqrt(sum(varcorr)) / sum(means))^2)
}

# no function for beta 1, it's simply 1/phi

beta2_wang_loreau <- function(varcorr, means) {
  sqrt_w_bar <- sum(sqrt(diag(varcorr))) / 9 ## because 9 patches

  return((sqrt_w_bar^2 - (sum(varcorr)) / (9^2)) / mean(means)^2)
}



## now we have these functions, let's extract what we need from the model

## first the predicted posterior mean patch and metapop-level abundances


patch_by_patch_Nmean <- fitted(mod, re_formula = NA, summary = FALSE)[, 1, ]
metapopNmean <- rowSums(patch_by_patch_Nmean)
patchNmean <- metapopNmean / 9

### there are three types of patch in a metapop differing in level of connectedness
### we can calculate mean patch size for each of them and compare them
Nmean_corners <- rowMeans(patch_by_patch_Nmean[, c(1, 3, 7, 9)])
Nmean_center <- patch_by_patch_Nmean[, 5]
Nmean_sides <- rowMeans(patch_by_patch_Nmean[, c(2, 4, 6, 8)])


## second the temporal variance covariance matrix **on the observed scale**
## this *on the observed scale* is VERY important, because a GLMM output gives it by definition on the latent scale
## How do we go from latent to observed???
## We use the QGglmm package and guidelines in this paper (http://www.genetics.org/content/204/3/1281)

## first we extract what we need on the latent scale

## posterior of the latent scale temperal VCV matrix
vcv_latent_time <- VarCorr(mod, summary = FALSE)$WEEK$cov

## posterior of the latent scale replicate VCV matrix
vcv_latent_metapop <- VarCorr(mod, summary = FALSE)$REPLICATE_fullID$cov

vcv_latent_tot <- vcv_latent_time + vcv_latent_metapop

fixef_latent <- fixef(mod, summary = FALSE) ### fixed effects means on the latent scale



vcv_obs_time <- vcv_latent_time ## just to create a dummy array; values will be filled with correct values below

## see de Villemereuil et al Genetics paper and QGglmm packages for details
## but broad description is
## for each posterior sample, create a vector psi, then fill it with value for each patch
## then use this psi to get the observed VCV from the latent VCV
## repeat for all posterior samples

for (i in 1:nsamples(mod)) {
  psiM <- NA
  for (j in 1:9) {
    psiM[j] <- QGpsi(fixef_latent[i, j], vcv_latent_tot[i, j, j], d.link.inv = function(x) {
      exp(x)
    })
  }
  psiM <- diag(psiM)
  vcv_obs_time[i, , ] <- psiM %*% vcv_latent_time[i, , ] %*% t(psiM)
}


### calculate our posterior variability metrics
alpha <- NA
phi <- NA
beta2 <- NA
gamma <- NA
for (i in 1:nsamples(mod)) {
  alpha[i] <- alpha_wang_loreau(vcv_obs_time[i, , ], patch_by_patch_Nmean[i, ])
  phi[i] <- phi_wang_loreau(vcv_obs_time[i, , ])
  gamma[i] <- gamma_wang_loreau(vcv_obs_time[i, , ], patch_by_patch_Nmean[i, ])
  beta2[i] <- beta2_wang_loreau(vcv_obs_time[i, , ], patch_by_patch_Nmean[i, ])
}

varmetrics <- data.frame(
  metapopNmean = metapopNmean, patchNmean = patchNmean,
  Nmean_corners = Nmean_corners, Nmean_center = Nmean_center, Nmean_sides = Nmean_sides,
  alpha = alpha, beta1 = 1 / phi, beta2 = beta2, gamma = gamma, TREAT = "all"
)

varmetrics <- as_tibble(varmetrics)

varmetrics
####### END OF STEP 3





####### STEP 4: FIT MODEL ASSUMING DIFFERENCES AMONG TREATMENTS


## What if we want to test whether our metrics differ between levels of NO_R
## we can!!

## first we need to create a new temporal variable
tab$WEEK2 <- paste(tab$WEEK, tab$NO_R)
## this is just because when you have separate VCV matrices per treatment,
## the brms argument gr(,by=) does not want treatment to
## have the same random effect names (here weeks)

## then we fit the model as in the previous case
## the only differences are

## (1)there a fixed effect of treatment now ~0+ NO_R
## IMPORTANT: the treatment effect MUST be a factor variable
## why ~0+NO_R and not ~NO_R
## because the 2nd case would gives use intercept(level 1) and delta(leveln-level1) as fixed effects
## while the 1st syntax gives us directly intercept(lvl1) and intercept(lvl n) so save lines of code

## (2) see how the random effect syntax is modified using gr(...,by=...) to say that we want separate VCV matrices per treatment level
bf11a <- bf(P11 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf12a <- bf(P12 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf13a <- bf(P13 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf21a <- bf(P21 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf22a <- bf(P22 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf23a <- bf(P23 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf31a <- bf(P31 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf32a <- bf(P32 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)
bf33a <- bf(P33 ~ 0 + NO_R + (1 | P | REPLICATE_fullID) + (1 | Q | gr(WEEK2, by = NO_R)), family = negbinomial)

prior2 <- c(
  # set_prior("normal(0,5)",class="Intercept",resp=c("P11","P12","P13","P21","P22","P23","P31","P32","P33")), ## no actual intercept, so silenced out
  set_prior("normal(0,5)", class = "b", resp = c("P11", "P12", "P13", "P21", "P22", "P23", "P31", "P32", "P33")),
  set_prior("exponential(1)", class = "sd", resp = c("P11", "P12", "P13", "P21", "P22", "P23", "P31", "P32", "P33")),
  set_prior("lkj(2)", class = "cor")
)

mod2 <- brm(mvbf(bf11a + bf12a + bf13a +
                bf21a + bf22a + bf23a +
                bf31a + bf32a + bf33a),
data = tab, iter = Niter, chains=Nchains, prior = prior2
)
####### END OF STEP 4



####### STEP 5: TO DO: EXTRACT METRICS, MULTI TREATMENT CASE

##it's slightly different for first model, because structure is more complex
## but not hard, I just need to take the time to do it cleanly

## you can try in the meantime :)

# VarCorr(mod2,summary=FALSE)$WEEK$cov[1,,][1:9,1:9] ##access to varcorr submatrix first trt level




# note on how to extract some distributional parameters from a negbinomial model (not needed for now, but stored here just in case)
# phi=summary(mod)$spec_pars[,1]
# p_nbino= exp(fixef_latent)/(exp(fixef_latent)+phi)  ###normally it's parametrized with the mean so exp(intercept)
# theta=p_nbino/(1-p_nbino)
# conversion from shape (phi) to theta based on brms distrib vignettes and negbin wikipedia page
## shape of ngbin = phi
# p = mu/(mu+phi)
# theta = scale = p/(1-p)
