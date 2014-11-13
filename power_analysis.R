# This power analysis w/ lme4 follows Ben Bolker's RPubs tutorial at:
# http://rpubs.com/bbolker/11703
# 
# The scenario I am modeling is the F3 family phenotyping for SBE
# this will be done in F3 families of SBE//INV parental genotype:
#       HET//M82 & HET//PEN
# 
# The mixed model is:
# sugar or BRIX ~ genotype + (1 | individual / fruit)
# I could also / instead try this model:
# sugar or BRIX ~ genotype + (genotype | individual / fruit)
#       ...but it's prob. not necessary to assume diff't among indiv. var.
#
# I'll get the expected phenotypic effect sizes for Brix, fructose, and glucose
# from my regressions of the Zamir et al. data
#
# I need to extract mean and variance or std. dev., but the PEN SBE allele
# (in an otherwise M82 bkgd) roughly reduces fructose and glucose in the fruit
# by ~0.75 units => **investigate what units those are!!**
#
# I need to specify:
# - theta: std. dev. of each random effect
# - beta: fixed effects params. => this is where I would input the effect size
#                                  derived from the Zamir data
# - sigma: residual standard error
#
# For the simulations I need to alter:
# - beta -- just the sbe.geno effect, i.e. effect of 1 allele of SBE => -0.4 to -0.1
# - theta -- among fruit std. dev. => 0.01 to 0.3 (~0.05 or less seems reasonable to expect from the Zamir et al.data)
#      I could check what the typical among fruit variance is for fruits of processing tomato lines.
# - sigma -- residual error => 0.1 - 1 (~0.5 seems reasonable to expect from the Zamir et al.data)
# - # indiv. per genotype => 5 to 12
# - # fruits per indiv. => 3 to 6
#
# SHOULD I ADD BLOCKS TO THE EXPT. DESIGN??
#
#
#
# LASTLY, repeat for BRIX

library(lme4)
library(plyr)
library(ggplot2)
library(coefplot2)
library(pbkrtest)

# vectors for array's dimnames
betas <- c(-0.4, -0.3, -0.2) # coefficient values for sbe.geno
thetas <- c(0.01, 0.05, 0.1) # random effect std. dev. among individuals
sigmas <- c(0.1, 0.5, 1)
no.indiv <- c(6, 9, 12)
no.fruits <- c(3, 4, 5)
nsim <- 100
vals <- c("est", "stderr", "tval", "pval")

# set up an empty multi-dimensional array with the right dimensions and names
pow <- array(dim = c(3, 3, 3, 3, 3, 100, 4), dimnames = list(betas = betas, thetas = thetas, sigmas = sigmas, 
                                                              no.indiv = no.indiv, no.fruits = no.fruits, nsim = 1:nsim,
                                                              vals = vals) )

# function to refit the LMM for each simulated replicate
fitsim <- function(j) {
  rf <- refit(fit1, ss[[j]])
  tmp <- numeric(length=4)
  tmp[1:3] <- coef(summary(rf))["sbe.geno",]
  tmp[4] <- -2 * pt(abs(tmp[3]), get_ddf_Lb(rf, fixef(rf)), lower.tail = TRUE, log.p = TRUE)
  tmp
}

# Nested loop to iterate through all the paramater values to test
for(b in 1:length(betas)) {
  beta <- c(1, betas[b])
  for(t in 1:length(thetas)) {
    theta <- thetas[t]
    for(s in 1:length(sigmas)) {
      sigma <- sigmas[s]
      for(i in 1:length(no.indiv)) {
        Nindiv <- no.indiv[i]
        for(f in 1:length(no.fruits)) {
          Nfruit <- no.fruits[f]
          # set up data frame with expt. parameters. SBE genotype is modeled as co-dominant with 0=homoz.M82 and 2=homoz.PEN
          expdat <- expand.grid(sbe.geno = as.integer(c(0:2)), indiv = factor(1:Nindiv), fruit = factor(1:Nfruit)) # expand.grid() sets up a data frame with all combinations
          expdat$indiv <- factor(paste(expdat$sbe.geno, expdat$indiv, sep='.'))
          # simulate response data using lme4::simulate according to the desired LMM model
          ss <- simulate(~ sbe.geno + (1 | indiv), nsim=nsim, family=gaussian, newdata=expdat,
                         newparams=list(theta=theta, beta=beta, sigma=sigma))
          expdat$resp <- ss[, 1] # Take 1st simulated response values and stick them in 
          fit1 <- lmer(resp ~ sbe.geno + (1 | indiv), data=expdat) # fit LMM w/ 1st simulated rep
          fitAll <- laply(seq(nsim), function(j) fitsim(j)) # refit with all nsim datasets
          pow[b,t,s,i,f,,] <- fitAll # input fitAll into the correct 200 x 4 array (which is a matrix) at the bottom of the hierarchy
        }
      }
    }
  }
}
