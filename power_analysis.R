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

# set up expt. design data frame
expdat <- expand.grid(sbe.geno = as.integer(c(0:2)), indiv = factor(1:12), fruit = factor(1:5))
expdat$indiv <- factor(paste(expdat$sbe.geno, expdat$indiv, sep='.'))

#set.seed(101)
#rm(list=".Random.seed", envir=globalenv())
nsim <- 100
#beta <- c(1, -0.35)
beta <- c(1, -0.2)
#theta <- 0.05
theta <- 0.1
sigma <- 0.5

ss <- simulate(~ sbe.geno + (1 | indiv), nsim=nsim, family=gaussian, newdata=expdat,
               newparams=list(theta=theta, beta=beta, sigma=sigma))

expdat$resp <- ss[, 1] # Took simulation 1 and stuck it in 
fit1 <- lmer(resp ~ sbe.geno + (1 | indiv), data=expdat)

fitsim <- function(i) {
  rf <- refit(fit1, ss[[i]])
  tmp <- numeric(length=4)
  tmp[1:3] <- coef(summary(rf))["sbe.geno",]
  tmp[4] <- -2 * pt(abs(tmp[3]), get_ddf_Lb(rf, fixef(rf)), lower.tail = TRUE, log.p = TRUE)
  tmp
}

#t1 <- system.time(fitAll <- laply(seq(nsim), function(i) fitsim(i)))
#t1
fitAll <- laply(seq(nsim), function(i) fitsim(i))
fitAll <- setNames(as.data.frame(fitAll), c("est", "stderr", "tval", "pval"))
head(fitAll)

with(fitAll, mean(pval < 0.05)) # power calculation

## TO DO: wrap the above code in a loop with arrays as B. Bolker suggests and test a range of sample sizes,
# effect sizes, and variances / errors
# THEN, plot the results