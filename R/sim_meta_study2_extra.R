##################################################################################################
# (STUDY 2) SIMULATION OF META-ANALYSIS MODELS - Extra models to test for reviewer 1 comments 
# 
# Note: this is not included in final MS but it was used in the Reply letter to reviewer comments.
###################################################################################################

rm(list=ls())

# Load libraries
library(pacman) # checks if package is installed, if not installs it
p_load(metafor, MASS, clubSandwich, ape)


########## Load simulated datasets ---------------------------------------

### load 
load("data/simdat_extra.RDATA")

### parameters for current job
# name <- "study2_extra"
# scen <- 63
# k.studies <- 20
# sigma2.u <- 0.3  
# sigma2.s <- 0.3     
# rho <- 0.5   
# k.species <- 40
# sigma2.n <- 0.3
# sigma2.p <- 0.3
# mu <- 0.2


########### Set up sampling error matrix assuming rho=0.5 within study ------------------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.5)



###############
############### Notes for revision in MEE
############### 

## 1. Clustering is important - the performance is based on the cluster (we modified in Methods and Discussion etc). It can only account for part of it. 
## 2. We tried suggestion of Rev. 1. However this model does not provide useful biological information (which is what biologist want to obtain from meta-analysis on multiple species; the variance components for species/phylo/study).
## 3. In methodological literature it is uncommon (even if technically correct) + refer to James P. paper about higher level clusters to use with CRVE (but be careful when cluster has multiple single measures per cluster)
## 4. Alternative, analysis for the overall mean estimate could be run on a simpler model without the species and phylo levels and then obtain the overmean estimates uncertainty with robust CR2. Then fit phylogenetic model to get var estimates.  
## 5. CRVE, shouldn't be used as a separate method/replacement of other models. It helps improve precise/robust inference




########### Run model 1: multilevel model with sampling VCV (rho.hat=0.5) -------------------------------


# run model
pmod1 <- rma.mv(yi,
                V=V,
                random = list(~1 | stud.sp.id, ~ 1| id),
                test = "t",
                data = dat,
                dfs = k.studies - 1) #default "residuals" (k-1), other option "contain" (k.studies-1)

pmod1_cr0_study <- robust(pmod1, cluster = study, adjust=FALSE)
pmod1_cr1_study <- robust(pmod1, cluster = study, adjust=TRUE)

pmod1_cr0_species <- robust(pmod1, cluster = species, adjust=FALSE)
pmod1_cr1_species <- robust(pmod1, cluster = species, adjust=TRUE)

pmod1_cr0_stud.sp.id <- robust(pmod1, cluster = stud.sp.id, adjust=FALSE)
pmod1_cr1_stud.sp.id <- robust(pmod1, cluster = stud.sp.id, adjust=TRUE)



# model 1 results
pmod1_res <- data.frame(CR_methods = c("none", rep("CR0", 3), rep("CR1", 3)),
                        Cluster = c("none", "study", "species", "stud.sp.id",
                                    "study", "species", "stud.sp.id"),
                        mu_est = c(pmod1$b, pmod1_cr0_study$b, pmod1_cr0_species$b, pmod1_cr0_stud.sp.id$b,
                                   pmod1_cr1_study$b, pmod1_cr1_species$b, pmod1_cr1_stud.sp.id$b),
                        mu_se = c(pmod1$se, pmod1_cr0_study$se, pmod1_cr0_species$se, pmod1_cr0_stud.sp.id$se,
                                  pmod1_cr1_study$se, pmod1_cr1_species$se, pmod1_cr1_stud.sp.id$se),
                        mu_pval = c(pmod1$pval, pmod1_cr0_study$pval, pmod1_cr0_species$pval, pmod1_cr0_stud.sp.id$pval,
                                    pmod1_cr1_study$pval, pmod1_cr1_species$pval, pmod1_cr1_stud.sp.id$pval),
                        mu_cov = c(pmod1$ci.lb < mu & pmod1$ci.ub > mu,
                                   pmod1_cr0_study$ci.lb < mu & pmod1_cr0_study$ci.ub > mu,
                                   pmod1_cr0_species$ci.lb < mu & pmod1_cr0_species$ci.ub > mu,
                                   pmod1_cr0_stud.sp.id$ci.lb < mu & pmod1_cr0_stud.sp.id$ci.ub > mu,
                                   pmod1_cr1_study$ci.lb < mu & pmod1_cr1_study$ci.ub > mu,
                                   pmod1_cr1_species$ci.lb < mu & pmod1_cr1_species$ci.ub > mu,
                                   pmod1_cr1_stud.sp.id$ci.lb < mu & pmod1_cr1_stud.sp.id$ci.ub > mu)
)



########### Run model 2: phylogenetic multilevel model with sampling VCV (rho.hat=0.5) -------------------------------


pmod2 <- rma.mv(yi,
                V=V,
                random = list(~ 1 | species, ~ 1 | species.phylo, ~ 1 | study, ~ 1| id),
                R = list(species.phylo=P),
                test = "t",
                data = dat,
                dfs = k.studies - 1)

pmod2_cr0_study <- robust(pmod2, cluster = study, adjust=FALSE)
pmod2_cr1_study <- robust(pmod2, cluster = study, adjust=TRUE)

pmod2_cr0_species <- robust(pmod2, cluster = species, adjust=FALSE)
pmod2_cr1_species <- robust(pmod2, cluster = species, adjust=TRUE)

pmod2_cr0_stud.sp.id <- robust(pmod2, cluster = stud.sp.id, adjust=FALSE)
pmod2_cr1_stud.sp.id <- robust(pmod2, cluster = stud.sp.id, adjust=TRUE)


# model 2 full results
pmod2_res <- data.frame(CR_methods = c("none", rep("CR0", 3), rep("CR1", 3)),
                        Cluster = c("none", "study", "species", "stud.sp.id",
                                    "study", "species", "stud.sp.id"),
                        mu_est = c(pmod2$b, pmod2_cr0_study$b, pmod2_cr0_species$b, pmod2_cr0_stud.sp.id$b,
                                   pmod2_cr1_study$b, pmod2_cr1_species$b, pmod2_cr1_stud.sp.id$b),
                        mu_se = c(pmod2$se, pmod2_cr0_study$se, pmod2_cr0_species$se, pmod2_cr0_stud.sp.id$se,
                                  pmod2_cr1_study$se, pmod2_cr1_species$se, pmod2_cr1_stud.sp.id$se),
                        mu_pval = c(pmod2$pval, pmod2_cr0_study$pval, pmod2_cr0_species$pval, pmod2_cr0_stud.sp.id$pval,
                                    pmod2_cr1_study$pval, pmod2_cr1_species$pval, pmod2_cr1_stud.sp.id$pval),
                        mu_cov = c(pmod2$ci.lb < mu & pmod2$ci.ub > mu,
                                   pmod2_cr0_study$ci.lb < mu & pmod2_cr0_study$ci.ub > mu,
                                   pmod2_cr0_species$ci.lb < mu & pmod2_cr0_species$ci.ub > mu,
                                   pmod2_cr0_stud.sp.id$ci.lb < mu & pmod2_cr0_stud.sp.id$ci.ub > mu,
                                   pmod2_cr1_study$ci.lb < mu & pmod2_cr1_study$ci.ub > mu,
                                   pmod2_cr1_species$ci.lb < mu & pmod2_cr1_species$ci.ub > mu,
                                   pmod2_cr1_stud.sp.id$ci.lb < mu & pmod2_cr1_stud.sp.id$ci.ub > mu)
)



table(table(dat$stud.sp.id))
