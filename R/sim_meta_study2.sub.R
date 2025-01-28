###############################################################################
# (STUDY 2.sub) SIMULATION OF META-REGRESSION MODELS
###############################################################################

rm(list=ls())

# Load libraries
library(pacman) # checks if package is installed, if not installs it
p_load(metafor, MASS, clubSandwich, ape)


########## Load parameters conditions ---------------------------------------

### load job array
tab <- read.csv("job_array_study2.sub.csv")
#tab <-  scen.tab2.sub

### job number from pbs script
job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#job <- 5 # for testing locally

### parameters for current job
name <- tab$name[tab$job_number == job] 
scen <- tab$scenario[tab$job_number == job] 
seed <- tab$sim[tab$job_number == job] # sim is equal to job_number
k.studies <- tab$k.studies[tab$job_number == job]  
sigma2.u <- tab$sigma2.u[tab$job_number == job]   
sigma2.s <- tab$sigma2.s[tab$job_number == job]       
rho <- tab$rho[tab$job_number == job]    
k.species <- tab$k.species[tab$job_number == job]
sigma2.n <- tab$sigma2.n[tab$job_number == job]
sigma2.p <- tab$sigma2.p[tab$job_number == job]
mu <- tab$mu[tab$job_number == job]

### results directory
dir <- tab$save_location[tab$job_number == job]




########### Simulate parameters  -------------------------------------------

set.seed(seed)

# simulate effect sizes for each study (same settings as Cinar et al., 2022)
k.per.study <- round(rbeta(n = k.studies, shape1 = 1.5, shape2 = 3) * 39) + 1

# set up id variables
study <- rep(seq_len(k.studies), times=k.per.study) # study id 
k <- length(study)                                  # total number of estimates
id <- seq_len(k)                                    # id of between study effect sizes 
esid <- unlist(lapply(k.per.study, seq_len))        # id of within study effect sizes

# simulate species indices
species <- round(rbeta(k, 2, 2) * (k.species - 1)) + 1
species[sample(k, k.species)] <- seq_len(k.species)
species.phylo <- species

# simulate tree
tree.power <- 0.3 # set tree power - what is best here?
tree <- rtree(k.species, tip.label=seq_len(k.species))
tree <- compute.brlen(tree, power = tree.power)
P <- vcv(tree, corr=TRUE)
P <- P[order(as.numeric(rownames(P))), order(as.numeric(rownames(P)))]

### when the power parameter is 0, vcv inserts Nan to the off-diagonals. We
### want the correlation matrix to be an identity matrix in these conditions.
### So, we will replace Nans with 0s.
#P[is.nan(P)] <- 0

### simulate fixed effects
# set up parameters
b0 <- 0.2          # intercept
b1.1 <- 0.05       # measurement type (study level) => slope of measurement type effect difference ( between reference = 0 and 1)
b1.2 <- -0.02      # measurement type (study level) => slope of measurement type effect difference (between reference = 0 and 2)
b2 <- -0.02        # weight between species (species level) => continuous (log)
b3 <- 0.05         # sex (observation level) => slope of sex effect difference (reference = 0)

measurement <- rep(0:2, times=round(k.studies/3)) # create measurement variable

sigma2.ws <- 0.1   # within species variance in covariate
sigma2.bs <- 0.3   # between species variance in covariate
species_effect <- rnorm(k.species, 0, sqrt(sigma2.bs))   # (assume mean centered)     
species_effect <- species_effect[species] 
within_species_effect <- rnorm(k, 0, sqrt(sigma2.ws))  # (assume mean centered)


# set up fixed predictors
x1 <- measurement[study] # study level predictor
x2 <- species_effect + within_species_effect      # species-level predictor
x3 <- sample(0:1, size = k, replace = T)    # observation-level predictor


### simulate random effects
u.u <- rnorm(k, 0, sqrt(sigma2.u))
u.s <- rnorm(k.studies, 0, sqrt(sigma2.s))[study]
u.n <- rnorm(k.species, 0, sqrt(sigma2.n))[species]
u.p <- mvrnorm(1, mu=rep(0,k.species), Sigma=sigma2.p*P)[species]

### simulate sample variance 
vi <- rbeta(k, 2, 20)    # same settings as Cinar et al., 2022

#### simulate sampling error
#ei <- rnorm(k, 0, sqrt(vi))

### set up sampling error VCV matrix for within study dependent effect sizes
VCV <- matrix(0, nrow = k, ncol = k) 

# fill in the lower off-diagonal of each study cluster with covariance between effect sizes
for (i in 2:k) {
  for (j in 1:i) {
    if (study[i] == study[j]) {
      VCV[i,j] <- rho * sqrt(vi[i] * vi[j]) 
    }
  }
}

# VCV is a k x k matrix with sampling variance in diagonal
VCV[upper.tri(VCV)] <- t(VCV)[upper.tri(VCV)] # fill in upper diagonal
diag(VCV) <- vi

### simulate sampling error for dependent effect sizes
mi <- mvrnorm(n = 1, mu = rep(0, length(vi)), Sigma = VCV)



########### Get estimates  ------------------------------------------------------------------

yi <- mu + 
  b1.1 * (x1 == 1) + b1.2 * (x1 == 2) + ## measurement type (study-level)
  b2*x2 +                               ## weight (species-level)
  b3*x3 +                               ## sex (observation-level)
  u.u + u.s + u.n + u.p + mi


# get simulated data
dat <- data.frame(name = name, scenario = scen, seed = seed, job_number = job,
                  study.id = study, id = id, esid = esid, species.id = species,
                  yi = yi, vi = vi,
                  b1.1 = b1.1, b1.2 = b1.2, b2 = b2, b3 = b3, 
                  sigma2.ws = sigma2.ws, sigma2.bs = sigma2.bs,
                  x1 = x1, x2 = x2, x3 = x3,
                  u.u = u.u, u.s = u.s, u.n = u.n, u.p = u.p,  mi = mi)
# save simulated data in R file
save(list = "dat", file = paste0("data/simdat_", job, ".RDATA"))


########### Run model 1: phylogenetic multilevel model  ----------------------------------------------

# run model
ptm <- proc.time()
pmod1 <- try(rma.mv(yi,
                   vi,
                   random = list(~ 1 | species,
                                 ~ 1 | species.phylo, 
                                 ~ 1 | study,
                                 ~ 1| id),
                   mods = ~ factor(x1) + x2 + factor(x3),
                   R = list(species.phylo=P),
                   test = "t",
                   dfs = k.studies - 1), #default "residuals" (k-1), other option "contain" (k.studies-1)
            silent = TRUE)
pmod1_time <- (proc.time() - ptm)[3]

pmod1_cr0 <- robust(pmod1, cluster = study, adjust=FALSE)
pmod1_cr1 <- robust(pmod1, cluster = study, adjust=TRUE)




########### Run model 2: phylogenetic multilevel model with sampling VCV (rho.hat=0.2) -------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.2)

ptm <- proc.time()
pmod2 <- try(rma.mv(yi,
                   V=V,
                   random = list(~ 1 | species,
                                 ~ 1 | species.phylo,
                                 ~ 1 | study,
                                 ~ 1| id),
                   mods = ~ factor(x1) + x2 + factor(x3),
                   R = list(species.phylo=P),
                   test = "t",
                   dfs = k.studies - 1), 
            silent = TRUE)
pmod2_time <- (proc.time() - ptm)[3]

pmod2_cr0 <- robust(pmod2, cluster = study, adjust=FALSE)
pmod2_cr1 <- robust(pmod2, cluster = study, adjust=TRUE)




########### Run model 3: phylogenetic multilevel model with sampling VCV (rho.hat=0.5) -------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.5)

ptm <- proc.time()
pmod3 <- try(rma.mv(yi,
                    V=V,
                    random = list(~ 1 | species,
                                  ~ 1 | species.phylo,
                                  ~ 1 | study,
                                  ~ 1| id),
                    mods = ~ factor(x1) + x2 + factor(x3),
                    R = list(species.phylo=P),
                    test = "t",
                    dfs = k.studies - 1), 
             silent = TRUE)
pmod3_time <- (proc.time() - ptm)[3]

pmod3_cr0 <- robust(pmod3, cluster = study, adjust=FALSE)
pmod3_cr1 <- robust(pmod3, cluster = study, adjust=TRUE)





########### Run model 4: phylogenetic multilevel model with sampling VCV (rho.hat=0.8) -------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.8)

ptm <- proc.time()
pmod4 <- try(rma.mv(yi,
                    V=V,
                    random = list(~ 1 | species,
                                  ~ 1 | species.phylo,
                                  ~ 1 | study,
                                  ~ 1| id),
                    mods = ~ factor(x1) + x2 + factor(x3),
                    R = list(species.phylo=P),
                    test = "t",
                    dfs = k.studies - 1), 
             silent = TRUE)
pmod4_time <- (proc.time() - ptm)[3]

pmod4_cr0 <- robust(pmod4, cluster = study, adjust=FALSE)
pmod4_cr1 <- robust(pmod4, cluster = study, adjust=TRUE)






########### Extract model estimates  -------------------------------------------------------------------

# mu estimate 
pmod1.est <- pmod1$b[1]
pmod1.cr0 <- pmod1_cr0$b[1]
pmod1.cr1 <- pmod1_cr1$b[1]
pmod2.est <- pmod2$b[1]
pmod2.cr0 <- pmod2_cr0$b[1]
pmod2.cr1 <- pmod2_cr1$b[1]
pmod3.est <- pmod3$b[1]
pmod3.cr0 <- pmod3_cr0$b[1]
pmod3.cr1 <- pmod3_cr1$b[1]
pmod4.est <- pmod4$b[1]
pmod4.cr0 <- pmod4_cr0$b[1]
pmod4.cr1 <- pmod4_cr1$b[1]


# mu estimate MSE
pmod1.est.mse <- (mu - pmod1.est)^2  
pmod1.cr0.mse <- (mu - pmod1.cr0)^2
pmod1.cr1.mse <- (mu - pmod1.cr1)^2
pmod2.est.mse <- (mu - pmod2.est)^2
pmod2.cr0.mse <- (mu - pmod2.cr0)^2
pmod2.cr1.mse <- (mu - pmod2.cr1)^2
pmod3.est.mse <- (mu - pmod3.est)^2
pmod3.cr0.mse <- (mu - pmod3.cr0)^2
pmod3.cr1.mse <- (mu - pmod3.cr1)^2
pmod4.est.mse <- (mu - pmod4.est)^2
pmod4.cr0.mse <- (mu - pmod4.cr0)^2
pmod4.cr1.mse <- (mu - pmod4.cr1)^2


# mu confidence interval
pmod1.est.ci.ub <- pmod1$ci.ub[1]
pmod1.est.ci.lb <- pmod1$ci.lb[1]
pmod2.est.ci.ub <- pmod2$ci.ub[1]
pmod2.est.ci.lb <- pmod2$ci.lb[1]
pmod1.cr0.ci.ub <- pmod1_cr0$ci.ub[1]
pmod1.cr0.ci.lb <- pmod1_cr0$ci.lb[1]
pmod1.cr1.ci.ub <- pmod1_cr1$ci.ub[1]
pmod1.cr1.ci.lb <- pmod1_cr1$ci.lb[1]
pmod2.cr0.ci.ub <- pmod2_cr0$ci.ub[1]
pmod2.cr0.ci.lb <- pmod2_cr0$ci.lb[1]
pmod2.cr1.ci.ub <- pmod2_cr1$ci.ub[1]
pmod2.cr1.ci.lb <- pmod2_cr1$ci.lb[1]
pmod3.est.ci.ub <- pmod3$ci.ub[1]
pmod3.est.ci.lb <- pmod3$ci.lb[1]
pmod3.cr0.ci.ub <- pmod3_cr0$ci.ub[1]
pmod3.cr0.ci.lb <- pmod3_cr0$ci.lb[1]
pmod3.cr1.ci.ub <- pmod3_cr1$ci.ub[1]
pmod3.cr1.ci.lb <- pmod3_cr1$ci.lb[1]
pmod4.est.ci.ub <- pmod4$ci.ub[1]
pmod4.est.ci.lb <- pmod4$ci.lb[1]
pmod4.cr0.ci.ub <- pmod4_cr0$ci.ub[1]
pmod4.cr0.ci.lb <- pmod4_cr0$ci.lb[1]
pmod4.cr1.ci.ub <- pmod4_cr1$ci.ub[1]
pmod4.cr1.ci.lb <- pmod4_cr1$ci.lb[1]



# mu coverage
pmod1.cov <- pmod1.est.ci.lb < mu && pmod1.est.ci.ub > mu
pmod1.cr0.cov <- pmod1.cr0.ci.lb < mu && pmod1.cr0.ci.ub > mu
pmod1.cr1.cov <- pmod1.cr1.ci.lb < mu && pmod1.cr1.ci.ub > mu
pmod2.cov <- pmod2.est.ci.lb < mu && pmod2.est.ci.ub > mu
pmod2.cr0.cov <- pmod2.cr0.ci.lb < mu && pmod2.cr0.ci.ub > mu
pmod2.cr1.cov <- pmod2.cr1.ci.lb < mu && pmod2.cr1.ci.ub > mu
pmod3.cov <- pmod3.est.ci.lb < mu && pmod3.est.ci.ub > mu
pmod3.cr0.cov <- pmod3.cr0.ci.lb < mu && pmod3.cr0.ci.ub > mu
pmod3.cr1.cov <- pmod3.cr1.ci.lb < mu && pmod3.cr1.ci.ub > mu
pmod4.cov <- pmod4.est.ci.lb < mu && pmod4.est.ci.ub > mu
pmod4.cr0.cov <- pmod4.cr0.ci.lb < mu && pmod4.cr0.ci.ub > mu
pmod4.cr1.cov <- pmod4.cr1.ci.lb < mu && pmod4.cr1.ci.ub > mu



# sigma.n estimate (species level variance)
# in model: sigma^2.1
pmod1.sigma2.n <- pmod1$sigma2[1]
pmod1.cr0.sigma2.n <- pmod1_cr0$sigma2[1]
pmod1.cr1.sigma2.n <- pmod1_cr1$sigma2[1]
pmod2.sigma2.n <- pmod2$sigma2[1]
pmod2.cr0.sigma2.n <- pmod2_cr0$sigma2[1]
pmod2.cr1.sigma2.n <- pmod2_cr1$sigma2[1]
pmod3.sigma2.n <- pmod3$sigma2[1]
pmod3.cr0.sigma2.n <- pmod3_cr0$sigma2[1]
pmod3.cr1.sigma2.n <- pmod3_cr1$sigma2[1]
pmod4.sigma2.n <- pmod4$sigma2[1]
pmod4.cr0.sigma2.n <- pmod4_cr0$sigma2[1]
pmod4.cr1.sigma2.n <- pmod4_cr1$sigma2[1]



# sigma.n estimate MSE
pmod1.sigma2.n.mse <- (sigma2.n - pmod1.sigma2.n)^2
pmod1.cr0.sigma2.n.mse <- (sigma2.n - pmod1.cr0.sigma2.n)^2
pmod1.cr1.sigma2.n.mse <- (sigma2.n - pmod1.cr1.sigma2.n)^2
pmod2.sigma2.n.mse <- (sigma2.n - pmod2.sigma2.n)^2
pmod2.cr0.sigma2.n.mse <- (sigma2.n - pmod2.cr0.sigma2.n)^2
pmod2.cr1.sigma2.n.mse <- (sigma2.n - pmod2.cr1.sigma2.n)^2
pmod3.sigma2.n.mse <- (sigma2.n - pmod3.sigma2.n)^2
pmod3.cr0.sigma2.n.mse <- (sigma2.n - pmod3.cr0.sigma2.n)^2
pmod3.cr1.sigma2.n.mse <- (sigma2.n - pmod3.cr1.sigma2.n)^2
pmod4.sigma2.n.mse <- (sigma2.n - pmod4.sigma2.n)^2
pmod4.cr0.sigma2.n.mse <- (sigma2.n - pmod4.cr0.sigma2.n)^2
pmod4.cr1.sigma2.n.mse <- (sigma2.n - pmod4.cr1.sigma2.n)^2


# sigma.p estimate (phylogenetic level variance)
# in model: sigma^2.2
pmod1.sigma2.p <- pmod1$sigma2[2]
pmod1.cr0.sigma2.p <- pmod1_cr0$sigma2[2]
pmod1.cr1.sigma2.p <- pmod1_cr1$sigma2[2]
pmod2.sigma2.p <- pmod2$sigma2[2]
pmod2.cr0.sigma2.p <- pmod2_cr0$sigma2[2]
pmod2.cr1.sigma2.p <- pmod2_cr1$sigma2[2]
pmod3.sigma2.p <- pmod3$sigma2[2]
pmod3.cr0.sigma2.p <- pmod3_cr0$sigma2[2]
pmod3.cr1.sigma2.p <- pmod3_cr1$sigma2[2]
pmod4.sigma2.p <- pmod4$sigma2[2]
pmod4.cr0.sigma2.p <- pmod4_cr0$sigma2[2]
pmod4.cr1.sigma2.p <- pmod4_cr1$sigma2[2]


# sigma.p estimate MSE
pmod1.sigma2.p.mse <- (sigma2.p - pmod1.sigma2.p)^2
pmod1.cr0.sigma2.p.mse <- (sigma2.p - pmod1.cr0.sigma2.p)^2
pmod1.cr1.sigma2.p.mse <- (sigma2.p - pmod1.cr1.sigma2.p)^2
pmod2.sigma2.p.mse <- (sigma2.p - pmod2.sigma2.p)^2
pmod2.cr0.sigma2.p.mse <- (sigma2.p - pmod2.cr0.sigma2.p)^2
pmod2.cr1.sigma2.p.mse <- (sigma2.p - pmod2.cr1.sigma2.p)^2
pmod3.sigma2.p.mse <- (sigma2.p - pmod3.sigma2.p)^2
pmod3.cr0.sigma2.p.mse <- (sigma2.p - pmod3.cr0.sigma2.p)^2
pmod3.cr1.sigma2.p.mse <- (sigma2.p - pmod3.cr1.sigma2.p)^2
pmod4.sigma2.p.mse <- (sigma2.p - pmod4.sigma2.p)^2
pmod4.cr0.sigma2.p.mse <- (sigma2.p - pmod4.cr0.sigma2.p)^2
pmod4.cr1.sigma2.p.mse <- (sigma2.p - pmod4.cr1.sigma2.p)^2


# sigma.s estimate (study level variance)
# in model: # in model: sigma^2.3
pmod1.sigma2.s <- pmod1$sigma2[3]
pmod2.sigma2.s <- pmod2$sigma2[3]
pmod1.cr0.sigma2.s <- pmod1_cr0$sigma2[3]
pmod1.cr1.sigma2.s <- pmod1_cr1$sigma2[3]
pmod2.cr0.sigma2.s <- pmod2_cr0$sigma2[3]
pmod2.cr1.sigma2.s <- pmod2_cr1$sigma2[3]
pmod3.sigma2.s <- pmod3$sigma2[3]
pmod3.cr0.sigma2.s <- pmod3_cr0$sigma2[3]
pmod3.cr1.sigma2.s <- pmod3_cr1$sigma2[3]
pmod4.sigma2.s <- pmod4$sigma2[3]
pmod4.cr0.sigma2.s <- pmod4_cr0$sigma2[3]
pmod4.cr1.sigma2.s <- pmod4_cr1$sigma2[3]


# sigma.s estimate MSE
pmod1.sigma2.s.mse <- (sigma2.s - pmod1.sigma2.s)^2
pmod2.sigma2.s.mse <- (sigma2.s - pmod2.sigma2.s)^2
pmod1.cr0.sigma2.s.mse <- (sigma2.s - pmod1.cr0.sigma2.s)^2
pmod1.cr1.sigma2.s.mse <- (sigma2.s - pmod1.cr1.sigma2.s)^2
pmod2.cr0.sigma2.s.mse <- (sigma2.s - pmod2.cr0.sigma2.s)^2
pmod2.cr1.sigma2.s.mse <- (sigma2.s - pmod2.cr1.sigma2.s)^2
pmod3.sigma2.s.mse <- (sigma2.s - pmod3.sigma2.s)^2
pmod3.cr0.sigma2.s.mse <- (sigma2.s - pmod3.cr0.sigma2.s)^2
pmod3.cr1.sigma2.s.mse <- (sigma2.s - pmod3.cr1.sigma2.s)^2
pmod4.sigma2.s.mse <- (sigma2.s - pmod4.sigma2.s)^2
pmod4.cr0.sigma2.s.mse <- (sigma2.s - pmod4.cr0.sigma2.s)^2
pmod4.cr1.sigma2.s.mse <- (sigma2.s - pmod4.cr1.sigma2.s)^2



# sigma.u estimate (estimate level variance)
# in model: # in model: sigma^2.4
pmod1.sigma2.u <- pmod1$sigma2[4]
pmod2.sigma2.u <- pmod2$sigma2[4]
pmod1.cr0.sigma2.u <- pmod1_cr0$sigma2[4]
pmod1.cr1.sigma2.u <- pmod1_cr1$sigma2[4]
pmod2.cr0.sigma2.u <- pmod2_cr0$sigma2[4]
pmod2.cr1.sigma2.u <- pmod2_cr1$sigma2[4]
pmod3.sigma2.u <- pmod3$sigma2[4]
pmod3.cr0.sigma2.u <- pmod3_cr0$sigma2[4]
pmod3.cr1.sigma2.u <- pmod3_cr1$sigma2[4]
pmod4.sigma2.u <- pmod4$sigma2[4]
pmod4.cr0.sigma2.u <- pmod4_cr0$sigma2[4]
pmod4.cr1.sigma2.u <- pmod4_cr1$sigma2[4]


# sigma.u estimate MSE
pmod1.sigma2.u.mse <- (sigma2.u - pmod1.sigma2.u)^2
pmod2.sigma2.u.mse <- (sigma2.u - pmod2.sigma2.u)^2
pmod1.cr0.sigma2.u.mse <- (sigma2.u - pmod1.cr0.sigma2.u)^2
pmod1.cr1.sigma2.u.mse <- (sigma2.u - pmod1.cr1.sigma2.u)^2
pmod2.cr0.sigma2.u.mse <- (sigma2.u - pmod2.cr0.sigma2.u)^2
pmod2.cr1.sigma2.u.mse <- (sigma2.u - pmod2.cr1.sigma2.u)^2
pmod3.sigma2.u.mse <- (sigma2.u - pmod3.sigma2.u)^2
pmod3.cr0.sigma2.u.mse <- (sigma2.u - pmod3.cr0.sigma2.u)^2
pmod3.cr1.sigma2.u.mse <- (sigma2.u - pmod3.cr1.sigma2.u)^2
pmod4.sigma2.u.mse <- (sigma2.u - pmod4.sigma2.u)^2
pmod4.cr0.sigma2.u.mse <- (sigma2.u - pmod4.cr0.sigma2.u)^2
pmod4.cr1.sigma2.u.mse <- (sigma2.u - pmod4.cr1.sigma2.u)^2



# b1.1 estimate
pmod1.b1.1 <- pmod1$b[2]
pmod1.cr0.b1.1 <- pmod1$b[2]
pmod1.cr1.b1.1 <- pmod1$b[2]
pmod2.b1.1 <- pmod2$b[2]
pmod2.cr0.b1.1 <- pmod2$b[2]
pmod2.cr1.b1.1 <- pmod2$b[2]
pmod3.b1.1 <- pmod3$b[2]
pmod3.cr0.b1.1 <- pmod3$b[2]
pmod3.cr1.b1.1 <- pmod3$b[2]
pmod4.b1.1 <- pmod4$b[2]
pmod4.cr0.b1.1 <- pmod4$b[2]
pmod4.cr1.b1.1 <- pmod4$b[2]


# b1.1 confidence intervals
pmod1.b1.1.ci.lb <- pmod1$ci.lb[2]
pmod1.b1.1.ci.ub <- pmod1$ci.ub[2]
pmod1.cr0.b1.1.ci.lb <- pmod1_cr0$ci.lb[2]
pmod1.cr0.b1.1.ci.ub <- pmod1_cr0$ci.ub[2]
pmod1.cr1.b1.1.ci.lb <- pmod1_cr1$ci.lb[2]
pmod1.cr1.b1.1.ci.ub <- pmod1_cr1$ci.ub[2]
pmod2.b1.1.ci.lb <- pmod2$ci.lb[2]
pmod2.b1.1.ci.ub <- pmod2$ci.ub[2]
pmod2.cr0.b1.1.ci.lb <- pmod2_cr0$ci.lb[2]
pmod2.cr0.b1.1.ci.ub <- pmod2_cr0$ci.ub[2]
pmod2.cr1.b1.1.ci.lb <- pmod2_cr1$ci.lb[2]
pmod2.cr1.b1.1.ci.ub <- pmod2_cr1$ci.ub[2]
pmod3.b1.1.ci.lb <- pmod3$ci.lb[2]
pmod3.b1.1.ci.ub <- pmod3$ci.ub[2]
pmod3.cr0.b1.1.ci.lb <- pmod3_cr0$ci.lb[2]
pmod3.cr0.b1.1.ci.ub <- pmod3_cr0$ci.ub[2]
pmod3.cr1.b1.1.ci.lb <- pmod3_cr1$ci.lb[2]
pmod3.cr1.b1.1.ci.ub <- pmod3_cr1$ci.ub[2]
pmod4.b1.1.ci.lb <- pmod4$ci.lb[2]
pmod4.b1.1.ci.ub <- pmod4$ci.ub[2]
pmod4.cr0.b1.1.ci.lb <- pmod4_cr0$ci.lb[2]
pmod4.cr0.b1.1.ci.ub <- pmod4_cr0$ci.ub[2]
pmod4.cr1.b1.1.ci.lb <- pmod4_cr1$ci.lb[2]
pmod4.cr1.b1.1.ci.ub <- pmod4_cr1$ci.ub[2]



# b1.2 estimate
pmod1.b1.2 <- pmod1$b[3]
pmod1.cr0.b1.2 <- pmod1$b[3]
pmod1.cr1.b1.2 <- pmod1$b[3]
pmod2.b1.2 <- pmod2$b[3]
pmod2.cr0.b1.2 <- pmod2$b[3]
pmod2.cr1.b1.2 <- pmod2$b[3]
pmod3.b1.2 <- pmod3$b[3]
pmod3.cr0.b1.2 <- pmod3$b[3]
pmod3.cr1.b1.2 <- pmod3$b[3]
pmod4.b1.2 <- pmod4$b[3]
pmod4.cr0.b1.2 <- pmod4$b[3]
pmod4.cr1.b1.2 <- pmod4$b[3]


# b1.2 confidence intervals
pmod1.b1.2.ci.lb <- pmod1$ci.lb[3]
pmod1.b1.2.ci.ub <- pmod1$ci.ub[3]
pmod1.cr0.b1.2.ci.lb <- pmod1_cr0$ci.lb[3]
pmod1.cr0.b1.2.ci.ub <- pmod1_cr0$ci.ub[3]
pmod1.cr1.b1.2.ci.lb <- pmod1_cr1$ci.lb[3]
pmod1.cr1.b1.2.ci.ub <- pmod1_cr1$ci.ub[3]
pmod2.b1.2.ci.lb <- pmod2$ci.lb[3]
pmod2.b1.2.ci.ub <- pmod2$ci.ub[3]
pmod2.cr0.b1.2.ci.lb <- pmod2_cr0$ci.lb[3]
pmod2.cr0.b1.2.ci.ub <- pmod2_cr0$ci.ub[3]
pmod2.cr1.b1.2.ci.lb <- pmod2_cr1$ci.lb[3]
pmod2.cr1.b1.2.ci.ub <- pmod2_cr1$ci.ub[3]
pmod3.b1.2.ci.lb <- pmod3$ci.lb[3]
pmod3.b1.2.ci.ub <- pmod3$ci.ub[3]
pmod3.cr0.b1.2.ci.lb <- pmod3_cr0$ci.lb[3]
pmod3.cr0.b1.2.ci.ub <- pmod3_cr0$ci.ub[3]
pmod3.cr1.b1.2.ci.lb <- pmod3_cr1$ci.lb[3]
pmod3.cr1.b1.2.ci.ub <- pmod3_cr1$ci.ub[3]
pmod4.b1.2.ci.lb <- pmod4$ci.lb[3]
pmod4.b1.2.ci.ub <- pmod4$ci.ub[3]
pmod4.cr0.b1.2.ci.lb <- pmod4_cr0$ci.lb[3]
pmod4.cr0.b1.2.ci.ub <- pmod4_cr0$ci.ub[3]
pmod4.cr1.b1.2.ci.lb <- pmod4_cr1$ci.lb[3]
pmod4.cr1.b1.2.ci.ub <- pmod4_cr1$ci.ub[3]



# b2 estimate
pmod1.b2 <- pmod1$b[4]
pmod1.cr0.b2 <- pmod1$b[4]
pmod1.cr1.b2 <- pmod1$b[4]
pmod2.b2 <- pmod2$b[4]
pmod2.cr0.b2 <- pmod2$b[4]
pmod2.cr1.b2 <- pmod2$b[4]
pmod3.b2 <- pmod3$b[4]
pmod3.cr0.b2 <- pmod3$b[4]
pmod3.cr1.b2 <- pmod3$b[4]
pmod4.b2 <- pmod4$b[4]
pmod4.cr0.b2 <- pmod4$b[4]
pmod4.cr1.b2 <- pmod4$b[4]


# b2 confidence intervals
pmod1.b2.ci.lb <- pmod1$ci.lb[4]
pmod1.b2.ci.ub <- pmod1$ci.ub[4]
pmod1.cr0.b2.ci.lb <- pmod1_cr0$ci.lb[4]
pmod1.cr0.b2.ci.ub <- pmod1_cr0$ci.ub[4]
pmod1.cr1.b2.ci.lb <- pmod1_cr1$ci.lb[4]
pmod1.cr1.b2.ci.ub <- pmod1_cr1$ci.ub[4]
pmod2.b2.ci.lb <- pmod2$ci.lb[4]
pmod2.b2.ci.ub <- pmod2$ci.ub[4]
pmod2.cr0.b2.ci.lb <- pmod2_cr0$ci.lb[4]
pmod2.cr0.b2.ci.ub <- pmod2_cr0$ci.ub[4]
pmod2.cr1.b2.ci.lb <- pmod2_cr1$ci.lb[4]
pmod2.cr1.b2.ci.ub <- pmod2_cr1$ci.ub[4]
pmod3.b2.ci.lb <- pmod3$ci.lb[4]
pmod3.b2.ci.ub <- pmod3$ci.ub[4]
pmod3.cr0.b2.ci.lb <- pmod3_cr0$ci.lb[4]
pmod3.cr0.b2.ci.ub <- pmod3_cr0$ci.ub[4]
pmod3.cr1.b2.ci.lb <- pmod3_cr1$ci.lb[4]
pmod3.cr1.b2.ci.ub <- pmod3_cr1$ci.ub[4]
pmod4.b2.ci.lb <- pmod4$ci.lb[4]
pmod4.b2.ci.ub <- pmod4$ci.ub[4]
pmod4.cr0.b2.ci.lb <- pmod4_cr0$ci.lb[4]
pmod4.cr0.b2.ci.ub <- pmod4_cr0$ci.ub[4]
pmod4.cr1.b2.ci.lb <- pmod4_cr1$ci.lb[4]
pmod4.cr1.b2.ci.ub <- pmod4_cr1$ci.ub[4]





# b3 estimate
pmod1.b3 <- pmod1$b[5]
pmod1.cr0.b3 <- pmod1$b[5]
pmod1.cr1.b3 <- pmod1$b[5]
pmod2.b3 <- pmod2$b[5]
pmod2.cr0.b3 <- pmod2$b[5]
pmod2.cr1.b3 <- pmod2$b[5]
pmod3.b3 <- pmod3$b[5]
pmod3.cr0.b3 <- pmod3$b[5]
pmod3.cr1.b3 <- pmod3$b[5]
pmod4.b3 <- pmod4$b[5]
pmod4.cr0.b3 <- pmod4$b[5]
pmod4.cr1.b3 <- pmod4$b[5]


# b3 confidence intervals
pmod1.b3.ci.lb <- pmod1$ci.lb[5]
pmod1.b3.ci.ub <- pmod1$ci.ub[5]
pmod1.cr0.b3.ci.lb <- pmod1_cr0$ci.lb[5]
pmod1.cr0.b3.ci.ub <- pmod1_cr0$ci.ub[5]
pmod1.cr1.b3.ci.lb <- pmod1_cr1$ci.lb[5]
pmod1.cr1.b3.ci.ub <- pmod1_cr1$ci.ub[5]
pmod2.b3.ci.lb <- pmod2$ci.lb[5]
pmod2.b3.ci.ub <- pmod2$ci.ub[5]
pmod2.cr0.b3.ci.lb <- pmod2_cr0$ci.lb[5]
pmod2.cr0.b3.ci.ub <- pmod2_cr0$ci.ub[5]
pmod2.cr1.b3.ci.lb <- pmod2_cr1$ci.lb[5]
pmod2.cr1.b3.ci.ub <- pmod2_cr1$ci.ub[5]
pmod3.b3.ci.lb <- pmod3$ci.lb[5]
pmod3.b3.ci.ub <- pmod3$ci.ub[5]
pmod3.cr0.b3.ci.lb <- pmod3_cr0$ci.lb[5]
pmod3.cr0.b3.ci.ub <- pmod3_cr0$ci.ub[5]
pmod3.cr1.b3.ci.lb <- pmod3_cr1$ci.lb[5]
pmod3.cr1.b3.ci.ub <- pmod3_cr1$ci.ub[5]
pmod4.b3.ci.lb <- pmod4$ci.lb[5]
pmod4.b3.ci.ub <- pmod4$ci.ub[5]
pmod4.cr0.b3.ci.lb <- pmod4_cr0$ci.lb[5]
pmod4.cr0.b3.ci.ub <- pmod4_cr0$ci.ub[5]
pmod4.cr1.b3.ci.lb <- pmod4_cr1$ci.lb[5]
pmod4.cr1.b3.ci.ub <- pmod4_cr1$ci.ub[5]






########### Save results  ----------------------------------------------

# save results 
res <- data.frame(name = rep(name, 12),
                  scenario = rep(scen, 12), 
                  sim = rep(seed, 12),
                  CR_method = rep(c("none", "CR0", "CR1"), each = 4),
                  model = c("PML", "PML-VCV-0.2", "PML-VCV-0.5", "PML-VCV-0.8", 
                            "PML-CR0", "PML-VCV-0.2-CR0", "PML-VCV-0.5-CR0", "PML-VCV-0.8-CR0", 
                            "PML-CR1", "PML-VCV-0.2-CR1", "PML-VCV-0.5-CR1", "PML-VCV-0.8-CR1"), 
                  comp.time = rep(c(pmod1_time, pmod2_time, pmod3_time, pmod4_time), 3),
                  k.studies = rep(k.studies, 12),
                  k.species = rep(k.species, 12),
                  b1.1 = rep(b1.1, 12),
                  b1.2 = rep(b1.2, 12),
                  b2 = rep(b2, 12),
                  b3 = rep(b3, 12),
                  sigma2.n = rep(sigma2.n, 12),
                  sigma2.p = rep(sigma2.p, 12),
                  sigma2.s = rep(sigma2.s, 12),
                  sigma2.u = rep(sigma2.u, 12),
                  mu = rep(mu, 12),
                  rho = rep(rho, 12),
                  mu_est = c(pmod1.est, pmod2.est, pmod3.est, pmod4.est, 
                             pmod1.cr0, pmod2.cr0, pmod3.cr0, pmod4.cr0,
                             pmod1.cr1, pmod2.cr1, pmod3.cr1, pmod4.cr1), 
                  mu_mse = c(pmod1.est.mse, pmod2.est.mse, pmod3.est.mse, pmod4.est.mse, 
                             pmod1.cr0.mse, pmod2.cr0.mse, pmod3.cr0.mse, pmod4.cr0.mse, 
                             pmod1.cr1.mse, pmod2.cr1.mse, pmod3.cr1.mse, pmod4.cr1.mse), 
                  mu_ci_ub = c(pmod1.est.ci.ub, pmod2.est.ci.ub, pmod3.est.ci.ub, pmod4.est.ci.ub, 
                               pmod1.cr0.ci.ub, pmod2.cr0.ci.ub, pmod3.cr0.ci.ub, pmod4.cr0.ci.ub,
                               pmod1.cr1.ci.ub, pmod2.cr1.ci.ub, pmod3.cr1.ci.ub, pmod4.cr1.ci.ub), 
                  mu_ci_lb = c(pmod1.est.ci.lb, pmod2.est.ci.lb, pmod3.est.ci.lb, pmod4.est.ci.lb,
                               pmod1.cr0.ci.lb, pmod2.cr0.ci.lb, pmod3.cr0.ci.lb, pmod4.cr0.ci.lb, 
                               pmod1.cr1.ci.lb, pmod2.cr1.ci.lb, pmod3.cr1.ci.lb, pmod4.cr1.ci.lb),
                  mu_cov = c(pmod1.cov, pmod2.cov, pmod3.cov, pmod4.cov, 
                             pmod1.cr0.cov, pmod2.cr0.cov, pmod3.cr0.cov, pmod4.cr0.cov,
                             pmod1.cr1.cov, pmod2.cr1.cov, pmod3.cr1.cov, pmod4.cr1.cov), 
                  sigma.n_est = c(pmod1.sigma2.n, pmod2.sigma2.n, pmod3.sigma2.n, pmod4.sigma2.n,
                                  pmod1.cr0.sigma2.n, pmod2.cr0.sigma2.n, pmod3.cr0.sigma2.n, pmod4.cr0.sigma2.n,
                                 pmod1.cr1.sigma2.n, pmod2.cr1.sigma2.n, pmod3.cr1.sigma2.n, pmod4.cr1.sigma2.n), 
                  sigma.n_mse = c(pmod1.sigma2.n.mse, pmod2.sigma2.n.mse, pmod3.sigma2.n.mse, pmod4.sigma2.n.mse, 
                                  pmod1.cr0.sigma2.n.mse, pmod2.cr0.sigma2.n.mse, pmod3.cr0.sigma2.n.mse, pmod4.cr0.sigma2.n.mse, 
                                  pmod1.cr1.sigma2.n.mse, pmod2.cr1.sigma2.n.mse, pmod3.cr1.sigma2.n.mse, pmod4.cr1.sigma2.n.mse), 
                  sigma.p_est = c(pmod1.sigma2.p, pmod2.sigma2.p, pmod3.sigma2.p, pmod4.sigma2.p,
                                  pmod1.cr0.sigma2.p, pmod2.cr0.sigma2.p, pmod3.cr0.sigma2.p, pmod4.cr0.sigma2.p, 
                                  pmod1.cr1.sigma2.p, pmod2.cr1.sigma2.p, pmod3.cr1.sigma2.p, pmod4.cr1.sigma2.p), 
                  sigma.p_mse = c(pmod1.sigma2.p.mse, pmod2.sigma2.p.mse, pmod3.sigma2.p.mse, pmod4.sigma2.p.mse, 
                                  pmod1.cr0.sigma2.p.mse, pmod2.cr0.sigma2.p.mse, pmod3.cr0.sigma2.p.mse, pmod4.cr0.sigma2.p.mse, 
                                  pmod1.cr1.sigma2.p.mse, pmod2.cr1.sigma2.p.mse, pmod3.cr1.sigma2.p.mse, pmod4.cr1.sigma2.p.mse), 
                  sigma.s_est = c(pmod1.sigma2.s, pmod2.sigma2.s, pmod3.sigma2.s, pmod4.sigma2.s, 
                                  pmod1.cr0.sigma2.s, pmod2.cr0.sigma2.s, pmod3.cr0.sigma2.s, pmod4.cr0.sigma2.s, 
                                  pmod1.cr1.sigma2.s, pmod2.cr1.sigma2.s, pmod3.cr1.sigma2.s, pmod4.cr1.sigma2.s), 
                  sigma.s_mse = c(pmod1.sigma2.s.mse, pmod2.sigma2.s.mse, pmod3.sigma2.s.mse, pmod4.sigma2.s.mse, 
                                  pmod1.cr0.sigma2.s.mse, pmod2.cr0.sigma2.s.mse, pmod3.cr0.sigma2.s.mse, pmod4.cr0.sigma2.s.mse, 
                                  pmod1.cr1.sigma2.s.mse, pmod2.cr1.sigma2.s.mse, pmod3.cr1.sigma2.s.mse, pmod4.cr1.sigma2.s.mse), 
                  sigma.u_est = c(pmod1.sigma2.u, pmod2.sigma2.u, pmod3.sigma2.u, pmod4.sigma2.u, 
                                  pmod1.cr0.sigma2.u, pmod2.cr0.sigma2.u, pmod3.cr0.sigma2.u, pmod4.cr0.sigma2.u, 
                                  pmod1.cr1.sigma2.u, pmod2.cr1.sigma2.u, pmod3.cr1.sigma2.u, pmod4.cr1.sigma2.u), 
                  sigma.u_mse = c(pmod1.sigma2.u.mse, pmod2.sigma2.u.mse, pmod3.sigma2.u.mse, pmod4.sigma2.u.mse, 
                                  pmod1.cr0.sigma2.u.mse, pmod2.cr0.sigma2.u.mse, pmod3.cr0.sigma2.u.mse, pmod4.cr0.sigma2.u.mse, 
                                  pmod1.cr1.sigma2.u.mse, pmod2.cr1.sigma2.u.mse, pmod3.cr1.sigma2.u.mse, pmod4.cr1.sigma2.u.mse), 
                  b1.1_est = c(pmod1.b1.1, pmod2.b1.1, pmod3.b1.1, pmod4.b1.1,
                             pmod1.cr0.b1.1, pmod2.cr0.b1.1, pmod3.cr0.b1.1, pmod4.cr0.b1.1,
                             pmod1.cr1.b1.1, pmod2.cr1.b1.1, pmod3.cr1.b1.1, pmod4.cr1.b1.1),
                  b1.1_ci_lb = c(pmod1.b1.1.ci.lb, pmod2.b1.1.ci.lb, pmod3.b1.1.ci.lb, pmod4.b1.1.ci.lb,
                               pmod1.cr0.b1.1.ci.lb, pmod2.cr0.b1.1.ci.lb, pmod3.cr0.b1.1.ci.lb, pmod4.cr0.b1.1.ci.lb,
                               pmod1.cr1.b1.1.ci.lb, pmod2.cr1.b1.1.ci.lb, pmod3.cr1.b1.1.ci.lb, pmod4.cr1.b1.1.ci.lb),
                  b1.1_ci_ub = c(pmod1.b1.1.ci.ub, pmod2.b1.1.ci.ub, pmod3.b1.1.ci.ub, pmod4.b1.1.ci.ub,
                               pmod1.cr0.b1.1.ci.ub, pmod2.cr0.b1.1.ci.ub, pmod3.cr0.b1.1.ci.ub, pmod4.cr0.b1.1.ci.ub,
                               pmod1.cr1.b1.1.ci.ub, pmod2.cr1.b1.1.ci.ub, pmod3.cr1.b1.1.ci.ub, pmod4.cr1.b1.1.ci.ub),
                  b1.2_est = c(pmod1.b1.2, pmod2.b1.2, pmod3.b1.2, pmod4.b1.2,
                               pmod1.cr0.b1.2, pmod2.cr0.b1.2, pmod3.cr0.b1.2, pmod4.cr0.b1.2,
                               pmod1.cr1.b1.2, pmod2.cr1.b1.2, pmod3.cr1.b1.2, pmod4.cr1.b1.2),
                  b1.2_ci_lb = c(pmod1.b1.2.ci.lb, pmod2.b1.2.ci.lb, pmod3.b1.2.ci.lb, pmod4.b1.2.ci.lb,
                                 pmod1.cr0.b1.2.ci.lb, pmod2.cr0.b1.2.ci.lb, pmod3.cr0.b1.2.ci.lb, pmod4.cr0.b1.2.ci.lb,
                                 pmod1.cr1.b1.2.ci.lb, pmod2.cr1.b1.2.ci.lb, pmod3.cr1.b1.2.ci.lb, pmod4.cr1.b1.2.ci.lb),
                  b1.2_ci_ub = c(pmod1.b1.2.ci.ub, pmod2.b1.2.ci.ub, pmod3.b1.2.ci.ub, pmod4.b1.2.ci.ub,
                                 pmod1.cr0.b1.2.ci.ub, pmod2.cr0.b1.2.ci.ub, pmod3.cr0.b1.2.ci.ub, pmod4.cr0.b1.2.ci.ub,
                                 pmod1.cr1.b1.2.ci.ub, pmod2.cr1.b1.2.ci.ub, pmod3.cr1.b1.2.ci.ub, pmod4.cr1.b1.2.ci.ub),
                  b2_est = c(pmod1.b2, pmod2.b2, pmod3.b2, pmod4.b2,
                               pmod1.cr0.b2, pmod2.cr0.b2, pmod3.cr0.b2, pmod4.cr0.b2,
                               pmod1.cr1.b2, pmod2.cr1.b2, pmod3.cr1.b2, pmod4.cr1.b2),
                  b2_ci_lb = c(pmod1.b2.ci.lb, pmod2.b2.ci.lb, pmod3.b2.ci.lb, pmod4.b2.ci.lb,
                                 pmod1.cr0.b2.ci.lb, pmod2.cr0.b2.ci.lb, pmod3.cr0.b2.ci.lb, pmod4.cr0.b2.ci.lb,
                                 pmod1.cr1.b2.ci.lb, pmod2.cr1.b2.ci.lb, pmod3.cr1.b2.ci.lb, pmod4.cr1.b2.ci.lb),
                  b2_ci_ub = c(pmod1.b2.ci.ub, pmod2.b2.ci.ub, pmod3.b2.ci.ub, pmod4.b2.ci.ub,
                                 pmod1.cr0.b2.ci.ub, pmod2.cr0.b2.ci.ub, pmod3.cr0.b2.ci.ub, pmod4.cr0.b2.ci.ub,
                                 pmod1.cr1.b2.ci.ub, pmod2.cr1.b2.ci.ub, pmod3.cr1.b2.ci.ub, pmod4.cr1.b2.ci.ub),
                  b3_est = c(pmod1.b3, pmod2.b3, pmod3.b3, pmod4.b3,
                               pmod1.cr0.b3, pmod2.cr0.b3, pmod3.cr0.b3, pmod4.cr0.b3,
                               pmod1.cr1.b3, pmod2.cr1.b3, pmod3.cr1.b3, pmod4.cr1.b3),
                  b3_ci_lb = c(pmod1.b3.ci.lb, pmod2.b3.ci.lb, pmod3.b3.ci.lb, pmod4.b3.ci.lb,
                                 pmod1.cr0.b3.ci.lb, pmod2.cr0.b3.ci.lb, pmod3.cr0.b3.ci.lb, pmod4.cr0.b3.ci.lb,
                                 pmod1.cr1.b3.ci.lb, pmod2.cr1.b3.ci.lb, pmod3.cr1.b3.ci.lb, pmod4.cr1.b3.ci.lb),
                  b3_ci_ub = c(pmod1.b3.ci.ub, pmod2.b3.ci.ub, pmod3.b3.ci.ub, pmod4.b3.ci.ub,
                               pmod1.cr0.b3.ci.ub, pmod2.cr0.b3.ci.ub, pmod3.cr0.b3.ci.ub, pmod4.cr0.b3.ci.ub,
                               pmod1.cr1.b3.ci.ub, pmod2.cr1.b3.ci.ub, pmod3.cr1.b3.ci.ub, pmod4.cr1.b3.ci.ub)
)


# save the output according to the job array:
save(list = "res", file = paste0("results/res_", job, ".RDATA"))
