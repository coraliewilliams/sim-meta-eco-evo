###############################################################################
# (STUDY 2) SIMULATION OF META-ANALYSIS MODELS
###############################################################################

rm(list=ls())

# Load libraries
library(pacman) # checks if package is installed, if not installs it
p_load(metafor, MASS, clubSandwich, ape)


########## Load parameters conditions ---------------------------------------

### load job array
tab <- read.csv("job_array_study2.csv")
# tab <-  scen.tab2

### job number from pbs script
job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#job <- 1717 # for testing locally

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
tree <- rtree(k.species, tip.label=seq_len(k.species))
tree <- compute.brlen(tree, power = 1) # based on Cinar et al., 2022
P <- vcv(tree, corr=TRUE)
P <- P[order(as.numeric(rownames(P))), order(as.numeric(rownames(P)))]

### when the power parameter is 0, vcv inserts Nan to the off-diagonals. We
### want the correlation matrix to be an identity matrix in these conditions.
### So, we will replace Nans with 0s.
#P[is.nan(P)] <- 0

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

yi <- mu + u.u + u.s + u.n + u.p + mi


# get simulated data
dat <- data.frame(name = name, scenario = scen, seed = seed, job_number = job,
                  study = study, id = id, esid = esid, yi = yi, vi = vi,
                  u.u = u.u, u.s = u.s, u.n = u.n, u.p = u.p,  mi = mi)
# save simulated data in R file
save(list = "dat", file = paste0("data/simdat_", job, ".RDATA"))


########### Run model 1: phylogenetic multilevel model  ----------------------------------------------

# run model
ptm <- proc.time()
pmod1 <- try(rma.mv(yi,
                   vi,
                   random = list(~ 1 | species, ~ 1 | species.phylo, ~ 1 | study, ~ 1| id),
                   R = list(species.phylo=P),
                   test = "t",
                   dfs = k.studies - 1), #default "residuals" (k-1), other option "contain" (k.studies-1)
            silent = TRUE)
pmod1_time <- (proc.time() - ptm)[3]

# if (inherits(pmod1, "try-error")){
#   system(paste0("echo 'PROBLEM (l = ", l, ", m = ", m, "): model 3 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   pmod1.comp.time[l] <- mod3_time
# }


pmod1_cr0 <- robust(pmod1, cluster = study, adjust=FALSE)
pmod1_cr1 <- robust(pmod1, cluster = study, adjust=TRUE)



######## testing: degrees of freedom

# using residuals (k-1) = "residuals"
# pmod.res <- rma.mv(yi,
#                     vi,
#                     random = list(~ 1 | species, ~ 1 | species.phylo, ~ 1 | study/id),
#                     R = list(species.phylo=P),
#                     test = "t",
#                     dfs = "residual") #default , other option "contain" (k.studies-1)
# 
# # using contains (k.studies) = "contain"
# pmod.contain <- rma.mv(yi,
#                 vi,
#                 random = list(~ 1 | species, ~ 1 | species.phylo, ~ 1 | study/id),
#                 R = list(species.phylo=P),
#                 test = "t",
#                 dfs = "contain") #default "residuals" (k-1), other option "contain" (k.studies-1)
# 




########### Run model 2: phylogenetic multilevel model with sampling VCV (rho.hat=0.2) -------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.2)

ptm <- proc.time()
pmod2 <- try(rma.mv(yi,
                   V=V,
                   random = list(~ 1 | species, ~ 1 | species.phylo, ~ 1 | study, ~ 1| id),
                   R = list(species.phylo=P),
                   test = "t",
                   dfs = k.studies - 1), 
            silent = TRUE)
pmod2_time <- (proc.time() - ptm)[3]

# if (inherits(pmod2, "try-error")){
#   system(paste0("echo 'PROBLEM (job = ", job, ", m = ", m, "): model 4 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   pmod2.comp.time[l] <- pmod2_time
# }

pmod2_cr0 <- robust(pmod2, cluster = study, adjust=FALSE)
pmod2_cr1 <- robust(pmod2, cluster = study, adjust=TRUE)




########### Run model 3: phylogenetic multilevel model with sampling VCV (rho.hat=0.5) -------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.5)

ptm <- proc.time()
pmod3 <- try(rma.mv(yi,
                    V=V,
                    random = list(~ 1 | species, ~ 1 | species.phylo, ~ 1 | study, ~ 1| id),
                    R = list(species.phylo=P),
                    test = "t",
                    dfs = k.studies - 1), 
             silent = TRUE)
pmod3_time <- (proc.time() - ptm)[3]

# if (inherits(pmod2, "try-error")){
#   system(paste0("echo 'PROBLEM (job = ", job, ", m = ", m, "): model 4 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   pmod2.comp.time[l] <- pmod2_time
# }

pmod3_cr0 <- robust(pmod3, cluster = study, adjust=FALSE)
pmod3_cr1 <- robust(pmod3, cluster = study, adjust=TRUE)





########### Run model 4: phylogenetic multilevel model with sampling VCV (rho.hat=0.8) -------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.8)

ptm <- proc.time()
pmod4 <- try(rma.mv(yi,
                    V=V,
                    random = list(~ 1 | species, ~ 1 | species.phylo, ~ 1 | study, ~ 1| id),
                    R = list(species.phylo=P),
                    test = "t",
                    dfs = k.studies - 1), 
             silent = TRUE)
pmod4_time <- (proc.time() - ptm)[3]

# if (inherits(pmod2, "try-error")){
#   system(paste0("echo 'PROBLEM (job = ", job, ", m = ", m, "): model 4 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   pmod2.comp.time[l] <- pmod2_time
# }

pmod4_cr0 <- robust(pmod4, cluster = study, adjust=FALSE)
pmod4_cr1 <- robust(pmod4, cluster = study, adjust=TRUE)



########### Run model 5: multilevel model with sampling VCV (rho.hat=0.5) -------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.5)

ptm <- proc.time()
pmod5 <- try(rma.mv(yi,
                    V=V,
                    random = list(~ 1 | species, ~ 1 | study, ~ 1| id),
                    test = "t",
                    dfs = k.studies - 1), 
             silent = TRUE)
pmod5_time <- (proc.time() - ptm)[3]

# if (inherits(pmod2, "try-error")){
#   system(paste0("echo 'PROBLEM (job = ", job, ", m = ", m, "): model 4 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   pmod2.comp.time[l] <- pmod2_time
# }

pmod5_cr0 <- robust(pmod5, cluster = study, adjust=FALSE)
pmod5_cr1 <- robust(pmod5, cluster = study, adjust=TRUE)





########### Extract model estimates  -------------------------------------------------------------------

# mu estimate 
pmod1.est <- pmod1$b
pmod1.cr0 <- pmod1_cr0$b
pmod1.cr1 <- pmod1_cr1$b
pmod2.est <- pmod2$b
pmod2.cr0 <- pmod2_cr0$b
pmod2.cr1 <- pmod2_cr1$b
pmod3.est <- pmod3$b
pmod3.cr0 <- pmod3_cr0$b
pmod3.cr1 <- pmod3_cr1$b
pmod4.est <- pmod4$b
pmod4.cr0 <- pmod4_cr0$b
pmod4.cr1 <- pmod4_cr1$b
pmod5.est <- pmod5$b
pmod5.cr0 <- pmod5_cr0$b
pmod5.cr1 <- pmod5_cr1$b

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
pmod5.est.mse <- (mu - pmod5.est)^2
pmod5.cr0.mse <- (mu - pmod5.cr0)^2
pmod5.cr1.mse <- (mu - pmod5.cr1)^2

# mu confidence interval
pmod1.est.ci.ub <- pmod1$ci.ub
pmod1.est.ci.lb <- pmod1$ci.lb
pmod2.est.ci.ub <- pmod2$ci.ub
pmod2.est.ci.lb <- pmod2$ci.lb
pmod1.cr0.ci.ub <- pmod1_cr0$ci.ub
pmod1.cr0.ci.lb <- pmod1_cr0$ci.lb
pmod1.cr1.ci.ub <- pmod1_cr1$ci.ub
pmod1.cr1.ci.lb <- pmod1_cr1$ci.lb
pmod2.cr0.ci.ub <- pmod2_cr0$ci.ub
pmod2.cr0.ci.lb <- pmod2_cr0$ci.lb
pmod2.cr1.ci.ub <- pmod2_cr1$ci.ub
pmod2.cr1.ci.lb <- pmod2_cr1$ci.lb
pmod3.est.ci.ub <- pmod3$ci.ub
pmod3.est.ci.lb <- pmod3$ci.lb
pmod3.cr0.ci.ub <- pmod3_cr0$ci.ub
pmod3.cr0.ci.lb <- pmod3_cr0$ci.lb
pmod3.cr1.ci.ub <- pmod3_cr1$ci.ub
pmod3.cr1.ci.lb <- pmod3_cr1$ci.lb
pmod4.est.ci.ub <- pmod4$ci.ub
pmod4.est.ci.lb <- pmod4$ci.lb
pmod4.cr0.ci.ub <- pmod4_cr0$ci.ub
pmod4.cr0.ci.lb <- pmod4_cr0$ci.lb
pmod4.cr1.ci.ub <- pmod4_cr1$ci.ub
pmod4.cr1.ci.lb <- pmod4_cr1$ci.lb
pmod5.est.ci.ub <- pmod5$ci.ub
pmod5.est.ci.lb <- pmod5$ci.lb
pmod5.cr0.ci.ub <- pmod5_cr0$ci.ub
pmod5.cr0.ci.lb <- pmod5_cr0$ci.lb
pmod5.cr1.ci.ub <- pmod5_cr1$ci.ub
pmod5.cr1.ci.lb <- pmod5_cr1$ci.lb


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
pmod5.cov <- pmod5.est.ci.lb < mu && pmod5.est.ci.ub > mu
pmod5.cr0.cov <- pmod5.cr0.ci.lb < mu && pmod5.cr0.ci.ub > mu
pmod5.cr1.cov <- pmod5.cr1.ci.lb < mu && pmod5.cr1.ci.ub > mu


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
## ML model
pmod5.sigma2.n <- pmod5$sigma2[1]
pmod5.cr0.sigma2.n <- pmod5_cr0$sigma2[1]
pmod5.cr1.sigma2.n <- pmod5_cr1$sigma2[1]


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
pmod5.sigma2.n.mse <- (sigma2.n - pmod5.sigma2.n)^2
pmod5.cr0.sigma2.n.mse <- (sigma2.n - pmod5.cr0.sigma2.n)^2
pmod5.cr1.sigma2.n.mse <- (sigma2.n - pmod5.cr1.sigma2.n)^2


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
## ML model has no phylogenetic random effect
pmod5.sigma2.p <- NA
pmod5.cr0.sigma2.p <- NA
pmod5.cr1.sigma2.p <- NA

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
## ML model has no phylogenetic random effect
pmod5.sigma2.p.mse <- NA
pmod5.cr0.sigma2.p.mse <- NA
pmod5.cr1.sigma2.p.mse <- NA


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
## ML model 
pmod5.sigma2.s <- pmod5$sigma2[2]
pmod5.cr0.sigma2.s <- pmod5_cr0$sigma2[2]
pmod5.cr1.sigma2.s <- pmod5_cr1$sigma2[2]

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
pmod5.sigma2.s.mse <- (sigma2.s - pmod5.sigma2.s)^2
pmod5.cr0.sigma2.s.mse <- (sigma2.s - pmod5.cr0.sigma2.s)^2
pmod5.cr1.sigma2.s.mse <- (sigma2.s - pmod5.cr1.sigma2.s)^2


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
## ML model 
pmod5.sigma2.u <- pmod5$sigma2[3]
pmod5.cr0.sigma2.u <- pmod5_cr0$sigma2[3]
pmod5.cr1.sigma2.u <- pmod5_cr1$sigma2[3]


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
pmod5.sigma2.u.mse <- (sigma2.u - pmod5.sigma2.u)^2
pmod5.cr0.sigma2.u.mse <- (sigma2.u - pmod5.cr0.sigma2.u)^2
pmod5.cr1.sigma2.u.mse <- (sigma2.u - pmod5.cr1.sigma2.u)^2



########### Save results  ----------------------------------------------

# save results 
res <- data.frame(name = rep(name, 15),
                  scenario = rep(scen, 15), 
                  sim = rep(seed, 15),
                  CR_method = rep(c("none", "CR0", "CR1"), each = 5),
                  model = c("PML", "PML-VCV-02", "PML-VCV-05", "PML-VCV-08", "ML-VCV-0.5",
                            "PML-CR0", "PML-VCV-0.2-CR0", "PML-VCV-0.5-CR0", "PML-VCV-0.8-CR0", "ML-VCV-0.5-CR0",
                            "PML-CR1", "PML-VCV-0.2-CR1", "PML-VCV-0.5-CR1", "PML-VCV-0.8-CR1","ML-VCV-0.5-CR1"), 
                  comp.time = rep(c(pmod1_time, pmod2_time, pmod3_time, pmod4_time, pmod5_time), 3),
                  mu_est = c(pmod1.est, pmod2.est, pmod3.est, pmod4.est, pmod5.est,
                             pmod1.cr0, pmod2.cr0, pmod3.cr0, pmod4.cr0, pmod5.cr0,
                             pmod1.cr1, pmod2.cr1, pmod3.cr1, pmod4.cr1, pmod5.cr1), 
                  mu_mse = c(pmod1.est.mse, pmod2.est.mse, pmod3.est.mse, pmod4.est.mse, pmod5.est.mse,
                             pmod1.cr0.mse, pmod2.cr0.mse, pmod3.cr0.mse, pmod4.cr0.mse, pmod5.cr0.mse,
                             pmod1.cr1.mse, pmod2.cr1.mse, pmod3.cr1.mse, pmod4.cr1.mse, pmod5.cr1.mse), 
                  mu_ci_ub = c(pmod1.est.ci.ub, pmod2.est.ci.ub, pmod3.est.ci.ub, pmod4.est.ci.ub, pmod5.est.ci.ub,
                               pmod1.cr0.ci.ub, pmod2.cr0.ci.ub, pmod3.cr0.ci.ub, pmod4.cr0.ci.ub, pmod5.cr0.ci.ub,
                               pmod1.cr1.ci.ub, pmod2.cr1.ci.ub, pmod3.cr1.ci.ub, pmod4.cr1.ci.ub, pmod5.cr1.ci.ub), 
                  mu_ci_lb = c(pmod1.est.ci.lb, pmod2.est.ci.lb, pmod3.est.ci.lb, pmod4.est.ci.lb, pmod5.est.ci.lb,
                               pmod1.cr0.ci.lb, pmod2.cr0.ci.lb, pmod3.cr0.ci.lb, pmod4.cr0.ci.lb, pmod5.cr0.ci.lb,
                               pmod1.cr1.ci.lb, pmod2.cr1.ci.lb, pmod3.cr1.ci.lb, pmod4.cr1.ci.lb, pmod5.cr1.ci.lb),
                  mu_cov = c(pmod1.cov, pmod2.cov, pmod3.cov, pmod4.cov, pmod5.cov,
                             pmod1.cr0.cov, pmod2.cr0.cov, pmod3.cr0.cov, pmod4.cr0.cov, pmod5.cr0.cov,
                             pmod1.cr1.cov, pmod2.cr1.cov, pmod3.cr1.cov, pmod4.cr1.cov, pmod5.cr1.cov), 
                  sigma.n_est = c(pmod1.sigma2.n, pmod2.sigma2.n, pmod3.sigma2.n, pmod4.sigma2.n, pmod5.sigma2.n,
                                  pmod1.cr0.sigma2.n, pmod2.cr0.sigma2.n, pmod3.cr0.sigma2.n, pmod4.cr0.sigma2.n, pmod5.cr0.sigma2.n,
                                 pmod1.cr1.sigma2.n, pmod2.cr1.sigma2.n, pmod3.cr1.sigma2.n, pmod4.cr1.sigma2.n, pmod5.sigma2.n), 
                  sigma.n_mse = c(pmod1.sigma2.n.mse, pmod2.sigma2.n.mse, pmod3.sigma2.n.mse, pmod4.sigma2.n.mse, pmod5.sigma2.n.mse,
                                  pmod1.cr0.sigma2.n.mse, pmod2.cr0.sigma2.n.mse, pmod3.cr0.sigma2.n.mse, pmod4.cr0.sigma2.n.mse, pmod5.cr0.sigma2.n.mse,
                                  pmod1.cr1.sigma2.n.mse, pmod2.cr1.sigma2.n.mse, pmod3.cr1.sigma2.n.mse, pmod4.cr1.sigma2.n.mse, pmod5.sigma2.n.mse), 
                  sigma.p_est = c(pmod1.sigma2.p, pmod2.sigma2.p, pmod3.sigma2.p, pmod4.sigma2.p, pmod5.sigma2.p,
                                  pmod1.cr0.sigma2.p, pmod2.cr0.sigma2.p, pmod3.cr0.sigma2.p, pmod4.cr0.sigma2.p, pmod5.cr0.sigma2.p,
                                  pmod1.cr1.sigma2.p, pmod2.cr1.sigma2.p, pmod3.cr1.sigma2.p, pmod4.cr1.sigma2.p, pmod5.sigma2.p), 
                  sigma.p_mse = c(pmod1.sigma2.p.mse, pmod2.sigma2.p.mse, pmod3.sigma2.p.mse, pmod4.sigma2.p.mse, pmod5.sigma2.p.mse,
                                  pmod1.cr0.sigma2.p.mse, pmod2.cr0.sigma2.p.mse, pmod3.cr0.sigma2.p.mse, pmod4.cr0.sigma2.p.mse, pmod5.cr0.sigma2.p.mse,
                                  pmod1.cr1.sigma2.p.mse, pmod2.cr1.sigma2.p.mse, pmod3.cr1.sigma2.p.mse, pmod4.cr1.sigma2.p.mse, pmod5.sigma2.p.mse), 
                  sigma.s_est = c(pmod1.sigma2.s, pmod2.sigma2.s, pmod3.sigma2.s, pmod4.sigma2.s, pmod5.sigma2.s,
                                  pmod1.cr0.sigma2.s, pmod2.cr0.sigma2.s, pmod3.cr0.sigma2.s, pmod4.cr0.sigma2.s, pmod5.sigma2.s,
                                  pmod1.cr1.sigma2.s, pmod2.cr1.sigma2.s, pmod3.cr1.sigma2.s, pmod4.cr1.sigma2.s, pmod5.sigma2.s), 
                  sigma.s_mse = c(pmod1.sigma2.s.mse, pmod2.sigma2.s.mse, pmod3.sigma2.s.mse, pmod4.sigma2.s.mse, pmod5.sigma2.s.mse,
                                  pmod1.cr0.sigma2.s.mse, pmod2.cr0.sigma2.s.mse, pmod3.cr0.sigma2.s.mse, pmod4.cr0.sigma2.s.mse, pmod5.sigma2.s.mse,
                                  pmod1.cr1.sigma2.s.mse, pmod2.cr1.sigma2.s.mse, pmod3.cr1.sigma2.s.mse, pmod4.cr1.sigma2.s.mse, pmod5.sigma2.s.mse), 
                  sigma.u_est = c(pmod1.sigma2.u, pmod2.sigma2.u, pmod3.sigma2.u, pmod4.sigma2.u, pmod5.sigma2.u,
                                  pmod1.cr0.sigma2.u, pmod2.cr0.sigma2.u, pmod3.cr0.sigma2.u, pmod4.cr0.sigma2.u, pmod5.sigma2.u,
                                  pmod1.cr1.sigma2.u, pmod2.cr1.sigma2.u, pmod3.cr1.sigma2.u, pmod4.cr1.sigma2.u, pmod5.sigma2.u), 
                  sigma.u_mse = c(pmod1.sigma2.u.mse, pmod2.sigma2.u.mse, pmod3.sigma2.u.mse, pmod4.sigma2.u.mse, pmod5.sigma2.u.mse,
                                  pmod1.cr0.sigma2.u.mse, pmod2.cr0.sigma2.u.mse, pmod3.cr0.sigma2.u.mse, pmod4.cr0.sigma2.u.mse, pmod5.sigma2.u.mse,
                                  pmod1.cr1.sigma2.u.mse, pmod2.cr1.sigma2.u.mse, pmod3.cr1.sigma2.u.mse, pmod4.cr1.sigma2.u.mse, pmod5.sigma2.u.mse), 
                  k.studies = rep(k.studies, 15),
                  sigma2.n = rep(sigma2.n, 15),
                  sigma2.p = rep(sigma2.p, 15),
                  sigma2.s = rep(sigma2.s, 15),
                  sigma2.u = rep(sigma2.u, 15),
                  mu = rep(mu, 15),
                  rho = rep(rho, 15)
)


# save the output according to the job array:
save(list = "res", file = paste0("results/res_", job, ".RDATA"))
