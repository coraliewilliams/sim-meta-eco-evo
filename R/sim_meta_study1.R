###############################################################################
# (STUDY 1) SIMULATION OF META-ANALYSIS MODELS
###############################################################################

rm(list=ls())

# Load libraries
library(pacman) # checks if package is installed, if not installs it
p_load(metafor, MASS, clubSandwich)



########## Load parameters conditions ---------------------------------------

### load job array
tab <- read.csv("job_array_study1.csv")
#tab <- scen.tab1

### job number from pbs script
job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#job <- 359999 # for testing
  
### parameters for current job
name <- tab$name[tab$job_number == job] 
scen <- tab$scenario[tab$job_number == job] 
seed <- tab$sim[tab$job_number == job]
k.studies <- tab$k.studies[tab$job_number == job]  
mu <- tab$mu[tab$job_number == job]
sigma2.u <- tab$sigma2.u[tab$job_number == job]   
sigma2.s <- tab$sigma2.s[tab$job_number == job]       
rho <- tab$rho[tab$job_number == job]    
rho.hat <- tab$rho.hat[tab$job_number == job]

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

### simulate random effects
u.u <- rnorm(k, 0, sqrt(sigma2.u))
u.s <- rnorm(k.studies, 0, sqrt(sigma2.s))[study]

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

yi <- mu + u.u + u.s + mi

# get simulated data
dat <- data.frame(name = name, scenario = scen, seed = seed, job_number = job,
                  study = study, id = id, esid = esid,
                  yi = yi, vi = vi, u.u = u.u, u.s = u.s, mi = mi)
# save simulated data in R file
save(list = "dat", file = paste0("data/simdat_", job, ".RDATA"))


########### Run model 1: fixed effect  ------------------------------------------------------------------

ptm <- proc.time()
mod1 <- try(rma.mv(yi, 
                   vi,
                   test = "t",
                   dfs = k.studies - 1), 
            silent = TRUE)
mod1_time <- (proc.time() - ptm)[3]

mod1_cr0 <- robust(mod1, cluster = study, adjust=FALSE)
mod1_cr1 <- robust(mod1, cluster = study, adjust=TRUE)
mod1_cr2 <- robust(mod1, cluster = study, clubSandwich=TRUE)



########### Run model 2: simple random effect  ------------------------------------------------------

# run model
ptm <- proc.time()
mod2 <- try(rma.mv(yi,
                   vi,
                   random = list(~ 1 | id),
                   test = "t",
                   dfs = k.studies - 1),
            silent = TRUE)
mod2_time <- (proc.time() - ptm)[3]

# if (inherits(mod2, "try-error")){
#   system(paste0("echo 'PROBLEM (l = ", l, ", m = ", m, "): model 2 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   mod2.comp.time[l] <- mod2_time
# }


mod2_cr0 <- robust(mod2, cluster = study, adjust=FALSE)
mod2_cr1 <- robust(mod2, cluster = study, adjust=TRUE)
mod2_cr2 <- robust(mod2, cluster = study, clubSandwich=TRUE)




########### Run model 3: multilevel model  ----------------------------------------------------------

# run model
ptm <- proc.time()
mod3 <- try(rma.mv(yi,
                   vi,
                   random = list(~ 1 | study, ~ 1 | id),
                   test = "t",
                   dfs = k.studies - 1),
            silent = TRUE)
mod3_time <- (proc.time() - ptm)[3]

# if (inherits(mod2, "try-error")){
#   system(paste0("echo 'PROBLEM (l = ", l, ", m = ", m, "): model 3 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   mod3.comp.time[l] <- mod3_time
# }


mod3_cr0 <- robust(mod3, cluster = study, adjust=FALSE)
mod3_cr1 <- robust(mod3, cluster = study, adjust=TRUE)
mod3_cr2 <- robust(mod3, cluster = study, clubSandwich=TRUE)



########### Run model 4: multilevel model with sampling VCV  -----------------------------------------

# get V matrix based on rho.hat 
V <- vcalc(vi, cluster=study, obs=esid, rho=rho.hat)

ptm <- proc.time()
mod4 <- try(rma.mv(yi,
                   V=V,
                   random = list(~ 1 | study, ~ 1 | id),
                   test = "t",
                   dfs = k.studies - 1),
            silent = TRUE)
mod4_time <- (proc.time() - ptm)[3]

# if (inherits(mod4, "try-error")){
#   system(paste0("echo 'PROBLEM (job = ", job, ", m = ", m, "): model 4 did not fit' >> ./logs/outputs_scen", j, ".txt"))
#   next
# } else {
#   mod4.comp.time[l] <- mod4_time
# }

mod4_cr0 <- robust(mod4, cluster = study, adjust=FALSE)
mod4_cr1 <- robust(mod4, cluster = study, adjust=TRUE)
mod4_cr2 <- robust(mod4, cluster = study, clubSandwich=TRUE)





########### Extract model estimates  -------------------------------------------------------------------

# mu estimate 
mod1.est <- mod1$b
mod2.est <- mod2$b
mod3.est <- mod3$b
mod4.est <- mod4$b
mod1.cr0 <- mod1_cr0$b
mod1.cr1 <- mod1_cr1$b
mod1.cr2 <- mod1_cr2$b
mod2.cr0 <- mod2_cr0$b
mod2.cr1 <- mod2_cr1$b
mod2.cr2 <- mod2_cr2$b
mod3.cr0 <- mod3_cr0$b
mod3.cr1 <- mod3_cr1$b
mod3.cr2 <- mod3_cr2$b
mod4.cr0 <- mod4_cr0$b
mod4.cr1 <- mod4_cr1$b
mod4.cr2 <- mod4_cr2$b

# mu estimate MSE
mod1.est.mse <- (mu - mod1.est)^2  
mod2.est.mse <- (mu - mod2.est)^2
mod3.est.mse <- (mu - mod3.est)^2
mod4.est.mse <- (mu - mod4.est)^2
mod1.cr0.mse <- (mu - mod1.cr0)^2
mod1.cr1.mse <- (mu - mod1.cr1)^2
mod1.cr2.mse <- (mu - mod1.cr2)^2
mod2.cr0.mse <- (mu - mod2.cr0)^2
mod2.cr1.mse <- (mu - mod2.cr1)^2
mod2.cr2.mse <- (mu - mod2.cr2)^2
mod3.cr0.mse <- (mu - mod3.cr0)^2
mod3.cr1.mse <- (mu - mod3.cr1)^2
mod3.cr2.mse <- (mu - mod3.cr2)^2
mod4.cr0.mse <- (mu - mod4.cr0)^2
mod4.cr1.mse <- (mu - mod4.cr1)^2
mod4.cr2.mse <- (mu - mod4.cr2)^2


# mu confidence interval
mod1.est.ci.ub <- mod1$ci.ub
mod1.est.ci.lb <- mod1$ci.lb
mod2.est.ci.ub <- mod2$ci.ub
mod2.est.ci.lb <- mod2$ci.lb
mod3.est.ci.ub <- mod3$ci.ub
mod3.est.ci.lb <- mod3$ci.lb
mod4.est.ci.ub <- mod4$ci.ub
mod4.est.ci.lb <- mod4$ci.lb
mod1.cr0.ci.ub <- mod1_cr0$ci.ub
mod1.cr0.ci.lb <- mod1_cr0$ci.lb
mod1.cr1.ci.ub <- mod1_cr1$ci.ub
mod1.cr1.ci.lb <- mod1_cr1$ci.lb
mod1.cr2.ci.ub <- mod1_cr2$ci.ub
mod1.cr2.ci.lb <- mod1_cr2$ci.lb
mod2.cr0.ci.ub <- mod2_cr0$ci.ub
mod2.cr0.ci.lb <- mod2_cr0$ci.lb
mod2.cr1.ci.ub <- mod2_cr1$ci.ub
mod2.cr1.ci.lb <- mod2_cr1$ci.lb
mod2.cr2.ci.ub <- mod2_cr2$ci.ub
mod2.cr2.ci.lb <- mod2_cr2$ci.lb
mod3.cr0.ci.ub <- mod3_cr0$ci.ub
mod3.cr0.ci.lb <- mod3_cr0$ci.lb
mod3.cr1.ci.ub <- mod3_cr1$ci.ub
mod3.cr1.ci.lb <- mod3_cr1$ci.lb
mod3.cr2.ci.ub <- mod3_cr2$ci.ub
mod3.cr2.ci.lb <- mod3_cr2$ci.lb
mod4.cr0.ci.ub <- mod4_cr0$ci.ub
mod4.cr0.ci.lb <- mod4_cr0$ci.lb
mod4.cr1.ci.ub <- mod4_cr1$ci.ub
mod4.cr1.ci.lb <- mod4_cr1$ci.lb
mod4.cr2.ci.ub <- mod4_cr2$ci.ub
mod4.cr2.ci.lb <- mod4_cr2$ci.lb


# mu coverage
mod1.cov <- mod1.est.ci.lb < mu && mod1.est.ci.ub > mu
mod2.cov <- mod2.est.ci.lb < mu && mod2.est.ci.ub > mu
mod3.cov <- mod3.est.ci.lb < mu && mod3.est.ci.ub > mu
mod4.cov <- mod4.est.ci.lb < mu && mod4.est.ci.ub > mu
mod1.cr0.cov <- mod1.cr0.ci.lb < mu && mod1.cr0.ci.ub > mu
mod1.cr1.cov <- mod1.cr1.ci.lb < mu && mod1.cr1.ci.ub > mu
mod1.cr2.cov <- mod1.cr2.ci.lb < mu && mod1.cr2.ci.ub > mu
mod2.cr0.cov <- mod2.cr0.ci.lb < mu && mod2.cr0.ci.ub > mu
mod2.cr1.cov <- mod2.cr1.ci.lb < mu && mod2.cr1.ci.ub > mu
mod2.cr2.cov <- mod2.cr2.ci.lb < mu && mod2.cr2.ci.ub > mu
mod3.cr0.cov <- mod3.cr0.ci.lb < mu && mod3.cr0.ci.ub > mu
mod3.cr1.cov <- mod3.cr1.ci.lb < mu && mod3.cr1.ci.ub > mu
mod3.cr2.cov <- mod3.cr2.ci.lb < mu && mod3.cr2.ci.ub > mu
mod4.cr0.cov <- mod4.cr0.ci.lb < mu && mod4.cr0.ci.ub > mu
mod4.cr1.cov <- mod4.cr1.ci.lb < mu && mod4.cr1.ci.ub > mu
mod4.cr2.cov <- mod4.cr2.ci.lb < mu && mod4.cr2.ci.ub > mu


# sigma.s estimate
mod1.sigma.s <- NA
mod2.sigma.s <- NA
mod3.sigma.s <- mod3$sigma2[1]
mod4.sigma.s <- mod4$sigma2[1]
mod1.cr0.sigma.s <- NA
mod1.cr1.sigma.s <- NA
mod1.cr2.sigma.s <- NA
mod2.cr0.sigma.s <- NA
mod2.cr1.sigma.s <- NA
mod2.cr2.sigma.s <- NA
mod3.cr0.sigma.s <- mod3_cr0$sigma2[1]
mod3.cr1.sigma.s <- mod3_cr1$sigma2[1]
mod3.cr2.sigma.s <- mod3_cr2$sigma2[1]
mod4.cr0.sigma.s <- mod4_cr0$sigma2[1]
mod4.cr1.sigma.s <- mod4_cr1$sigma2[1]
mod4.cr2.sigma.s <- mod4_cr2$sigma2[1]


# sigma.s estimate MSE
mod1.sigma.s.mse <- NA
mod2.sigma.s.mse <- NA
mod3.sigma.s.mse <- (sigma2.s - mod3.sigma.s)^2
mod4.sigma.s.mse <- (sigma2.s - mod4.sigma.s)^2
mod1.cr0.sigma.s.mse <- NA
mod1.cr1.sigma.s.mse <- NA
mod1.cr2.sigma.s.mse <- NA
mod2.cr0.sigma.s.mse <- NA
mod2.cr1.sigma.s.mse <- NA
mod2.cr2.sigma.s.mse <- NA
mod3.cr0.sigma.s.mse <- (sigma2.s - mod3.cr0.sigma.s)^2
mod3.cr1.sigma.s.mse <- (sigma2.s - mod3.cr1.sigma.s)^2
mod3.cr2.sigma.s.mse <- (sigma2.s - mod3.cr2.sigma.s)^2
mod4.cr0.sigma.s.mse <- (sigma2.s - mod4.cr0.sigma.s)^2
mod4.cr1.sigma.s.mse <- (sigma2.s - mod4.cr1.sigma.s)^2
mod4.cr2.sigma.s.mse <- (sigma2.s - mod4.cr2.sigma.s)^2



# sigma.u estimate 
mod1.sigma.u <- mod1$sigma2 #----> we assume zero heterogeneity
mod2.sigma.u <- mod2$sigma2 
mod3.sigma.u <- mod3$sigma2[2]
mod4.sigma.u <- mod4$sigma2[2]
mod1.cr0.sigma.u <- mod1_cr0$sigma2
mod1.cr1.sigma.u <- mod1_cr1$sigma2
mod1.cr2.sigma.u <- mod1_cr2$sigma2
mod2.cr0.sigma.u <- mod2_cr0$sigma2
mod2.cr1.sigma.u <- mod2_cr1$sigma2
mod2.cr2.sigma.u <- mod2_cr2$sigma2
mod3.cr0.sigma.u <- mod3_cr0$sigma2[2]
mod3.cr1.sigma.u <- mod3_cr1$sigma2[2]
mod3.cr2.sigma.u <- mod3_cr2$sigma2[2]
mod4.cr0.sigma.u <- mod4_cr0$sigma2[2]
mod4.cr1.sigma.u <- mod4_cr1$sigma2[2]
mod4.cr2.sigma.u <- mod4_cr2$sigma2[2]


# sigma.u estimate MSE
mod1.sigma.u.mse <- (0 - mod1.sigma.u)^2 #-----> this should be zero!!
mod2.sigma.u.mse <- (sigma2.u - mod2.sigma.u)^2
mod3.sigma.u.mse <- (sigma2.u - mod3.sigma.u)^2
mod4.sigma.u.mse <- (sigma2.u - mod4.sigma.u)^2
mod1.cr0.sigma.u.mse <- (sigma2.u - mod1.cr0.sigma.u)^2
mod1.cr1.sigma.u.mse <- (sigma2.u - mod1.cr1.sigma.u)^2
mod1.cr2.sigma.u.mse <- (sigma2.u - mod1.cr2.sigma.u)^2
mod2.cr0.sigma.u.mse <- (sigma2.u - mod2.cr0.sigma.u)^2
mod2.cr1.sigma.u.mse <- (sigma2.u - mod2.cr1.sigma.u)^2
mod2.cr2.sigma.u.mse <- (sigma2.u - mod2.cr2.sigma.u)^2
mod3.cr0.sigma.u.mse <- (sigma2.u - mod3.cr0.sigma.u)^2
mod3.cr1.sigma.u.mse <- (sigma2.u - mod3.cr1.sigma.u)^2
mod3.cr2.sigma.u.mse <- (sigma2.u - mod3.cr2.sigma.u)^2
mod4.cr0.sigma.u.mse <- (sigma2.u - mod4.cr0.sigma.u)^2
mod4.cr1.sigma.u.mse <- (sigma2.u - mod4.cr1.sigma.u)^2
mod4.cr2.sigma.u.mse <- (sigma2.u - mod4.cr2.sigma.u)^2




########### Save results  ----------------------------------------------

# save results 
res <- data.frame(name = rep(name, 16),
                  scenario = rep(scen, 16), 
                  sim = rep(seed, 4),
                  CR_method = rep(c("none", "CR0", "CR1", "CR2"), each = 4),
                  model = c("FE", "RE", "ML", "ML-VCV", 
                            "FE-CR0", "RE-CR0", "ML-CR0", "ML-VCV-CR0",
                            "FE-CR1", "RE-CR1", "ML-CR1", "ML-VCV-CR1",
                            "FE-CR2", "RE-CR2", "ML-CR2", "ML-VCV-CR2"), 
                  comp.time = rep(c(mod1_time, mod2_time, mod3_time, mod4_time), 4),
                  mu_est = c(mod1.est, mod2.est, mod3.est, mod4.est,
                             mod1.cr0, mod2.cr0, mod3.cr0, mod4.cr0,
                             mod1.cr1, mod2.cr1, mod3.cr1, mod4.cr1,
                             mod1.cr2, mod2.cr2, mod3.cr2, mod4.cr2), 
                  mu_mse = c(mod1.est.mse, mod2.est.mse, mod3.est.mse, mod4.est.mse,
                             mod1.cr0.mse, mod2.cr0.mse, mod3.cr0.mse, mod4.cr0.mse,
                             mod1.cr1.mse, mod2.cr1.mse, mod3.cr1.mse, mod4.cr1.mse,
                             mod1.cr2.mse, mod2.cr2.mse, mod3.cr2.mse, mod4.cr2.mse), 
                  mu_ci_ub = c(mod1.est.ci.ub, mod2.est.ci.ub, mod3.est.ci.ub, mod4.est.ci.ub,
                               mod1.cr0.ci.ub, mod2.cr0.ci.ub, mod3.cr0.ci.ub, mod4.cr0.ci.ub,
                               mod1.cr1.ci.ub, mod2.cr1.ci.ub, mod3.cr1.ci.ub, mod4.cr1.ci.ub,
                               mod1.cr2.ci.ub, mod2.cr2.ci.ub, mod3.cr2.ci.ub, mod4.cr2.ci.ub), 
                  mu_ci_lb = c(mod1.est.ci.lb, mod2.est.ci.lb, mod3.est.ci.lb, mod4.est.ci.lb,
                               mod1.cr0.ci.lb, mod2.cr0.ci.lb, mod3.cr0.ci.lb, mod4.cr0.ci.lb,
                               mod1.cr1.ci.lb, mod2.cr1.ci.lb, mod3.cr1.ci.lb, mod4.cr1.ci.lb,
                               mod1.cr2.ci.lb, mod2.cr2.ci.lb, mod3.cr2.ci.lb, mod4.cr2.ci.lb),
                  mu_cov = c(mod1.cov, mod2.cov, mod3.cov, mod4.cov,
                             mod1.cr0.cov, mod2.cr0.cov, mod3.cr0.cov, mod4.cr0.cov,
                             mod1.cr1.cov, mod2.cr1.cov, mod3.cr1.cov, mod4.cr1.cov,
                             mod1.cr2.cov, mod2.cr2.cov, mod3.cr2.cov, mod4.cr2.cov), 
                  sigma.s_est = c(mod1.sigma.s, mod2.sigma.s, mod3.sigma.s, mod4.sigma.s,
                                  mod1.cr0.sigma.s, mod2.cr0.sigma.s, mod3.cr0.sigma.s, mod4.cr0.sigma.s,
                                  mod1.cr1.sigma.s, mod2.cr1.sigma.s, mod3.cr1.sigma.s, mod4.cr1.sigma.s,
                                  mod1.cr2.sigma.s, mod2.cr2.sigma.s, mod3.cr2.sigma.s, mod4.cr2.sigma.s), 
                  sigma.s_mse = c(mod1.sigma.s.mse, mod2.sigma.s.mse, mod3.sigma.s.mse, mod4.sigma.s.mse,
                                  mod1.cr0.sigma.s.mse, mod2.cr0.sigma.s.mse, mod3.cr0.sigma.s.mse, mod4.cr0.sigma.s.mse,
                                  mod1.cr1.sigma.s.mse, mod2.cr1.sigma.s.mse, mod3.cr1.sigma.s.mse, mod4.cr1.sigma.s.mse,
                                  mod1.cr2.sigma.s.mse, mod2.cr2.sigma.s.mse, mod3.cr2.sigma.s.mse, mod4.cr2.sigma.s.mse), 
                  sigma.u_est = c(mod1.sigma.u, mod2.sigma.u, mod3.sigma.u, mod4.sigma.u,
                                  mod1.cr0.sigma.u, mod2.cr0.sigma.u, mod3.cr0.sigma.u, mod4.cr0.sigma.u,
                                  mod1.cr1.sigma.u, mod2.cr1.sigma.u, mod3.cr1.sigma.u, mod4.cr1.sigma.u,
                                  mod1.cr2.sigma.u, mod2.cr2.sigma.u, mod3.cr2.sigma.u, mod4.cr2.sigma.u), 
                  sigma.u_mse = c(mod1.sigma.u.mse, mod2.sigma.u.mse, mod3.sigma.u.mse, mod4.sigma.u.mse,
                                  mod1.cr0.sigma.u.mse, mod2.cr0.sigma.u.mse, mod3.cr0.sigma.u.mse, mod4.cr0.sigma.u.mse,
                                  mod1.cr1.sigma.u.mse, mod2.cr1.sigma.u.mse, mod3.cr1.sigma.u.mse, mod4.cr1.sigma.u.mse,
                                  mod1.cr2.sigma.u.mse, mod2.cr2.sigma.u.mse, mod3.cr2.sigma.u.mse, mod4.sigma.u.mse), 
                  k.studies = rep(k.studies, 16),
                  sigma2.s = rep(sigma2.s, 16),
                  sigma2.u = rep(sigma2.u, 16),
                  mu = rep(mu, 16),
                  rho = rep(rho,16),
                  rho.hat = rep(rho.hat, 16)
)


# save the output according to the job array:
save(list = "res", file = paste0("results/raw/res_", job, ".RDATA"))
