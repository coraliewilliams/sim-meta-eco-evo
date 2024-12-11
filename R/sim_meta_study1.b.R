###############################################################################
# (STUDY 1.B) SIMULATION OF META-ANALYSIS MODELS
###############################################################################

rm(list=ls())

# Load libraries
library(pacman) # checks if package is installed, if not installs it
p_load(metafor, MASS, clubSandwich)



########## Load parameters conditions ---------------------------------------

### load job array
#tab <- read.csv("job_array_study1.csv")
tab <- read.csv("output/study1/job_array_study1.csv")

### job number from pbs script
#job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
job <- 359999 # for testing
  
### parameters for current job
name <- tab$name[tab$job_number == job] 
scen <- tab$scenario[tab$job_number == job] 
seed <- tab$sim[tab$job_number == job]
k.studies <- tab$k.studies[tab$job_number == job]  
mu <- tab$mu[tab$job_number == job]
sigma2.u <- tab$sigma2.u[tab$job_number == job]   
sigma2.s <- tab$sigma2.s[tab$job_number == job]       
rho <- tab$rho[tab$job_number == job]    


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
vi <- rbeta(k, 2, 20)*4
mean(vi)
median(vi)
max(vi)
min(vi)

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
# VCV_test <- vcalc(vi, cluster=study, obs=id, rho=rho)
# cov2cor(VCV[1:10, 1:10])


### simulate sampling error for dependent effect sizes
mi <- mvrnorm(n = 1, mu = rep(0, length(vi)), Sigma = VCV)



########### Get estimates  ------------------------------------------------------------------

yi <- mu + u.u + u.s + mi

# yi <- mu + x1 + x2 + u.u + u.s + mi

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

# get V matrix based on assumed rho of 0.2
V <- vcalc(vi, cluster=study, obs=id, rho=0.2)

ptm <- proc.time()
mod4 <- try(rma.mv(yi,
                   V=V,
                   random = list(~ 1 | study, ~ 1 | id),
                   test = "t",
                   dfs = k.studies - 1),
            silent = TRUE)
mod4_time <- (proc.time() - ptm)[3]

mod4_cr0 <- robust(mod4, cluster = study, adjust=FALSE)
mod4_cr1 <- robust(mod4, cluster = study, adjust=TRUE)
mod4_cr2 <- robust(mod4, cluster = study, clubSandwich=TRUE)




########### Run model 5: multilevel model with sampling VCV  -----------------------------------------

# get V matrix based on assumed rho of 0.5
V <- vcalc(vi, cluster=study, obs=id, rho=0.5)

ptm <- proc.time()
mod5 <- try(rma.mv(yi,
                   V=V,
                   random = list(~ 1 | study, ~ 1 | id),
                   test = "t",
                   dfs = k.studies - 1),
            silent = TRUE)
mod5_time <- (proc.time() - ptm)[3]


mod5_cr0 <- robust(mod5, cluster = study, adjust=FALSE)
mod5_cr1 <- robust(mod5, cluster = study, adjust=TRUE)
mod5_cr2 <- robust(mod5, cluster = study, clubSandwich=TRUE)



########### Run model 6: multilevel model with sampling VCV  -----------------------------------------

# get V matrix based on assumed rho of 0.8
V <- vcalc(vi, cluster=study, obs=id, rho=0.8)

ptm <- proc.time()
mod6 <- try(rma.mv(yi,
                   V=V,
                   random = list(~ 1 | study, ~ 1 | id),
                   test = "t",
                   dfs = k.studies - 1),
            silent = TRUE)
mod6_time <- (proc.time() - ptm)[3]


mod6_cr0 <- robust(mod6, cluster = study, adjust=FALSE)
mod6_cr1 <- robust(mod6, cluster = study, adjust=TRUE)
mod6_cr2 <- robust(mod6, cluster = study, clubSandwich=TRUE)






########### Extract model estimates  -------------------------------------------------------------------

# mu estimate 
mod1.est <- mod1$b
mod2.est <- mod2$b
mod3.est <- mod3$b
mod4.est <- mod4$b
mod5.est <- mod5$b
mod6.est <- mod6$b
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
mod5.cr0 <- mod5_cr0$b
mod5.cr1 <- mod5_cr1$b
mod5.cr2 <- mod5_cr2$b
mod6.cr0 <- mod6_cr0$b
mod6.cr1 <- mod6_cr1$b
mod6.cr2 <- mod6_cr2$b


# mu standard error
mod1.se <- mod1$se
mod2.se <- mod2$se
mod3.se <- mod3$se
mod4.se <- mod4$se
mod5.se <- mod5$se
mod6.se <- mod6$se
mod1.cr0.se <- mod1_cr0$se
mod1.cr1.se <- mod1_cr1$se
mod1.cr2.se <- mod1_cr2$se
mod2.cr0.se <- mod2_cr0$se
mod2.cr1.se <- mod2_cr1$se
mod2.cr2.se <- mod2_cr2$se
mod3.cr0.se <- mod3_cr0$se
mod3.cr1.se <- mod3_cr1$se
mod3.cr2.se <- mod3_cr2$se
mod4.cr0.se <- mod4_cr0$se
mod4.cr1.se <- mod4_cr1$se
mod4.cr2.se <- mod4_cr2$se
mod5.cr0.se <- mod5_cr0$se
mod5.cr1.se <- mod5_cr1$se
mod5.cr2.se <- mod5_cr2$se
mod6.cr0.se <- mod6_cr0$se
mod6.cr1.se <- mod6_cr1$se
mod6.cr2.se <- mod6_cr2$se


# mu t-statistic
mod1.zval <- mod1$zval
mod2.zval <- mod2$zval
mod3.zval <- mod3$zval
mod4.zval <- mod4$zval
mod5.zval <- mod5$zval
mod6.zval <- mod6$zval
mod1.cr0.zval <- mod1_cr0$zval
mod1.cr1.zval <- mod1_cr1$zval
mod1.cr2.zval <- mod1_cr2$zval
mod2.cr0.zval <- mod2_cr0$zval
mod2.cr1.zval <- mod2_cr1$zval
mod2.cr2.zval <- mod2_cr2$zval
mod3.cr0.zval <- mod3_cr0$zval
mod3.cr1.zval <- mod3_cr1$zval
mod3.cr2.zval <- mod3_cr2$zval
mod4.cr0.zval <- mod4_cr0$zval
mod4.cr1.zval <- mod4_cr1$zval
mod4.cr2.zval <- mod4_cr2$zval
mod5.cr0.zval <- mod5_cr0$zval
mod5.cr1.zval <- mod5_cr1$zval
mod5.cr2.zval <- mod5_cr2$zval
mod6.cr0.zval <- mod6_cr0$zval
mod6.cr1.zval <- mod6_cr1$zval
mod6.cr2.zval <- mod6_cr2$zval

# mu p-value
mod1.pval <- mod1$pval
mod2.pval <- mod2$pval
mod3.pval <- mod3$pval
mod4.pval <- mod4$pval
mod5.pval <- mod5$pval
mod6.pval <- mod6$pval
mod1.cr0.pval <- mod1_cr0$pval
mod1.cr1.pval <- mod1_cr1$pval
mod1.cr2.pval <- mod1_cr2$pval
mod2.cr0.pval <- mod2_cr0$pval
mod2.cr1.pval <- mod2_cr1$pval
mod2.cr2.pval <- mod2_cr2$pval
mod3.cr0.pval <- mod3_cr0$pval
mod3.cr1.pval <- mod3_cr1$pval
mod3.cr2.pval <- mod3_cr2$pval
mod4.cr0.pval <- mod4_cr0$pval
mod4.cr1.pval <- mod4_cr1$pval
mod4.cr2.pval <- mod4_cr2$pval
mod5.cr0.pval <- mod5_cr0$pval
mod5.cr1.pval <- mod5_cr1$pval
mod5.cr2.pval <- mod5_cr2$pval
mod6.cr0.pval <- mod6_cr0$pval
mod6.cr1.pval <- mod6_cr1$pval
mod6.cr2.pval <- mod6_cr2$pval



# mu estimate MSE
mod1.est.mse <- (mu - mod1.est)^2  
mod2.est.mse <- (mu - mod2.est)^2
mod3.est.mse <- (mu - mod3.est)^2
mod4.est.mse <- (mu - mod4.est)^2
mod5.est.mse <- (mu - mod5.est)^2
mod6.est.mse <- (mu - mod6.est)^2
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
mod5.cr0.mse <- (mu - mod5.cr0)^2
mod5.cr1.mse <- (mu - mod5.cr1)^2
mod5.cr2.mse <- (mu - mod5.cr2)^2
mod6.cr0.mse <- (mu - mod6.cr0)^2
mod6.cr1.mse <- (mu - mod6.cr1)^2
mod6.cr2.mse <- (mu - mod6.cr2)^2


# mu confidence interval
mod1.est.ci.ub <- mod1$ci.ub
mod1.est.ci.lb <- mod1$ci.lb
mod2.est.ci.ub <- mod2$ci.ub
mod2.est.ci.lb <- mod2$ci.lb
mod3.est.ci.ub <- mod3$ci.ub
mod3.est.ci.lb <- mod3$ci.lb
mod4.est.ci.ub <- mod4$ci.ub
mod4.est.ci.lb <- mod4$ci.lb
mod5.est.ci.ub <- mod5$ci.ub
mod5.est.ci.lb <- mod5$ci.lb
mod6.est.ci.ub <- mod6$ci.ub
mod6.est.ci.lb <- mod6$ci.lb
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
mod5.cr0.ci.ub <- mod5_cr0$ci.ub
mod5.cr0.ci.lb <- mod5_cr0$ci.lb
mod5.cr1.ci.ub <- mod5_cr1$ci.ub
mod5.cr1.ci.lb <- mod5_cr1$ci.lb
mod5.cr2.ci.ub <- mod5_cr2$ci.ub
mod5.cr2.ci.lb <- mod5_cr2$ci.lb
mod6.cr0.ci.ub <- mod6_cr0$ci.ub
mod6.cr0.ci.lb <- mod6_cr0$ci.lb
mod6.cr1.ci.ub <- mod6_cr1$ci.ub
mod6.cr1.ci.lb <- mod6_cr1$ci.lb
mod6.cr2.ci.ub <- mod6_cr2$ci.ub
mod6.cr2.ci.lb <- mod6_cr2$ci.lb


# mu coverage
mod1.cov <- mod1.est.ci.lb < mu && mod1.est.ci.ub > mu
mod2.cov <- mod2.est.ci.lb < mu && mod2.est.ci.ub > mu
mod3.cov <- mod3.est.ci.lb < mu && mod3.est.ci.ub > mu
mod4.cov <- mod4.est.ci.lb < mu && mod4.est.ci.ub > mu
mod5.cov <- mod5.est.ci.lb < mu && mod5.est.ci.ub > mu
mod6.cov <- mod6.est.ci.lb < mu && mod6.est.ci.ub > mu
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
mod5.cr0.cov <- mod5.cr0.ci.lb < mu && mod5.cr0.ci.ub > mu
mod5.cr1.cov <- mod5.cr1.ci.lb < mu && mod5.cr1.ci.ub > mu
mod5.cr2.cov <- mod5.cr2.ci.lb < mu && mod5.cr2.ci.ub > mu
mod6.cr0.cov <- mod6.cr0.ci.lb < mu && mod6.cr0.ci.ub > mu
mod6.cr1.cov <- mod6.cr1.ci.lb < mu && mod6.cr1.ci.ub > mu
mod6.cr2.cov <- mod6.cr2.ci.lb < mu && mod6.cr2.ci.ub > mu


# sigma.s estimate
mod1.sigma.s <- NA
mod2.sigma.s <- NA
mod3.sigma.s <- mod3$sigma2[1]
mod4.sigma.s <- mod4$sigma2[1]
mod5.sigma.s <- mod5$sigma2[1]
mod6.sigma.s <- mod6$sigma2[1]
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
mod5.cr0.sigma.s <- mod5_cr0$sigma2[1]
mod5.cr1.sigma.s <- mod5_cr1$sigma2[1]
mod5.cr2.sigma.s <- mod5_cr2$sigma2[1]
mod6.cr0.sigma.s <- mod6_cr0$sigma2[1]
mod6.cr1.sigma.s <- mod6_cr1$sigma2[1]
mod6.cr2.sigma.s <- mod6_cr2$sigma2[1]


# sigma.s estimate MSE
mod1.sigma.s.mse <- NA
mod2.sigma.s.mse <- NA
mod3.sigma.s.mse <- (sigma2.s - mod3.sigma.s)^2
mod4.sigma.s.mse <- (sigma2.s - mod4.sigma.s)^2
mod5.sigma.s.mse <- (sigma2.s - mod5.sigma.s)^2
mod6.sigma.s.mse <- (sigma2.s - mod6.sigma.s)^2
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
mod5.cr0.sigma.s.mse <- (sigma2.s - mod5.cr0.sigma.s)^2
mod5.cr1.sigma.s.mse <- (sigma2.s - mod5.cr1.sigma.s)^2
mod5.cr2.sigma.s.mse <- (sigma2.s - mod5.cr2.sigma.s)^2
mod6.cr0.sigma.s.mse <- (sigma2.s - mod6.cr0.sigma.s)^2
mod6.cr1.sigma.s.mse <- (sigma2.s - mod6.cr1.sigma.s)^2
mod6.cr2.sigma.s.mse <- (sigma2.s - mod6.cr2.sigma.s)^2



# sigma.u estimate 
mod1.sigma.u <- mod1$sigma2 #----> we assume zero heterogeneity
mod2.sigma.u <- mod2$sigma2 
mod3.sigma.u <- mod3$sigma2[2]
mod4.sigma.u <- mod4$sigma2[2]
mod5.sigma.u <- mod5$sigma2[2]
mod6.sigma.u <- mod6$sigma2[2]
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
mod5.cr0.sigma.u <- mod5_cr0$sigma2[2]
mod5.cr1.sigma.u <- mod5_cr1$sigma2[2]
mod5.cr2.sigma.u <- mod5_cr2$sigma2[2]
mod6.cr0.sigma.u <- mod6_cr0$sigma2[2]
mod6.cr1.sigma.u <- mod6_cr1$sigma2[2]
mod6.cr2.sigma.u <- mod6_cr2$sigma2[2]

# sigma.u estimate MSE
mod1.sigma.u.mse <- (0 - mod1.sigma.u)^2 #-----> this should be zero!!
mod2.sigma.u.mse <- (sigma2.u - mod2.sigma.u)^2
mod3.sigma.u.mse <- (sigma2.u - mod3.sigma.u)^2
mod4.sigma.u.mse <- (sigma2.u - mod4.sigma.u)^2
mod5.sigma.u.mse <- (sigma2.u - mod5.sigma.u)^2
mod6.sigma.u.mse <- (sigma2.u - mod6.sigma.u)^2
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
mod5.cr0.sigma.u.mse <- (sigma2.u - mod5.cr0.sigma.u)^2
mod5.cr1.sigma.u.mse <- (sigma2.u - mod5.cr1.sigma.u)^2
mod5.cr2.sigma.u.mse <- (sigma2.u - mod5.cr2.sigma.u)^2
mod6.cr0.sigma.u.mse <- (sigma2.u - mod6.cr0.sigma.u)^2
mod6.cr1.sigma.u.mse <- (sigma2.u - mod6.cr1.sigma.u)^2
mod6.cr2.sigma.u.mse <- (sigma2.u - mod6.cr2.sigma.u)^2




########### Save results  ----------------------------------------------

# save results 
res <- data.frame(name = rep(name, 24),
                  scenario = rep(scen, 24), 
                  sim = rep(seed, 6),
                  CR_method = rep(c("none", "CR0", "CR1", "CR2"), each = 6),
                  model = c("FE", "RE", "ML", "ML-VCV-0.2", "ML-VCV-0.5", "ML-VCV-0.8", 
                            "FE-CR0", "RE-CR0", "ML-CR0", "ML-VCV-02-CR0", "ML-VCV-05-CR0", "ML-VCV-08-CR0",
                            "FE-CR1", "RE-CR1", "ML-CR1",  "ML-VCV-02-CR1", "ML-VCV-05-CR1", "ML-VCV-08-CR1",
                            "FE-CR2", "RE-CR2", "ML-CR2", "ML-VCV-02-CR2", "ML-VCV-05-CR2", "ML-VCV-08-CR2"), 
                  comp.time = rep(c(mod1_time, mod2_time, mod3_time, mod4_time, mod5_time, mod6_time), 4),
                  mu_est = c(mod1.est, mod2.est, mod3.est, mod4.est, mod5.est, mod6.est,
                             mod1.cr0, mod2.cr0, mod3.cr0, mod4.cr0, mod5.cr0, mod6.cr0,
                             mod1.cr1, mod2.cr1, mod3.cr1, mod4.cr1, mod5.cr1, mod6.cr1,
                             mod1.cr2, mod2.cr2, mod3.cr2, mod4.cr2, mod5.cr2, mod6.cr2), 
                  mu_se = c(mod1.se, mod2.se, mod3.se, mod4.se, mod5.se, mod6.se,
                            mod1.cr0.se, mod2.cr0.se, mod3.cr0.se, mod4.cr0.se, mod5.cr0.se, mod6.cr0.se,
                            mod1.cr1.se, mod2.cr1.se, mod3.cr1.se, mod4.cr1.se, mod5.cr1.se, mod6.cr1.se,
                            mod1.cr2.se, mod2.cr2.se, mod3.cr2.se, mod4.cr2.se, mod5.cr2.se, mod6.cr2.se),
                  mu_tstat = c(mod1.zval, mod2.zval, mod3.zval, mod4.zval, mod5.zval, mod6.zval,
                             mod1.cr0.zval, mod2.cr0.zval, mod3.cr0.zval, mod4.cr0.zval, mod5.cr0.zval, mod6.cr0.zval,
                             mod1.cr1.zval, mod2.cr1.zval, mod3.cr1.zval, mod4.cr1.zval, mod5.cr1.zval, mod6.cr1.zval,
                             mod1.cr2.zval, mod2.cr2.zval, mod3.cr2.zval, mod4.cr2.zval, mod5.cr2.zval, mod6.cr2.zval),
                  mu_pval = c(mod1.pval, mod2.pval, mod3.pval, mod4.pval, mod5.pval, mod6.pval,
                             mod1.cr0.pval, mod2.cr0.pval, mod3.cr0.pval, mod4.cr0.pval, mod5.cr0.pval, mod6.cr0.pval,
                             mod1.cr1.pval, mod2.cr1.pval, mod3.cr1.pval, mod4.cr1.pval, mod5.cr1.pval, mod6.cr1.pval,
                             mod1.cr2.pval, mod2.cr2.pval, mod3.cr2.pval, mod4.cr2.pval, mod5.cr2.pval, mod6.cr2.pval),
                  mu_mse = c(mod1.est.mse, mod2.est.mse, mod3.est.mse, mod4.est.mse, mod5.est.mse, mod6.est.mse,
                             mod1.cr0.mse, mod2.cr0.mse, mod3.cr0.mse, mod4.cr0.mse, mod5.cr0.mse, mod6.cr0.mse,
                             mod1.cr1.mse, mod2.cr1.mse, mod3.cr1.mse, mod4.cr1.mse, mod5.cr1.mse, mod6.cr1.mse,
                             mod1.cr2.mse, mod2.cr2.mse, mod3.cr2.mse, mod4.cr2.mse, mod5.cr2.mse, mod6.cr2.mse), 
                  mu_ci_ub = c(mod1.est.ci.ub, mod2.est.ci.ub, mod3.est.ci.ub, mod4.est.ci.ub, mod5.est.ci.ub, mod6.est.ci.ub,
                               mod1.cr0.ci.ub, mod2.cr0.ci.ub, mod3.cr0.ci.ub, mod4.cr0.ci.ub, mod5.cr0.ci.ub, mod6.cr0.ci.ub,
                               mod1.cr1.ci.ub, mod2.cr1.ci.ub, mod3.cr1.ci.ub, mod4.cr1.ci.ub, mod5.cr1.ci.ub, mod6.cr1.ci.ub,
                               mod1.cr2.ci.ub, mod2.cr2.ci.ub, mod3.cr2.ci.ub, mod4.cr2.ci.ub, mod5.cr2.ci.ub, mod6.cr2.ci.ub), 
                  mu_ci_lb = c(mod1.est.ci.lb, mod2.est.ci.lb, mod3.est.ci.lb, mod4.est.ci.lb, mod5.est.ci.lb, mod6.est.ci.lb,
                               mod1.cr0.ci.lb, mod2.cr0.ci.lb, mod3.cr0.ci.lb, mod4.cr0.ci.lb, mod5.cr0.ci.lb, mod6.cr0.ci.lb,
                               mod1.cr1.ci.lb, mod2.cr1.ci.lb, mod3.cr1.ci.lb, mod4.cr1.ci.lb, mod5.cr1.ci.lb, mod6.cr1.ci.lb,
                               mod1.cr2.ci.lb, mod2.cr2.ci.lb, mod3.cr2.ci.lb, mod4.cr2.ci.lb, mod5.cr2.ci.lb, mod6.cr2.ci.lb),
                  mu_cov = c(mod1.cov, mod2.cov, mod3.cov, mod4.cov, mod5.cov, mod6.cov,
                             mod1.cr0.cov, mod2.cr0.cov, mod3.cr0.cov, mod4.cr0.cov, mod5.cr0.cov, mod6.cr0.cov,
                             mod1.cr1.cov, mod2.cr1.cov, mod3.cr1.cov, mod4.cr1.cov, mod5.cr1.cov, mod6.cr1.cov,
                             mod1.cr2.cov, mod2.cr2.cov, mod3.cr2.cov, mod4.cr2.cov, mod5.cr2.cov, mod6.cr2.cov), 
                  sigma.s_est = c(mod1.sigma.s, mod2.sigma.s, mod3.sigma.s, mod4.sigma.s, mod5.sigma.s, mod6.sigma.s,
                                  mod1.cr0.sigma.s, mod2.cr0.sigma.s, mod3.cr0.sigma.s, mod4.cr0.sigma.s, mod5.cr0.sigma.s, mod6.cr0.sigma.s,
                                  mod1.cr1.sigma.s, mod2.cr1.sigma.s, mod3.cr1.sigma.s, mod4.cr1.sigma.s, mod5.cr1.sigma.s, mod6.cr1.sigma.s,
                                  mod1.cr2.sigma.s, mod2.cr2.sigma.s, mod3.cr2.sigma.s, mod4.cr2.sigma.s, mod5.cr2.sigma.s, mod6.cr2.sigma.s), 
                  sigma.s_mse = c(mod1.sigma.s.mse, mod2.sigma.s.mse, mod3.sigma.s.mse, mod4.sigma.s.mse, mod5.sigma.s.mse, mod6.sigma.s.mse,
                                  mod1.cr0.sigma.s.mse, mod2.cr0.sigma.s.mse, mod3.cr0.sigma.s.mse, mod4.cr0.sigma.s.mse, mod5.cr0.sigma.s.mse, mod6.cr0.sigma.s.mse,
                                  mod1.cr1.sigma.s.mse, mod2.cr1.sigma.s.mse, mod3.cr1.sigma.s.mse, mod4.cr1.sigma.s.mse, mod5.cr1.sigma.s.mse, mod6.cr1.sigma.s.mse,
                                  mod1.cr2.sigma.s.mse, mod2.cr2.sigma.s.mse, mod3.cr2.sigma.s.mse, mod4.cr2.sigma.s.mse, mod5.cr2.sigma.s.mse, mod6.cr2.sigma.s.mse), 
                  sigma.u_est = c(mod1.sigma.u, mod2.sigma.u, mod3.sigma.u, mod4.sigma.u, mod5.sigma.u, mod6.sigma.u,
                                  mod1.cr0.sigma.u, mod2.cr0.sigma.u, mod3.cr0.sigma.u, mod4.cr0.sigma.u, mod5.cr0.sigma.u, mod6.cr0.sigma.u,
                                  mod1.cr1.sigma.u, mod2.cr1.sigma.u, mod3.cr1.sigma.u, mod4.cr1.sigma.u, mod5.cr1.sigma.u, mod6.cr1.sigma.u,
                                  mod1.cr2.sigma.u, mod2.cr2.sigma.u, mod3.cr2.sigma.u, mod4.cr2.sigma.u, mod5.cr2.sigma.u, mod6.cr2.sigma.u), 
                  sigma.u_mse = c(mod1.sigma.u.mse, mod2.sigma.u.mse, mod3.sigma.u.mse, mod4.sigma.u.mse, mod5.sigma.u.mse, mod6.sigma.u.mse,
                                  mod1.cr0.sigma.u.mse, mod2.cr0.sigma.u.mse, mod3.cr0.sigma.u.mse, mod4.cr0.sigma.u.mse, mod5.cr0.sigma.u.mse, mod6.cr0.sigma.u.mse,
                                  mod1.cr1.sigma.u.mse, mod2.cr1.sigma.u.mse, mod3.cr1.sigma.u.mse, mod4.cr1.sigma.u.mse, mod5.cr1.sigma.u.mse, mod6.cr1.sigma.u.mse,
                                  mod1.cr2.sigma.u.mse, mod2.cr2.sigma.u.mse, mod3.cr2.sigma.u.mse, mod4.sigma.u.mse, mod5.sigma.u.mse, mod6.sigma.u.mse), 
                  k.studies = rep(k.studies, 24),
                  sigma2.s = rep(sigma2.s, 24),
                  sigma2.u = rep(sigma2.u, 24),
                  mu = rep(mu, 24),
                  rho = rep(rho, 24)
)



# save the output according to the job array:
save(list = "res", file = paste0("results/raw/res_", job, ".RDATA"))
write.csv(res, file = paste0("output/study1/res_rho_0.8.csv"), row.names = FALSE)

