###############################################################################
# SIMULATION PARAMETERS AND JOB ARRAY SET UP
###############################################################################


#############
## STUDY 1 ## --------------------------------------------------------------------------
#############

name.1 <- "study1"

# model parameters
k.studies <- c(20, 50)        # number of studies
mu <- 0.2                     # overall mean effect size (assume a Z value)
sigma2.s <- c(0.05, 0.3)      # study level variance component: two conditions=0.05, 0.30
sigma2.u <- c(0.05, 0.3)      # estimate level variance component: two conditions=0.05, 0.30
rho <- c(0.2, 0.5, 0.8)       # true correlation between effect sizes within a study
rho.hat <- c(0.2, 0.5, 0.8)   # true correlation between effect sizes within a study

# number of replications per scenario
repl <- 5000
sim <- rep(1:repl)

# make table of all scenarios
scen.tab1 <- expand.grid(sim = sim, name = name.1, mu = mu, k.studies = k.studies,
                        sigma2.s = sigma2.s, sigma2.u = sigma2.u, rho = rho, rho.hat = rho.hat)

# add columns for to store results and job number
scen.tab1$save_location <- rep("/srv/scratch/z5394590/phylo_meta_sandwich/", each=nrow(scen.tab1))
scen.tab1$job_number <- c(1:nrow(scen.tab1))
conds.1 <- length(k.studies) * length(sigma2.s) * length(sigma2.u) * length(rho) * length(rho.hat) 
scen.tab1$scenario <- rep(1:conds.1, each = repl)

# save as csv file
write.csv(scen.tab1, "output/study1/job_array_study1.csv", row.names = FALSE)



#############
## STUDY 2 ## -----------------------------------------------------------------------
#############


name.2 <- "study2"

# model parameters
k.species <- c(40, 100)   # number of species (must be <= k.studies)
sigma2.n <- c(0.05, 0.3)  # variance of the species level random effect
sigma2.p <- c(0.05, 0.3)  # variance of phylogenetic random effect

# make table of all scenarios
scen.tab2 <- expand.grid(sim = sim, name = name.2, k.studies = k.studies, 
                         k.species = k.species, mu = mu, sigma2.n = sigma2.n,
                         sigma2.p = sigma2.p, sigma2.s = sigma2.s, sigma2.u = sigma2.u,
                         rho = rho)

# only keep rows with the following c(k.studies, k.species) combination: c(20, 40), c(50, 100)
library(dplyr)
scen.tab2 <- scen.tab2 %>% 
  filter((k.studies == 20 & k.species == 40) | (k.studies == 50 & k.species == 100))


# add columns for to store results and job number
scen.tab2$save_location <- rep("/srv/scratch/z5394590/phylo_meta_sandwich/", each=nrow(scen.tab2))
scen.tab2$job_number <- c(1:nrow(scen.tab2))
conds.2 <- 2 * length(sigma2.s) * length(sigma2.u) * length(sigma2.n) * length(sigma2.p)  * length(rho)
scen.tab2$scenario <- rep(1:conds.2, each = repl)

# save as csv file
write.csv(scen.tab2, "output/study2/job_array_study2.csv", row.names = FALSE)






###############
## STUDY 2.sub ## -----------------------------------------------------------------------
###############


name.2.sub <- "study2.sub"

# model parameters
k.studies <- c(20)        # number of studies
k.species <- c(40)            # number of species (must be <= k.studies)
sigma2.s <- c(0.3)      # study level variance component: two conditions=0.05, 0.30
sigma2.u <- c(0.3)      # estimate level variance component: two conditions=0.05, 0.30
sigma2.n <- c(0.3)  # variance of the species level random effect
sigma2.p <- c(0.3)  # variance of phylogenetic random effect
rho <- c(0.2, 0.5, 0.8)       # true correlation between effect sizes within a study

repl <- 1000
sim <- rep(1:repl)

# make table of all scenarios
scen.tab2.sub <- expand.grid(sim = sim, name = name.2.sub, 
                             k.studies = k.studies, k.species = k.species,
                             sigma2.n = sigma2.n, sigma2.p = sigma2.p,
                             sigma2.s = sigma2.s, sigma2.u = sigma2.u, rho = rho)

# only keep rows with the following c(k.studies, k.species) combination: c(20, 40), c(50, 100)
# scen.tab2 <- scen.tab2 %>% 
#   filter((k.studies == 20 & k.species == 40) | (k.studies == 50 & k.species == 100))


# add columns for to store results and job number
scen.tab2.sub$save_location <- rep("/srv/scratch/z5394590/phylo_meta_sandwich/", each=nrow(scen.tab2.sub))
scen.tab2.sub$job_number <- c(1:nrow(scen.tab2.sub))
conds.2 <- length(sigma2.s) * length(sigma2.u) * length(sigma2.n) * length(sigma2.p)  * length(rho)
scen.tab2.sub$scenario <- rep(1:conds.2, each = repl)

# save as csv file
write.csv(scen.tab2.sub, "output/study2/job_array_study2.sub.csv", row.names = FALSE)

