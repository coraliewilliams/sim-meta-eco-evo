################################################################################
# Combine all results
################################################################################

# Load libraries
library(dplyr)

# Load the datasets
load("output/study2/collated_sim_study2_results_set1.RDATA")
dat1 <- dat

load("output/study2/collated_sim_study2_results_set2.RDATA")
dat2 <- dat

load("output/study2/collated_sim_study2_results_set3.RDATA")
dat3 <- dat

# Combine the datasets
combined_dat <- rbind(dat1, dat2, dat3)

# save 
save(combined_dat, file="output/study2/all_results_study2.RDATA")

# Optionally, clean up the environment by removing the individual datasets
rm(dat, dat1, dat2, dat3)


## check: 
## table(results$rho.hat, results$rho)
## table(results$model_type, results$rho.hat, useNA="ifany")
## table(results$model_type, results$rho, useNA="ifany")


# Table of convergence rates ----
(table(results$model)/480000)*100
