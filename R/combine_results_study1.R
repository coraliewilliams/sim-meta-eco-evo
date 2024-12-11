################################################################################
# Combine all results
################################################################################

library(dplyr)

# load result files
load("output/study1/collated_sim_study1_results_set1.RDATA")
dat1 <- dat
dat1$rho.hat <- rep(0.2, length(dat1$rho))
load("output/study1/collated_sim_study1_results_set2.RDATA")
dat2 <- dat
load("output/study1/collated_sim_study1_results_set3.RDATA")
dat3 <- dat

# combine all results
combined_dat <- rbind(dat1, dat2, dat3)


# create new variables for model type based on rhot.hat and CR method
combined_dat <- combined_dat |>
  mutate(model_type = case_when(
    str_detect(model, "FE") ~ "FE",
    str_detect(model, "RE") ~ "RE",
    str_detect(model, "ML-VCV") & rho.hat==0.2 ~ "ML-VCV-0.2",
    str_detect(model, "ML-VCV") & rho.hat==0.5 ~ "ML-VCV-0.5",
    str_detect(model, "ML-VCV") & rho.hat==0.8 ~ "ML-VCV-0.8",
    str_detect(model, "ML") ~ "ML"),
    model_name = paste(model_type, CR_method, sep="-"))

# select main results including m.vcv.0.2 models
m <- combined_dat |> 
  filter((scenario %in% c(1:24)))

# subset vcv.0.5 and vcv.0.8 models to match the main scenarios
m.vcv.05 <- combined_dat |> 
  filter((model_type %in% c("ML-VCV-0.5") & scenario %in% c(25:48))) |> 
  mutate(scenario = scenario - 24)
m.vcv.08 <- combined_dat |> 
  filter((model_type %in% c("ML-VCV-0.8") & scenario %in% c(49:72))) |> 
  mutate(scenario = scenario - 48)

# combine all 
results <- rbind(m, m.vcv.05, m.vcv.08)


# replace rho.hat for models other than ML with NA
results <- results |>
  mutate(rho.hat = ifelse(model_type %in% c("FE", "RE", "ML"), NA, rho.hat))

## check: 
## table(results$rho.hat, results$rho)
## table(results$model_type, results$rho.hat, useNA="ifany")
## table(results$model_type, results$rho, useNA="ifany")

# save 
save(results, file="output/study1/all_results_study1.RDATA")
rm(dat, dat1, dat2, dat3, combined_dat)

