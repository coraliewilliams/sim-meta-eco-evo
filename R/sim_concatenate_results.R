###############################################################################
# Combine simulation results
###############################################################################

library(dplyr)

home.wd <- "/srv/scratch/z5394590/phylo_meta_sandwich/main/study1/results/set3"

# get a list of all the individual result files
res.list <- list.files(path = home.wd, pattern = "^res_.*\\.RDATA$", full.names = TRUE)

# initialise a list to store the results
dat_list <- vector("list", length(res.list))

# loop through individual sim-by-model files
for (i in seq_along(res.list)) {
  # load the individual simulation results
  load(res.list[[i]])
  # store the result in the list
  dat_list[[i]] <- res
  # remove the current results to free memory
  rm(res)
}

# combine all data frames at once using dplyr's bind_rows
dat <- bind_rows(dat_list)

# save the single result data frame within the results folder
save(dat, file = file.path("/srv/scratch/z5394590/phylo_meta_sandwich/main/study1/results", "collated_sim_study1_results_set3.RDATA"))

# clean up
rm(dat_list)



