###############################################################################
# Combine simulation results
###############################################################################

home.wd <- setwd("./R/pilot/phylo/results")

# initialise the result storage
dat <- NULL

# change into appropriate result folder
setwd(paste(home.wd, "raw", sep = "/"))

# get a list of all the individual result files
res.list <- list.files()[grepl("res_", list.files(), fixed = T)]

# inner loop through individual sim-by-model files
for (job in res.list) {
  # get the job number
  job.no <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(job, ".", fixed = T), function(x){x[1]})), "_", fixed = T), function(x){x[length(x)]})))
  # load in the individual simulation results
  load(job)
  # rbind to results data frame
  dat <- rbind(dat, res)
  # remove the current results to ensure no duplicates
  rm(res)
}


# save the single result data frame within the results folder
setwd(home.wd)
save(list = "dat", file = "collated_sim_pilot_results.RDATA")


################################################################################
# Plot results
################################################################################

####### Notes
# Try to load Rdata of all sims from OSF repo 
# Load the Rdata file from URL


results <- dat
results$model <- factor(results$model, levels=c("PML", "PML-VCV", 
                                                "PML-CR0", "PML-VCV-CR0",
                                                "PML-CR1", "PML-VCV-CR1"))

# load libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggdark)
library(ggplot2)
library(ggdist)  # for half-eye plots
library(gridExtra)


# Derive coverage proportion per model
cov <- results %>%
  group_by(model, CR_method) %>%
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))


###############################################################################

# Set up plot labels 
plot_labels <- list(
  comp.time = "Runtime (seconds)",
  mu_est = "Mean μ Estimate",
  mu_bias = "Bias μ Estimate",
  mu_mse = "MSE μ Estimate",
  mu_cov = "Coverage",
  mu_ci_low = "μ CI Lower Bound",
  mu_ci_high = "μ CI Upper Bound",
  mu_ci_width = "μ CI Width",
  s2_sp_bias = "σ2 (species random effect) bias",
  s2_phylo_bias = "σ2 (phylo random effect) bias",
  s2_resid_bias = "Residual σ2 bias"
)


###############################################################################


# Plot of mu estimate vs models
mu_est_plot <- ggplot(results, aes(x=factor(model), y=mu_est, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2), 0.4)) +
  labs(title="μ Estimate",x="Model", y = "μ Estimate")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.2, colour="darkgray")+ # mu=0.2
  theme_bw()
mu_est_plot


# Plot of mu estimate MSE vs models
mu_mse_plot <- ggplot(results, aes(x=factor(model), y=mu_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2), 0.4)) +
  labs(title="Mean Squared Error of μ estimate",x="Model", y = "MSE of μ estimate")+
  geom_boxplot(width=0.4)+
  theme_bw()
mu_mse_plot



# Plot of sigma2.n estimate MSE vs models
sigma2.n_mse_plot <- ggplot(results, aes(x=factor(model), y=sigma.n_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2), 0.4)) +
  labs(title="Mean Squared Error of species level variance",x="Model", y = "MSE of variance estimate")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2.n_mse_plot


# Plot of sigma2.p estimate MSE vs models
sigma2.p_mse_plot <- ggplot(results, aes(x=factor(model), y=sigma.p_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2), 0.4)) +
  labs(title="Mean Squared Error of phylogenetic level variance ",x="Model", y = "MSE of variance estimate")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2.p_mse_plot



# Plot of sigma2.e estimate MSE vs models
sigma2.e_mse_plot <- ggplot(results, aes(x=factor(model), y=sigma.e_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2), 0.4)) +
  labs(title="Mean Squared Error of effect size level variance estimate",x="Model", y = "MSE of variance estimate of effect size level random effect")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2.e_mse_plot



# Plot of sigma.s estimate MSE vs models
sigma2.s_mse_plot <- ggplot(results, aes(x=factor(model), y=sigma.s_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2), 0.4)) +
  labs(title="Mean Squared Error of study level variance estimate",x="Model", y = "MSE of variance estimate of study level random effect")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2.s_mse_plot




# Plot the coverage proportion for each model
cov_plot <- ggplot(cov, aes(x=model, y=cov_prop, color=CR_method, fill=CR_method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), alpha=0.6) +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 2), 0.4)) +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0, 1)) +
  labs(title="Coverage proportion of μ estimate",
       x="Model",
       y="Coverage proportion") +
  theme_bw()
cov_plot




################################################################################

# Function to plot performance measures ------------

plot_results <- function(res, variable_to_plot, name="res", save=TRUE) {
  
  # Reorder models for plots
  res$model <- factor(res$model, levels=c("brms", "MCMCglmm", "pglmm", "glmmTMB"))
  
  # Create ggplot object
  gg <- ggplot(res, aes(x=model, y=get(variable_to_plot), color=model, fill=model)) +
    stat_halfeye() +  # Replace with half-eye plot from ggdist
    geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=0.6) +  # Add points
    geom_hline(aes(yintercept=0), color="white") + # Add line at zero
    geom_boxplot(width=0.1, alpha=0.4) +
    scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 3))+
    scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 3), 0.4)) + ### IMPROVE: with semi-transparent pastel
    labs(title=sprintf("Plot of %s", variable_to_plot),
         x="Package",
         y=variable_to_plot) +
    scale_shape_manual(name="Species Size", values=c(16, 17, 18, 19)) +
    facet_wrap(~ factor(species_size),
               labeller = labeller(species_size = function(x) paste("Species size:", x)))
  
  if (save) {
    filename <- sprintf("output/sim_simple/%s_%s_%diter_sim_simple.png",
                        name, variable_to_plot, iters)
    ggsave(filename, plot = gg, width = 8, height = 4)
  }
  
  return(gg)
}


# Plot each performance measure
plot_results(results, "run_time")
plot_results(results, "mu_bias")
plot_results(results, "mu_mse")
plot_results(results, "mu_ci_width")
plot_results(results, "s2_sp_bias")
plot_results(results, "s2_phylo_bias")
plot_results(results, "s2_resid_bias")





################################################################################
# STUDY 2
################################################################################

# Function to plot performance measures ------------

plot_results <- function(res, variable_to_plot, name="res", save=TRUE) {
  
  # Reorder models for plots
  res$model <- factor(res$model, levels=c("brms", "MCMCglmm", "pglmm", "glmmTMB"))
  
  # Create ggplot object
  gg <- ggplot(res, aes(x=model, y=get(variable_to_plot), color=model, fill=model)) +
    stat_halfeye() +  # Replace with half-eye plot from ggdist
    geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=0.6) +  # Add points
    geom_hline(aes(yintercept=0), color="white") + # Add line at zero
    geom_boxplot(width=0.1, alpha=0.4) +
    scale_color_manual(values=c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B")) +
    scale_fill_manual(values=alpha(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 0.4)) + # Fill with semi-transparent pastels
    labs(title=sprintf("Plot of %s", variable_to_plot),
         x="Package",
         y=variable_to_plot) +
    scale_shape_manual(name="Species Size", values=c(16, 17, 18, 19)) +
    facet_wrap(~ factor(species_size),
               labeller = labeller(species_size = function(x) paste("Species size:", x)))
  
  if (save) {
    filename <- sprintf("output/sim_simple/%s_%s_%diter_sim_simple.png",
                        name, variable_to_plot, iters)
    ggsave(filename, plot = gg, width = 8, height = 4)
  }
  
  return(gg)
}


# Plot each performance measure
plot_results(results, "run_time")
plot_results(results, "mu_bias")
plot_results(results, "mu_mse")
plot_results(results, "mu_ci_width")
plot_results(results, "s2_sp_bias")
plot_results(results, "s2_phylo_bias")
plot_results(results, "s2_resid_bias")