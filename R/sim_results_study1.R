################################################################################
# Combine all results
################################################################################

# Load libraries
library(dplyr)

# Load scenarios table
scen <- read.csv("output/study1/job_array_study1.csv")
scen_rho.hat <- scen |> select(sim, scenario, rho.hat)

# Load the datasets
load("output/study1/collated_sim_study1_results_set1.RDATA")
dat1 <- dat
dat1 <- dat1 |> merge(scen_rho.hat)

load("output/study1/collated_sim_study1_results_set2.RDATA")
dat2 <- dat

load("output/study1/collated_sim_study1_results_set3.RDATA")
dat3 <- dat

# Combine the datasets
combined_dat <- rbind(dat1, dat2, dat3)

# save 
save(combined_dat, file="output/study1/all_results_study1.RDATA")

# Optionally, clean up the environment by removing the individual datasets
rm(dat, dat1, dat2, dat3)



################################################################################
# Plot results
################################################################################

####### Notes
# Try to load Rdata of all sims from OSF repo 
# Load the Rdata file from URL


results <- combined_dat

results$model <- factor(results$model, levels=c("FE", "RE", "ML", "ML-VCV",
                                                "FE-CR0", "RE-CR0", "ML-CR0", "ML-VCV-CR0",
                                                "FE-CR1", "RE-CR1", "ML-CR1", "ML-VCV-CR1",
                                                "FE-CR2", "RE-CR2", "ML-CR2", "ML-VCV-CR2"))

###############################################################################
# Filter by study size

results <- results %>%
  filter(k.studies == 20 & sigma2.s == 0.05 & sigma2.u == 0.05)

# Filter out all CR0 and CR2 models
results <- results %>%
  filter(CR_method != "CR0" & CR_method != "CR2")


###############################################################################

# load libraries
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggdark)
library(ggplot2)
library(ggdist)  # for half-eye plots
library(gridExtra)


# Derive coverage proportion per model ###### ----> separate by rho and rho.hat
cov <- results %>%
  group_by(model, CR_method) %>%
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))

# Derive variance of mu estimate per model
mu_var <- results %>%
  group_by(model, CR_method) %>%
  summarise(mu_var = var(mu_est, na.rm = TRUE))


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
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="μ Estimate",x="Model", y = "μ Estimate")+
  geom_boxplot(width=0.4)+
  facet_wrap(~ factor(rho))+
  geom_hline(yintercept=0.2, colour="darkgray")+ # mu=0.2
  theme_bw()
mu_est_plot


# Plot of mu estimate MSE vs models
mu_mse_plot <- ggplot(results, aes(x=factor(model), y=mu_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="Mean Squared Error of μ estimate",x="Model", y = "MSE of μ estimate")+
  facet_wrap(~ factor(rho))+
  geom_boxplot(width=0.4)+
  theme_bw()
mu_mse_plot

# Plot of variance of mu estimate vs models
mu_var_plot <- ggplot(mu_var, aes(x=model, y=mu_var, color=CR_method, fill=CR_method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), alpha=0.6) +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  scale_y_continuous(breaks=seq(0,0.05,0.01), limits=c(0, 0.05)) +
  labs(title="Variance of μ estimate (over 100 replications)",
       x="Model",
       y="Variance of μ estimate") +
  theme_bw()
mu_var_plot


# Plot the coverage proportion for each model
cov_plot <- ggplot(cov, aes(x=model, y=cov_prop, color=CR_method, fill=CR_method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), alpha=0.6) +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0, 1)) +
  labs(title="Coverage proportion of μ estimate (over 100 replications)",
       x="Model",
       y="Coverage proportion") +
  theme_bw()
cov_plot




# Plot of sigma2.e estimate MSE vs models
sigma2.e_mse_plot <- ggplot(results, aes(x=factor(model), y=sigma.u_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="Mean Squared Error of effect size level variance estimate",x="Model", y = "MSE of variance estimate of effect size level random effect")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2.e_mse_plot



# Plot of sigma.s estimate MSE vs models

results_ml <- results %>% 
  filter(model %in% c("ML", "ML-VCV",
                      "ML-CR0", "ML-VCV-CR0",
                      "ML-CR1", "ML-VCV-CR1",
                      "ML-CR2", "ML-VCV-CR2"))

sigma2.s_mse_plot <- ggplot(results_ml, aes(x=factor(model), y=sigma.s_mse, color=CR_method, fill=CR_method)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="Mean Squared Error of study level variance estimate",x="Model", y = "MSE of variance estimate of study level random effect")+
  geom_boxplot(width=0.4)+
  theme_bw()
sigma2.s_mse_plot



# Plot the coverage proportion for each model
# SPLIT BY STUDY NUMBER:
# CR1 and CR2 will be better for small sample sizes
cov_plot <- ggplot(cov, aes(x=model, y=cov_prop, color=CR_method, fill=CR_method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), alpha=0.6) +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
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
