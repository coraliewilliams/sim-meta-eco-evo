################################################################################
# Plot results
################################################################################

load("~/Projects/sim-meta-eco-evo/output/study1/all_results_study1.RDATA")

######## Notes #############################
# Try to load Rdata of all sims from OSF repo 
# Load the Rdata file from URL in final script
########################################### 

###
# load librairies
library(pacman)
p_load(tidyverse, metafor, simstudy, ggplot2, ggpubr,
       dplyr, broom, knitr, kableExtra, here, cowplot,
       ggdist, gridExtra)

###
# reorder factor levels for model_type and model
results$model_type <- factor(results$model_type, levels=c("FE", "RE", "ML", "ML-VCV-0.2", 
                                                          "ML-VCV-0.5", "ML-VCV-0.8"))

results$model <- factor(results$model, levels=c("FE", "RE", "ML", "ML-VCV",
                                                "FE-CR0", "RE-CR0", "ML-CR0", "ML-VCV-CR0",
                                                "FE-CR1", "RE-CR1", "ML-CR1", "ML-VCV-CR1",
                                                "FE-CR2", "RE-CR2", "ML-CR2", "ML-VCV-CR2"))

results$model_name <- factor(results$model_name, 
                             levels=c("FE-none", "FE-CR0", "FE-CR1", "FE-CR2",
                                      "RE-none", "RE-CR0", "RE-CR1", "RE-CR2",
                                      "ML-none", "ML-CR0", "ML-CR1", "ML-CR2",
                                      "ML-VCV-0.2-none", "ML-VCV-0.2-CR0", "ML-VCV-0.2-CR1", "ML-VCV-0.2-CR2",
                                      "ML-VCV-0.5-none", "ML-VCV-0.5-CR0", "ML-VCV-0.5-CR1", "ML-VCV-0.5-CR2",
                                      "ML-VCV-0.8-none", "ML-VCV-0.8-CR0", "ML-VCV-0.8-CR1", "ML-VCV-0.8-CR2")) 

###
# derive confidence interval width
results <- results |>
  group_by(model, scenario) |>
  mutate(mu_ci_width = mu_ci_ub - mu_ci_lb)

# derive RMSE
results <- results |>
  group_by(model, scenario) |>
  mutate(mu_rmse = sqrt(mean((mu_est - mu)^2)))


# derive sample variance of mu estimate
sample_var <- results |> 
  group_by(model, model_type, CR_method) |>
  summarise(S2 = sum((mu_est - mean(mu_est))^2) / (length(mu_est) - 1)) |> 
  ungroup()



################################################################################
# ------------------------ PLOTS: All conditions ----------------------------
################################################################################

###
res <- results |> 
  filter(CR_method %in% c("none", "CR2")) |> 
  filter(scenario %in% c(7,8,15,16,23,24))


# derive coverage
cov_all <- results |>
  group_by(model_name, rho, sigma2.s, sigma2.u, k.studies) |>
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))


####################
# (1) Bias 
####################

mu_est_plot_all <-
  ggplot(res, aes(x=factor(model_name), y=mu_est, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Overall mean estimate",x="Model", y = "mu estimate")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.2, colour="darkgray")+
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_est_plot_all
ggsave("output/study1/All/mu_est_all.png", plot = mu_est_plot_all, width = 14, height = 8)




####################
# (2) MSE
####################


mu_mse_plot_all <-
  ggplot(res, aes(x=factor(model_name), y=mu_mse, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Mean Squared Error (MSE) of overall mean estimate",x="Model", y = "MSE")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_mse_plot_all
ggsave("output/study1/All/mu_mse_all.png", plot = mu_mse_plot_all, width = 14, height = 8)





####################
# (3) RMSE
####################

mu_rmse_plot_all <-
  ggplot(res, aes(x=factor(model_name), y=mu_rmse, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Root Mean Squared Error (RMSE) of overall mean estimate",x="Model", y = "RMSE")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_rmse_plot_all
ggsave("output/study1/All/mu_rmse_all.png", plot = mu_rmse_plot_all, width = 14, height = 8)





####################
# (4) Coverage 
####################


cov_plot_all <-
  ggplot(cov_all, aes(x=factor(model_name), y=cov_prop, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 12))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 12), 0.4)) +
  labs(title="Coverage rate of mu",x="Model", y = "Coverage rate of mu")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cov_plot_all
ggsave("output/study1/All/coverage_all.png", plot = cov_plot_all, width = 14, height = 8)


cov_plot_ML <- cov_all |> 
  filter(model_name %in% c("ML-none", 
                           "ML-VCV-0.2-none",
                           "ML-VCV-0.5-none",
                           "ML-VCV-0.8-none")) |>
  ggplot(aes(x=factor(model_name), y=cov_prop, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 6), 0.4)) +
  labs(title="Confidence interval width of mu",x="Model", y = "CI width of mu")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cov_plot_ML
ggsave("output/study1/All/coverage_all.png", plot = cov_plot_all, width = 14, height = 8)


####################
# (5) Confidence interval widths
####################

ci_plot_all <-
  ggplot(results, aes(x=factor(model_name), y=mu_ci_width, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 6), 0.4)) +
  labs(title="Confidence interval width of mu",x="Model", y = "CI width")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ci_plot_all
ggsave("output/study1/All/ci_width_all.png", plot = ci_plot_all, width = 14, height = 8)


####################
# (6) Variance estimate u (within study)
####################

sigma.u_plot_all <-
  ggplot(results, aes(x=factor(model_name), y=sigma.u_est, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 6), 0.4)) +
  labs(title="Estimate of within-study variance (s2.u)",x="Model", y = "sigma2.u")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.3, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma.u_plot_all
ggsave("output/study1/All/sigma2.u_estimate.png", plot = sigma.u_plot_all, width = 10, height = 8)


####################
# (7) Variance estimate s (between study)
####################


sigma.s_plot_all <-
  ggplot(results, aes(x=factor(model_name), y=sigma.s_est, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Estimate of between-study variance (s2.s)",x="Model", y = "sigma2.s")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.3, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma.s_plot_all
ggsave("output/study1/All/sigma2.s_estimate.png", plot = sigma.s_plot_all, width = 10, height = 8)





################################################################################
# ------------------------ PLOTS: VCV specification ----------------------------
################################################################################


### 
rv <- results |> 
  filter(CR_method %in% c("none")) |> 
  filter(scenario %in% c(7,8,15,16,23,24))


## sigma2 u (effect size level - within study)
sigma.u_plot_all <-
  ggplot(rv, aes(x=factor(model_name), y=sigma.u_est, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Estimate of within-study variance (s2.u)",x="Model", y = "sigma2.u")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.3, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  scale_y_continuous(breaks = seq(0, 1.4, by = 0.1)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma.u_plot_all
ggsave("output/study1/VCV/sigma2.u_estimate.png", plot = sigma.u_plot_all, width = 10, height = 8)


## sigma2 s (study level - between study)
sigma.s_plot_all <-
  ggplot(rv, aes(x=factor(model_name), y=sigma.s_est, color=model_name, fill=model_name)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Estimate of between-study variance (s2.s)",x="Model", y = "sigma2.s")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.3, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma.s_plot_all
ggsave("output/study1/VCV/sigma2.s_estimate.png", plot = sigma.s_plot_all, width = 10, height = 8)




################################################################################
# ------------------------ PLOTS: CRVE ----------------------------
################################################################################













################################################################################
# Save all plots
################################################################################

save(list = ls(), file = "output/study1/plots.RDATA")




##### ML results 1
results_ML1 <- results %>% 
  filter(model_full %in% c("ML-none", 
                      "ML-VCV-0.2-none",
                      "ML-VCV-0.5-none",
                      "ML-VCV-0.8-none"))
results_ML1 <-  droplevels(results_ML1)


table(results_ML1$model_full, results_ML1$rho)
table(results_ML1$model_full, results_ML1$rho.hat)

### Plot of mu estimate vs models

mu_est_plot <- ggplot(results_ML1, aes(x=factor(model_full), y=mu_est, color=model_full, fill=model_full)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="μ Estimate",x="Model", y = "μ Estimate")+
  geom_boxplot(width=0.4)+
  facet_wrap(~ factor(rho))+
  geom_hline(yintercept=0.2, colour="darkgray")+ # mu=0.2
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_est_plot
ggsave("output/study1/ML/mu_est.png", plot = mu_est_plot, width = 11, height = 7)



# Plot of mu estimate MSE vs models
mu_mse_plot <- ggplot(results_ML1, aes(x=factor(model_full), y=mu_mse, color=model_full, fill=model_full)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  labs(title="Mean Squared Error of μ estimate",x="Model", y = "MSE of μ estimate")+
  facet_wrap(~ factor(rho))+
  geom_boxplot(width=0.4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_mse_plot
ggsave("output/study1/ML/mu_mse.png", plot = mu_mse_plot, width = 11, height = 7)


# Plot of variance of mu estimate vs models
mu_var_plot <- ggplot(mu_var, aes(x=model_full, y=mu_var, color=model_full, fill=model_full)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), alpha=0.6) +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  scale_y_continuous(breaks=seq(0,0.05,0.01), limits=c(0, 0.05)) +
  labs(title="Variance of μ estimate (over 100 replications)",
       x="Model",
       y="Variance of μ estimate") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_var_plot


#### ML dataset 2
results_ML2 <- results %>% 
  filter(model_full %in% c("ML-none", 
                           "ML-VCV-0.2-none",
                           "ML-VCV-0.5-none",
                           "ML-VCV-0.8-none",
                           "ML-CR2",
                           "ML-VCV-0.2-CR2",
                           "ML-VCV-0.5-CR2",
                           "ML-VCV-0.8-CR2"))
results_ML2 <-  droplevels(results_ML2)


cov_all_sub <- cov_all |> 
  filter(model_full %in% c("ML-none", 
                           "ML-VCV-0.2-none",
                           "ML-VCV-0.5-none",
                           "ML-VCV-0.8-none"))

# Plot the coverage proportion for each model
cov_plot <- ggplot(cov_all_sub, aes(x=model_full, y=cov_prop, color=model_full, fill=model_full)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), alpha=0.6) +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5", "#F4756B"), 4), 0.4)) +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0, 1)) +
  facet_wrap(~ factor(rho))+
  labs(title="Coverage proportion of μ estimate",
       x="Model",
       y="Coverage proportion") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
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