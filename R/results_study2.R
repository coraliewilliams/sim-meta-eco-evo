################################################################################
# Plot results
################################################################################

load("~/Projects/sim-meta-eco-evo/output/study2/all_results_study2.RDATA")
results <- combined_dat

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



# reorder factor levels for model_type and model
results$model <- factor(results$model, levels=c("ML-VCV-0.5", "ML-VCV-0.5-CR0", "ML-VCV-0.5-CR1",
                                                "PML", "PML-CR0", "PML-CR1",
                                                "PML-VCV-02", "PML-VCV-0.2-CR0", "PML-VCV-0.2-CR1",
                                                "PML-VCV-05", "PML-VCV-0.5-CR0", "PML-VCV-0.5-CR1",
                                                "PML-VCV-08", "PML-VCV-0.8-CR0", "PML-VCV-0.8-CR1"))

# check
# table(results$model)


########
# derive confidence interval width
results <- results |>
  group_by(model, scenario) |>
  mutate(mu_ci_width = mu_ci_ub - mu_ci_lb)

# derive RMSE
results <- results |>
  group_by(model, scenario) |>
  mutate(mu_rmse = sqrt(mean((mu_est - mu)^2)))


################################################################################
# ------------------------ PLOTS: All conditions ----------------------------
################################################################################

###
res <- results |> 
  filter(CR_method %in% c("none", "CR1")) 


# derive coverage
cov_all <- results |>
  group_by(model, rho, sigma2.n, sigma2.p, sigma2.s, sigma2.u) |>
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))



####################
# (1) Bias 
####################

mu_est_plot_all <-
  ggplot(res, aes(x=factor(model), y=mu_est, color=model, fill=model)) + 
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
ggsave("output/study2/All/mu_est_all.png", plot = mu_est_plot_all, width = 14, height = 8)




####################
# (2) MSE
####################


mu_mse_plot_all <-
  ggplot(res, aes(x=factor(model), y=mu_mse, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Mean Squared Error (MSE) of overall mean estimate",x="Model", y = "MSE")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_mse_plot_all
ggsave("output/study2/All/mu_mse_all.png", plot = mu_mse_plot_all, width = 14, height = 8)





####################
# (3) RMSE
####################

mu_rmse_plot_all <-
  ggplot(res, aes(x=factor(model), y=mu_rmse, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91"), 6), 0.4)) +
  labs(title="Root Mean Squared Error (RMSE) of overall mean estimate",x="Model", y = "RMSE")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mu_rmse_plot_all
ggsave("output/study2/All/mu_rmse_all.png", plot = mu_rmse_plot_all, width = 14, height = 8)





####################
# (4) Coverage 
####################


cov_plot_all <-
  ggplot(cov_all, aes(x=factor(model), y=cov_prop, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 12))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 12), 0.4)) +
  labs(title="Coverage rate of mu",x="Model", y = "Coverage rate of mu")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cov_plot_all
ggsave("output/study2/All/coverage_all.png", plot = cov_plot_all, width = 14, height = 8)




####################
# (5) Confidence interval widths
####################

ci_plot_all <-
  ggplot(results, aes(x=factor(model), y=mu_ci_width, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 6), 0.4)) +
  labs(title="Confidence interval width of mu",x="Model", y = "CI width")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ci_plot_all
ggsave("output/study2/All/ci_width_all.png", plot = ci_plot_all, width = 14, height = 8)



### Filter by 
res2 <- results |> 
  filter(sigma2.s == 0.3, sigma2.u == 0.3) 





####################
# (6) Variance estimate u (within study)
####################

sigma.u_plot_all <-
  ggplot(res2, aes(x=factor(model), y=sigma.u_est, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 6), 0.4)) +
  labs(title="Estimate of within-study variance (s2.u)",x="Model", y = "sigma2.u")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.3, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma.u_plot_all
ggsave("output/study2/All/sigma2.u_estimate.png", plot = sigma.u_plot_all, width = 10, height = 8)


####################
# (7) Variance estimate s (between study)
####################


sigma.s_plot_all <-
  ggplot(res2, aes(x=factor(model), y=sigma.s_est, color=model, fill=model)) + 
  stat_halfeye() +
  scale_color_manual(values=rep(c("#7297C7", "#FECA91", "#A5D9A5"), 6))+
  scale_fill_manual(values=alpha(rep(c("#7297C7", "#FECA91", "#A5D9A5"), 6), 0.4)) +
  labs(title="Estimate of between-study variance (s2.s)",x="Model", y = "sigma2.s")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.3, colour="darkgray")+ 
  facet_wrap(~factor(rho))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma.s_plot_all
ggsave("output/study2/All/sigma2.s_estimate.png", plot = sigma.s_plot_all, width = 10, height = 8)


