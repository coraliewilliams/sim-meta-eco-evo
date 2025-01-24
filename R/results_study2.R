################################################################################
# Plot results
################################################################################

load("~/PhD/2_MetaRVE/results/sims/all_results_study2.RDATA") ### change this to OSF or repo where i'll store large files
results <- combined_dat


###
# load librairies
library(pacman)
p_load(tidyverse, ggplot2, broom, here,
       cowplot, ggdist, patchwork,
       scales, latex2exp)

# reformat names (typo)
results$model <- gsub("PML-VCV-02", "PML-VCV-0.2", results$model)
results$model <- gsub("PML-VCV-05", "PML-VCV-0.5", results$model)
results$model <- gsub("PML-VCV-08", "PML-VCV-0.8", results$model)

# create new variable for plots: model_type
results$model_type <- gsub("-CR1|-CR0", "", results$model)


# reorder factor levels for model_type and model
results$model <- factor(results$model, levels=c("ML-VCV-0.5", "ML-VCV-0.5-CR0", "ML-VCV-0.5-CR1",
                                                "PML", "PML-CR0", "PML-CR1",
                                                "PML-VCV-0.2", "PML-VCV-0.2-CR0", "PML-VCV-0.2-CR1",
                                                "PML-VCV-0.5", "PML-VCV-0.5-CR0", "PML-VCV-0.5-CR1",
                                                "PML-VCV-0.8", "PML-VCV-0.8-CR0", "PML-VCV-0.8-CR1"))


results$model_type <- factor(results$model_type, levels=c("ML-VCV-0.5", 
                                                          "PML", "PML-VCV-0.2", "PML-VCV-0.5", "PML-VCV-0.8"))


results$CR_method <- factor(results$CR_method, levels=c("none", "CR0", "CR1", "CR2"))



########
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

# add labels for plots
results$rho_lab <- factor(results$rho, 
                          labels=c(
                            `0.2`=parse(text=TeX("$\\rho=0.2$")),
                            `0.5`=parse(text=TeX("$\\rho=0.5$")),
                            `0.8`=parse(text=TeX("$\\rho=0.8$"))
                          ))

results$k.studies_lab <- factor(results$k.studies, 
                                labels=c(
                                  `20`=parse(text=TeX("$N_{studies}=20$")),
                                  `50`=parse(text=TeX("$N_{studies}=50$"))
                                ))


# remove ML-VCV-0.5 model
results <- results |> 
  filter(!model_type %in% c("ML-VCV-0.5"))

# set up colors
col5mods <- c("#3CB371", "#FFD966", "#FFC845", "#FFB72A")
#col6mods <- c("#C39BD6", "#2E8B57", "#FFE066", "#FFDC4D", "#FFD733")



################################################################################
# ------------------------ a. Overall mean estimate ----------------------------
################################################################################

## data setup
# filter out CRVE methods and only keep VCV-0.5 for plots
res.a <- results |> 
  filter(CR_method %in% c("none") & 
           !model %in% c("ML-VCV-0.5")) 


# derive coverage for this subset across all varying conditions
cov.a <- res.a |>
  group_by(model_type, rho, sigma2.s, sigma2.u, sigma2.p, sigma2.n, k.studies) |> #here look at model_type as there is no CRVE methods
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))


####################
# (1) Bias 
####################

mu_est_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="", x="", y = TeX("$\\hat{\\mu}$"))+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.2, colour="darkgray", linewidth=0.5)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2/Fig1/mu_est_plot.png", plot = mu_est_plot.a, width = 10, height = 6)




####################
# (2) MSE
####################


mu_mse_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title="",x="", y = TeX("$\\hat{\\mu}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2/Fig1/mu_mse_all.png", plot = mu_mse_plot.a, width = 10, height = 6)



####################
# (3) RMSE
####################

mu_rmse_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_rmse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  # scale_y_continuous(limit=c(0, 0.2), breaks=round(seq(0, 0.2, 0.05), 3)) +
  labs(title="",x="", y = TeX("$\\hat{\\mu}$ RMSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2/Fig1/mu_rmse_all.png", plot = mu_rmse_plot.a, width = 10, height = 6)





####################
# (4) Coverage 
####################

cov_plot.a <-
  ggplot(cov.a, aes(x=factor(model_type), y=cov_prop, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2/Fig1/coverage_all.png", plot = cov_plot.a, width = 10, height = 6)



####################
# (5) Confidence interval widths
####################

ci_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2/Fig1/ci_width_all.png", plot = ci_plot.a, width = 10, height = 6)



### Combine into Figure 1
figure1 <- wrap_plots(mu_est_plot.a, mu_mse_plot.a, cov_plot.a, ci_plot.a) +
  plot_annotation(tag_levels='A')

ggsave("output/study2/SFig1/Sfigure1.png", plot = figure1, width = 10, height = 9)





################################################################################
# ------------------------ b. CRVE methods ----------------------------
################################################################################


# derive coverage for this subset across all varying conditions
cov.b <- results |>
  group_by(model_type, CR_method, rho, sigma2.s, sigma2.u, k.studies) |> 
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))


####################
# (1) Coverage 
####################

cov_plot.b <-
  cov.b |> 
  ggplot(aes(x=factor(rho), y=cov_prop, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col5mods,4, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("Coverage rate of $\\hat{\\mu}$"), x="", y = "")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  scale_y_continuous(limits=c(0.45, 0.98), breaks=seq(0, 0.95, 0.05))+
  facet_wrap(~factor(CR_method), ncol=3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/Fig2/coverage.png", plot = cov_plot.b, width = 4, height = 12)



####################
# (2) Confidence interval widths
####################

ci_plot.b <-
  results |> 
  ggplot(aes(x=factor(rho), y=mu_ci_width, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col5mods, name="model")+
  scale_fill_manual(values=alpha(col5mods, 0.4), name="model") +
  labs(title=TeX("Confidence interval width of $\\hat{\\mu}$"), x="", y ="")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(CR_method), ncol=3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/Fig2/ci_width.png", plot = ci_plot.b, width = 4, height = 12)


### Combine into Figure 4
figure4 <- cov_plot.b / ci_plot.b +
  plot_annotation(tag_levels='A')

ggsave("output/study2/Figures/figure4.png", plot = figure4, width = 11, height = 9)




### Supplementary figure: ML-VCV working models per CRVE method

####################
# (1) Coverage 
####################

cov_plot.b <-
  cov.b |> 
  ggplot(aes(x=factor(rho), y=cov_prop, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col6mods,4, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title=TeX("Coverage rate of $\\hat{\\mu}$"), x="", y = "")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~factor(CR_method), ncol=1, scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####################
# (2) Confidence interval widths
####################

ci_plot.b <-
  results |> 
  ggplot(aes(x=factor(rho), y=mu_ci_width, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col6mods, name="model")+
  scale_fill_manual(values=alpha(col6mods, 0.4), name="model") +
  labs(title=TeX("Confidence interval width of $\\hat{\\mu}$"), x="", y ="")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(CR_method), ncol=1, scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Combine into Figure 2
Sfig <- wrap_plots(cov_plot.b, ci_plot.b) +
  plot_annotation(tag_levels='A')

ggsave("output/study1/Fig2/Sfigure1_CRVE.png", plot = Sfig, width = 10, height = 12)



























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


