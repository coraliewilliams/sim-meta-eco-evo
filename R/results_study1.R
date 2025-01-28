################################################################################
# Plot results
################################################################################

load("~/PhD/2_MetaRVE/results/sims/all_results_study1.RDATA") ### change this to OSF or repo where i'll store large files


###
# load librairies
library(pacman)
p_load(tidyverse, ggplot2, broom, here,
       cowplot, ggdist, patchwork,
       scales, latex2exp)

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


results$CR_method <- factor(results$CR_method, levels=c("none", "CR0", "CR1", "CR2"))


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


# set up colors
col6mods <- c("#7297C7", "#FECA91","#A5D9A5","#D4A5E8","#C39BD6", "#A087CA")
  




################################################################################
# ------------------------ a. Overall mean estimate ----------------------------
################################################################################

## data setup
# filter out CRVE methods and only keep VCV-0.5 for plots
res.a <- results |> 
  filter(CR_method %in% c("none")) 


# derive coverage for this subset across all varying conditions
cov.a <- res.a |>
  group_by(model_type, rho, rho_lab, sigma2.s, sigma2.u, k.studies) |> #here look at model_type as there is no CRVE methods
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))


####################
# (1) Bias 
####################

mu_est_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title="", x="", y = TeX("$\\hat{\\mu}$"))+
  scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.2, colour="darkgray", linewidth=0.5)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study1/Fig1/mu_est_plot.png", plot = mu_est_plot.a, width = 10, height = 6)




####################
# (2) MSE
####################


mu_mse_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title="",x="", y = TeX("$\\hat{\\mu}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study1/Fig1/mu_mse_all.png", plot = mu_mse_plot.a, width = 10, height = 6)



####################
# (3) RMSE
####################

mu_rmse_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_rmse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
 # scale_y_continuous(limit=c(0, 0.2), breaks=round(seq(0, 0.2, 0.05), 3)) +
  labs(title="",x="", y = TeX("$\\hat{\\mu}$ RMSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study1/Fig1/mu_rmse_all.png", plot = mu_rmse_plot.a, width = 10, height = 6)





####################
# (4) Coverage 
####################

cov_plot.a <-
  ggplot(cov.a, aes(x=factor(model_type), y=cov_prop, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study1/Fig1/coverage_all.png", plot = cov_plot.a, width = 10, height = 6)



####################
# (5) Confidence interval widths
####################

ci_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=mu_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study1/Fig1/ci_width_all.png", plot = ci_plot.a, width = 10, height = 6)




### Combine into Figure 1
figure1 <- wrap_plots(mu_est_plot.a, mu_mse_plot.a, cov_plot.a, ci_plot.a) +
  plot_annotation(tag_levels='A')

ggsave("output/study1/Fig1/figure1.png", plot = figure1, width = 11, height = 8)




################################################################################
# ------------------------ b. CRVE methods ----------------------------
################################################################################

## keep only one working model VCV=0.5 for main text
## VCV-0.2 and VCV-0.8 for supplementary

# derive coverage for this subset across all varying conditions
cov.b <- results |>
  group_by(model_type, CR_method, rho, sigma2.s, sigma2.u, k.studies) |> 
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))


#### Main text figure (only VCV-0.5)

####################
# (1) Coverage 
####################

cov_plot.b <-
  cov.b |> 
  filter(!model_type %in% c("ML-VCV-0.2", "ML-VCV-0.8")) |> 
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
ggsave("output/study1/Fig2/coverage.png", plot = cov_plot.b, width = 4, height = 12)



####################
# (2) Confidence interval widths
####################

ci_plot.b <-
  results |> 
  filter(!model_type %in% c("ML-VCV-0.2", "ML-VCV-0.8")) |> 
  ggplot(aes(x=factor(rho), y=mu_ci_width, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col6mods, name="model")+
  scale_fill_manual(values=alpha(col6mods, 0.4), name="model") +
  labs(title=TeX("Confidence interval width of $\\hat{\\mu}$"), x="", y ="")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(CR_method), ncol=1, scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/Fig2/ci_width.png", plot = ci_plot.b, width = 4, height = 12)


### Combine into Figure 2
figure2 <- wrap_plots(cov_plot.b, ci_plot.b) +
  plot_annotation(tag_levels='A')

ggsave("output/study1/Fig2/figure2.png", plot = figure2, width = 9, height = 11)




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
# ------------------ c. Estimate of random effect variances----------------------
################################################################################

### the FE model is removed from total bias and MSE measures as it is not representative (I didn't extract the Tau2 measure (residual variance) from the models


### 
res.c <- results |> 
  filter(CR_method %in% c("none")) |> 
  mutate(sigma2.u_bias = sigma.u_est - sigma2.u,
         sigma2.s_bias = sigma.s_est - sigma2.s,
         sigma2.total = sigma2.u + sigma2.s,
         sigma2.total_est = replace_na(sigma.u_est, 0) + replace_na(sigma.s_est, 0),
         sigma2.total_bias = sigma2.total_est - sigma2.total,
         sigma2.total_mse = (sigma2.total_est - sigma2.total)^2)
         

# derive mean bias and MSE of variance components (if needed for effect size plots)
sigma.mean <- res.c |>
  group_by(model_type, CR_method, rho, sigma2.s, sigma2.u, k.studies) |> 
  summarise(sigma2.u_bias_mean = mean(sigma2.u_bias),
            sigma2.s_bias_mean = mean(sigma2.s_bias),
            sigma2.u_mse_mean = mean(sigma.u_mse),
            sigma2.s_mse_mean = mean(sigma.s_mse))



####################
# (1) Variance estimate u (within study)
####################

## Estimate
sigma2.u_plot_est_0.05 <-
  res.c |> 
  filter(sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.u), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.u_plot_est_0.3 <-
  res.c |> 
  filter(sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.3$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.u), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Bias 
sigma2.u_plot_bias_0.05 <-
  res.c |> 
  filter(sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.u_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.u_plot_bias_0.3 <-
  res.c |> 
  filter(sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.u_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.3$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### MSE
sigma2.u_plot_mse_0.05 <-
  res.c |> 
  filter(sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u$=0.05"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_u$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.u_plot_mse_0.3 <-
  res.c |> 
  filter(sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u$=0.3"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_u$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






####################
# (2) Variance estimates (among/between study) - only ML models 
####################


## Estimate
sigma2.s_plot_est_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_s$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.s), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.s_plot_est_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s = 0.3$"),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_s$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.s), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



### Bias
sigma2.s_plot_bias_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.s_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX(""),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_s$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sigma2.s_plot_bias_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.s_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX(""),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_s$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




## MSE
sigma2.s_plot_mse_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s$=0.05"), 
       x = "",
       y = TeX("$\\hat{\\sigma}^2_s$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.s_plot_mse_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s$=0.3"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_s$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






####################
# (3) Total variance estimate
####################

#### Subset by (0.3, 0.3) sigma2 conditions
sigma2.T_plot_bias_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE") & sigma2.u == 0.3 & sigma2.s == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.total_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[2:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[2:6], 0.4), guide="none") +
  labs(title = "", x = "", y = TeX("$\\hat{\\sigma}^2_{total}$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0, colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sigma2.T_plot_mse_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE") & sigma2.u == 0.3 & sigma2.s == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.total_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[2:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[2:6], 0.4), guide="none") +
  labs(title = "", x = "", y = TeX("$\\hat{\\sigma}^2_{total}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#### All sigma.2 conditions
sigma2.T_plot_bias_all <-
  res.c |> 
  filter(!model_type %in% c("FE")) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.total_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[2:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[2:6], 0.4), guide="none") +
  labs(title = "", x = "", y = TeX("$\\hat{\\sigma}^2_{total}$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0, colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma2.T_plot_bias_all

sigma2.T_plot_mse_all <-
  res.c |> 
  filter(!model_type %in% c("FE")) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.total_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[2:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[2:6], 0.4), guide="none") +
  labs(title = "", x = "", y = TeX("$\\hat{\\sigma}^2_{total}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma2.T_plot_mse_all



### Combine into Figure 3
figure3 <- sigma2.u_plot_est_0.3 / 
           sigma2.s_plot_est_0.3 / 
           sigma2.T_plot_mse_0.3 +
  plot_annotation(tag_levels='A')

ggsave("output/study1/Fig3/figure3.png", plot = figure3, width = 8.5, height = 12)


### Combine into supplementary figures 
### S2
supp.fig_sigma2.u.s_0.5 <- sigma2.u_plot_est_0.05 / sigma2.s_plot_est_0.05 + 
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure2_sigma2.u.s.png", 
       plot = supp.fig_sigma2.u.s_0.5, width = 9, height = 9.5)

### S3
supp.fig_sigma2.u.s_mse <- sigma2.u_plot_mse_0.05 / sigma2.u_plot_mse_0.3 +
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure3_sigma2.u_mse.png", 
       plot = supp.fig_sigma2.u.s_mse, width = 9, height = 9.5)


### S4
supp.fig_sigma2.u.s_mse <- sigma2.s_plot_mse_0.05 / sigma2.s_plot_mse_0.3 +
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure4_sigma2.s_mse.png", 
       plot = supp.fig_sigma2.u.s_mse, width = 9, height = 9.5)


### S5
supp.fig_sigma2.total <- sigma2.T_plot_bias_all / sigma2.T_plot_mse_all + 
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure5_sigma2.total.png", 
       plot = supp.fig_sigma2.total, width = 9, height = 9.5)




################################################################################
# ------------------ d. Overall mean for ML models (Extra) ----------------------
################################################################################


## data setup
# only keep VCV-0.5 and true rho=0.5 (but maybe keep all and combine more info?)
res.d <- results |> 
  filter(!model_type %in% c("FE", "RE")) 


# derive coverage for this subset across all varying conditions
cov.d <- res.d |>
  group_by(model_type, CR_method, rho, sigma2.s, sigma2.u, k.studies) |> 
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))

# set colors to have CR_method x ML models
col <- rep(col6mods[3:6],4)


####################
# (1) Coverage 
####################

# not sure yet how to present this... very small difference <5% between models

cov_plot <-
  cov.d |>
  ggplot(aes(x=factor(model_type), y=cov_prop, color=CR_method, fill=CR_method)) +
  scale_color_manual(values=col)+
  scale_fill_manual(values=alpha(col, 0.4)) +
  labs(title=TeX("by within study sampling correlation $\\rho$"), x="", y = TeX("Coverage rate of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+
  facet_grid(~ rho, scales = "free", labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






####################
# (2) Confidence interval widths
####################
## each boxplot represents 5,000 simulations for one rho condition, study size, model and CR method.


ci_plot_0.05.0.05 <-
  res.d |> 
  filter(sigma2.u == 0.05 & sigma2.s == 0.05) |>
  ggplot(aes(x=factor(CR_method), y=mu_ci_width, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col)+
  scale_fill_manual(values=alpha(col, 0.4)) +
  labs(title=TeX("$\\sigma^2_u = 0.05$ and $\\sigma^2_s = 0.05$ (by $\\rho$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_grid(factor(k.studies) ~ factor(rho), scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/xtra/ci_0.05_0.05.png", plot = ci_plot_0.05.0.05, width = 12.8, height = 8)



ci_plot_0.05.0.3 <-
  res.d |> 
  filter(sigma2.u == 0.05 & sigma2.s == 0.3) |>
  ggplot(aes(x=factor(CR_method), y=mu_ci_width, color=model_type, fill=model_type)) + 
  scale_color_manual(values=col)+
  scale_fill_manual(values=alpha(col, 0.4)) +
  labs(title=TeX("$\\sigma^2_u = 0.05$ and $\\sigma^2_s = 0.3$ (by $\\rho$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_grid(factor(k.studies) ~ factor(rho), scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/xtra/ci_plot_0.05.0.3.png", plot = ci_plot_0.05.0.3, width = 12.8, height = 8)



ci_plot_0.3.0.3 <-
  res.d |> 
  filter(sigma2.u == 0.3 & sigma2.s == 0.3) |>
  ggplot(aes(x=factor(CR_method), y=mu_ci_width, color=model_type, fill=model_type)) + 
  scale_color_manual(values=col)+
  scale_fill_manual(values=alpha(col, 0.4)) +
  labs(title=TeX("$\\sigma^2_u = 0.3$ and $\\sigma^2_s = 0.3$ (by $\\rho$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_grid(factor(k.studies) ~ factor(rho), scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/xtra/ci_plot_0.3.0.3.png", plot = ci_plot_0.3.0.3, width = 12.8, height = 8)


ci_plot_0.3.0.05 <-
  res.d |> 
  filter(sigma2.u == 0.3 & sigma2.s == 0.05) |>
  ggplot(aes(x=factor(CR_method), y=mu_ci_width, color=model_type, fill=model_type)) + 
  scale_color_manual(values=col)+
  scale_fill_manual(values=alpha(col, 0.4)) +
  labs(title=TeX("$\\sigma^2_u = 0.3$ and $\\sigma^2_s = 0.05$ (by $\\rho$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_grid(factor(k.studies) ~ factor(rho), scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/xtra/ci_plot_0.3.0.05.png", plot = ci_plot_0.3.0.05, width = 12.8, height = 8)









################################################################################
# Additional stuff 
################################################################################


## data setup
res.bb <- results |> 
  filter(!model_type %in% c("ML-VCV-0.2", "ML-VCV-0.8")) 

# derive coverage for this subset across all varying conditions
cov.bb <- res.bb |>
  group_by(model_type, CR_method, rho, sigma2.s, sigma2.u, k.studies) |> 
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))

# Coverage
cov_plot.bb <-
  ggplot(cov.bb, aes(x=factor(rho), y=cov_prop, color=model_type, fill=model_type)) + 
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title=TeX("Coverage rate of $\\hat{\\mu}$"), x="", y = "")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~factor(CR_method), ncol=1, scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# CI widths
ci_plot.bb <-
  ggplot(res.bb, aes(x=factor(rho), y=mu_ci_width, color=model_type, fill=model_type)) + 
  scale_color_manual(values=col6mods)+
  scale_fill_manual(values=alpha(col6mods, 0.4)) +
  labs(title=TeX("Confidence interval width of $\\hat{\\mu}$"), x="", y ="")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(CR_method), ncol=1, scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Combine into Figure 2.Bis
figure2_bis <- wrap_plots(cov_plot.bb, ci_plot.bb) +
  plot_annotation(tag_levels='A')

ggsave("output/study1/Fig2/figure2_bis.png", plot = figure2_bis, width = 11, height = 14)



