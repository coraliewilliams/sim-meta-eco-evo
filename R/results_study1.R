################################################################################
# Plot results
################################################################################


# expect to have 2880000 rows... (360k x 24 models) / 3 true rho values
load("~/PhD/2_MetaRVE/results/sims/all_results_study1.RDATA") ### change this to OSF or repo where i'll store large files


###
# load librairies
library(pacman)
p_load(tidyverse, ggplot2, broom, here,
       cowplot, ggdist, patchwork,
       scales, latex2exp, xtable)

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



# add labels for plots
results$rho_lab <- factor(results$rho, 
                      labels=c(
                            `0.2`=parse(text=TeX("$\\phi=0.2$")),
                            `0.5`=parse(text=TeX("$\\phi=0.5$")),
                            `0.8`=parse(text=TeX("$\\phi=0.8$"))
                            ))

results$k.studies_lab <- factor(results$k.studies, 
                            labels=c(
                      `20`=parse(text=TeX("$N_{studies}=20$")),
                      `50`=parse(text=TeX("$N_{studies}=50$"))
                    ))


# derive standard error of mu estimate from CI
results$mu_se <- (results$mu_ci_ub - results$mu_ci_lb) / (2 * qt(0.975, results$k.studies - 1))


######## For results for suppl. material
# derive CI widths with other degrees of freedom ('residual' = k - 1)
results$mu_ci_ub_df_r <- results$mu_est + (results$mu_se*qt(0.975, results$k - 1))
results$mu_ci_lb_df_r <- results$mu_est - (results$mu_se * qt(0.975, results$k - 1))
results$mu_ci_width_r <- results$mu_ci_ub_df_r - results$mu_ci_lb_df_r
results$mu_cov_r <- (results$mu_ci_lb_df_r <= results$mu) & (results$mu_ci_ub_df_r >= results$mu)

# derive CI widths with Z-test, Wald type (instead t-test)
results$mu_ci_ub_z <- results$mu_est + (results$mu_se * qnorm(0.975))
results$mu_ci_lb_z <- results$mu_est - (results$mu_se * qnorm(0.975))
results$mu_ci_width_z <- results$mu_ci_ub_z - results$mu_ci_lb_z
results$mu_cov_z <- (results$mu_ci_lb_z <= results$mu) & (results$mu_ci_ub_z >= results$mu)


## for type I error rates
# for z-test (wald CI)
results$z_stat <- (results$mu_est - results$mu) / results$mu_se
results$p_value_z <- 2 * (1 - pnorm(abs(results$z_stat)))

# for t-test (t CI)
results$t_stat <- (results$mu_est - results$mu) / results$mu_se
results$p_value_t <- 2 * (1 - pt(abs(results$t_stat), results$k.studies - 1))
results$p_value_t_df <- 2 * (1 - pt(abs(results$t_stat), results$k - 1))

# derive Type I error rates
results$mu_typeI_z <- results$p_value_z < 0.05
results$mu_typeI_t <- results$p_value_t < 0.05
results$mu_typeI_t_df <- results$p_value_t_df < 0.05


########
## Converged rate % and average comp.time of models per condition 
conv.summary <- results |> 
  group_by(model_type, k.studies) |> 
  summarise(n_models = n(), comp_time = mean(comp.time)) |> 
  mutate(perc_models = (n_models*100) / 240000, 
         comp_time = round(comp_time, 2)) |> 
  select(-n_models) |> 
  arrange(k.studies) ##to order ascendinng by k.studies
#xtable(conv.summary) ##save table for supporting information


## derive sample variance of estimates (for MCSE)
sample_var <- results |> 
  group_by(model_type, rho) |>
  summarise(mean_mu_est = mean(mu_est),
            mu_S2 = sum((mu_est - mean(mu_est))^2) / (n() - 1),
            u_S2 = sum((sigma.u_est - mean(sigma.u_est))^2) / (n() - 1),
            s_S2 = sum((sigma.s_est - mean(sigma.s_est))^2) / (n() - 1)) |> 
  ungroup()



###### set up colors
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


## derive the bias Monte Carlo SE (per model, method and condition) for overall mean
mu_mcse <- sample_var |> 
  group_by(model_type, rho) |> 
  summarise(mean_mu_est = mean_mu_est,  #get mean estimate
            mu_mcse = round(sqrt(mu_S2/n()),5)) |> 
  arrange(rho) 
#print(xtable(mu_mcse, digits=c(0,2,2,2,4)), include.rownames=FALSE) ##save table for supporting information



## derive the coverage Monte Carlo SE (per model, method and condition)
cov_mcse <- cov.a |> 
  group_by(model_type, rho, CR_method) |> 
  summarise(mean_cov = mean(cov_prop),  # Compute the mean coverage
            cov_mcse = round(sqrt((mean_cov * (1 - mean_cov)) / n()), 5)) |>  # Apply MCSE formula
  arrange(rho) 
#print(xtable(cov_mcse, digits=c(0,2,2,3,4)), include.rownames=FALSE) ##save table for supporting information



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





####################
# (1) Coverage 
####################


cov_plot.b <-
  cov.b |> 
  ggplot(aes(x=factor(rho), y=cov_prop, color=model_type, fill=model_type)) + 
  scale_color_manual(values=col6mods,4, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title=TeX("Coverage rate of $\\hat{\\mu}$"), x="", y = "")+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~factor(CR_method), ncol=1, scales="free_y")+
  geom_blank(data = subset(cov.b, CR_method %in% c("CR0", "CR1", "CR2")), aes(y = 0.92)) +
  geom_blank(data = subset(cov.b, CR_method %in% c("CR0", "CR1", "CR2")), aes(y = 0.97)) +
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
  facet_wrap(~factor(CR_method), ncol=1, scales="free_y")+
  geom_blank(data = subset(cov.b, CR_method %in% c("none", "CR0", "CR1")), aes(y = 0)) +
  geom_blank(data = subset(cov.b, CR_method %in% c("none", "CR0", "CR1")), aes(y = 3)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Combine into Figure 2
fig2 <- cov_plot.b | ci_plot.b + plot_annotation(tag_levels='A')

ggsave("output/study1/Fig2/Figure2.png", plot = fig2, width = 10, height = 12)






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
  filter(!model_type %in% c("FE", "RE") & sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.u), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.u_plot_est_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.3$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.u), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Bias 
sigma2.u_plot_bias_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.u_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.u_plot_bias_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.u_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.3$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### MSE
sigma2.u_plot_mse_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u$=0.05"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_u$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.u_plot_mse_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods[3:6], guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u$=0.3"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_u$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






####################
# (2) Variance estimates (among/between study) 
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
### S1
supp.fig_sigma2.u.s_0.5 <- sigma2.u_plot_est_0.05 / sigma2.s_plot_est_0.05 + 
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure_sigma2.u.s.png", 
       plot = supp.fig_sigma2.u.s_0.5, width = 9, height = 9.5)

### S2
supp.fig_sigma2.u.s_mse <- sigma2.u_plot_mse_0.05 / sigma2.u_plot_mse_0.3 +
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure_sigma2.u_mse.png", 
       plot = supp.fig_sigma2.u.s_mse, width = 9, height = 9.5)


### S3
supp.fig_sigma2.u.s_mse <- sigma2.s_plot_mse_0.05 / sigma2.s_plot_mse_0.3 +
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure_sigma2.s_mse.png", 
       plot = supp.fig_sigma2.u.s_mse, width = 9, height = 9.5)


### S4
supp.fig_sigma2.total <- sigma2.T_plot_bias_all / sigma2.T_plot_mse_all + 
  plot_annotation(tag_levels='A') 
ggsave("output/study1/Fig3/Sfigure_sigma2.total.png", 
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
  labs(title=TeX("by within study sampling correlation $\\phi$"), x="", y = TeX("Coverage rate of $\\hat{\\mu}$"))+
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
  labs(title=TeX("$\\sigma^2_u = 0.05$ and $\\sigma^2_s = 0.05$ (by $\\phi$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
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
  labs(title=TeX("$\\sigma^2_u = 0.05$ and $\\sigma^2_s = 0.3$ (by $\\phi$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
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
  labs(title=TeX("$\\sigma^2_u = 0.3$ and $\\sigma^2_s = 0.3$ (by $\\phi$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
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
  labs(title=TeX("$\\sigma^2_u = 0.3$ and $\\sigma^2_s = 0.05$ (by $\\phi$ and study size)"), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_grid(factor(k.studies) ~ factor(rho), scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study1/xtra/ci_plot_0.3.0.05.png", plot = ci_plot_0.3.0.05, width = 12.8, height = 8)







################################################################################
# (e) Standard error of the overall mean 
################################################################################


# plot SE per model type and CRVE method
se_plot <-
  results |> 
  filter(!model_type=="FE" & !model_type=="RE" & CR_method %in% c("none", "CR2")) |> 
  ggplot(aes(x=factor(CR_method), y=mu_se, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col6mods[3:6],4, guide="none")+
  scale_fill_manual(values=alpha(col6mods[3:6], 0.4), guide="none") +
  labs(title=TeX("Standard Error (SE) of $\\hat{\\mu}$"), x="", y = "SE")+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(rho), ncol=1, scales="free")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave("output/study1/xtra/se_plot.png", plot=se_plot, width=4, height=12)



#####################
# Derive difference in SE between models without CRVE and model with CRVE method

CR2_dat <- results |> 
  filter(CR_method=="CR2") |> 
  rename(mu_se.cr2 = mu_se)

none_dat <- results |> 
  filter(CR_method=="none")

# calculate difference in SE between none and CR2 method for each condition
se_diff <- none_dat |> 
  inner_join(CR2_dat, by=c("scenario", "sim", "model_type", "rho", "rho_lab", "sigma2.u", "sigma2.s", "k.studies")) |> 
  mutate(se_diff = mu_se - mu_se.cr2) 


# plot difference in SE per model type 
se_diff_plot <-
  se_diff |> 
  filter(!model_type=="FE" & !model_type=="RE") |> 
  ggplot(aes(x=factor(model_type), y=se_diff, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col6mods,4, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title=TeX("Difference in Standard Errors (SE) of $\\hat{\\mu}$"), x="", y = "Difference in SE (none - CR2)")+
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.9)+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed, ncol=1, scales="free")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave("output/study1/xtra/se_diff_plot.png", plot=se_diff_plot, width=4, height=12)






################################################################################
# (e) Confidence intervals of overall mean (with different settings)
################################################################################


col.3 <- c("#08519C", "#3182BD", "#6BAED6")

## Methods
# A - t-test and adjusted df
# B - z-test
# C - t-test and non-adjusted df


######## CI width

# set up CI dataset
res.ci <- results |>
  dplyr::select(model, scenario, model_type, k.studies,
                mu_ci_width, mu_ci_width_z, mu_ci_width_r,
                rho, rho_lab, sigma2.s, sigma2.u, k, CR_method) |>
  pivot_longer(
    cols = starts_with("mu_ci_width"),
    names_to = "method",
    values_to = "ci_value"
  ) |>
  mutate(method = case_when(
    method == "mu_ci_width"   ~ "A",
    method == "mu_ci_width_z" ~ "B",
    method == "mu_ci_width_r" ~ "C",
    TRUE ~ method
  ))


# plot for study=20
ci_per_study <-
  res.ci |> 
  ggplot(aes(x=factor(model_type), y=ci_value, color=method, fill=method)) + 
  scale_color_manual(values=col.3)+
  scale_fill_manual(values=alpha(col.3, 0.4)) +
  labs(title=TeX(""), x="", y = TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_grid(factor(CR_method) ~ rho_lab,
             labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/study1/inference/ci_plots.png",
       plot=ci_per_study, width=11, height=11)



######### Coverage rate

# set up coverage dataset
res.cov <- results |>
  dplyr::select(model, scenario, model_type, k.studies,
                rho, rho_lab, mu_cov, mu_cov_r, mu_cov_z, 
                sigma2.s, sigma2.u, k, CR_method) |>
  pivot_longer(
    cols = starts_with("mu_cov"),
    names_to = "method",
    values_to = "cov_value"
  ) |>
  mutate(method = case_when(
    method == "mu_cov"   ~ "A",
    method == "mu_cov_z" ~ "B",
    method == "mu_cov_r" ~ "C",
    TRUE ~ method
  ))

# get coverage rate for each model type and CRVE method
cov_prop <- res.cov |> 
  group_by(model_type, method, rho, rho_lab, k.studies,
           sigma2.u, sigma2.s, CR_method) |> 
  summarise(cov_prop = mean(cov_value, na.rm = TRUE))



# plot of coverage rate 
cov_per_method <-
  cov_prop |>
  filter(!model_type %in% c("FE", "RE")) |>
  ggplot(aes(x=factor(model_type), y=cov_prop, color=method, fill=method)) +
  scale_color_manual(values=col.3)+
  scale_fill_manual(values=alpha(col.3, 0.4)) +
  labs(title=TeX(""),
       x="",
       y = TeX("Coverage rate of $\\hat{\\mu}$"))+
  geom_boxplot()+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+
  facet_grid(cols=vars(rho_lab),
             rows=vars(CR_method),
             #scales = "free",
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/study1/inference/coverage_rate_MLonly.png",
       plot=cov_per_method, width=10, height=12)




######## Type I error rate (this is extra just for curiosity)

# set up typeI dataset
res.typeI <- results |>
  #filter(CR_method == "none") |>
  dplyr::select(model, scenario, model_type, k.studies,
                rho, rho_lab, sigma2.s, sigma2.u, k, 
                mu_typeI_t, mu_typeI_t_df, mu_typeI_z,
                CR_method) |>
  pivot_longer(
    cols = starts_with("mu_typeI"),
    names_to = "method",
    values_to = "typeI_value"
  ) |>
  mutate(method = case_when(
    method == "mu_typeI_t" ~ "A",
    method == "mu_typeI_t_df" ~ "B",
    method == "mu_typeI_z" ~ "C",
    TRUE ~ method
  ))

# get type I error rate for each model type and CRVE method
typeI_rate <- res.typeI |> 
  group_by(model_type, method, rho, rho_lab, k.studies, sigma2.u, sigma2.s, CR_method) |> 
  summarise(typeI_rate = mean(typeI_value, na.rm = TRUE)*100)



# plot type I error rate (the best would be Sattherwaite degrees of freedom)
typeI_per_method <-
  typeI_rate |>
  filter(!model_type %in% c("FE", "RE")) |>
  ggplot(aes(x=factor(model_type), y=typeI_rate, color=method, fill=method)) +
  scale_color_manual(values=col.3)+
  scale_fill_manual(values=alpha(col.3, 0.4)) +
  labs(title=TeX("by within study sampling correlation $\\phi$ and CRVE method"),
       x="",
       y = TeX("Type I error rate (%) of $\\hat{\\mu}$"))+
  geom_boxplot()+
  geom_hline(yintercept=5, colour="darkgray", linewidth=0.6)+
  facet_grid(cols=vars(rho_lab),
             rows=vars(CR_method),
             #scales = "free",
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/study1/inference/typeI_error.png",
       plot=typeI_per_method, width=10, height=12)

