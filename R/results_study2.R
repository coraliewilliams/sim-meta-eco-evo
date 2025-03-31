################################################################################
# Plot results
################################################################################

load("~/PhD/2_MetaRVE/results/sims/all_results_study2.RDATA") ### change this to OSF or repo where i'll store large files
results <- combined_dat

# expect to have 7,680,000 rows... (480k x 15 models)


###
# load librairies
library(pacman)
p_load(tidyverse, ggplot2, broom, here,
       cowplot, ggdist, patchwork,
       scales, latex2exp, xtable)

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
                            `0.2`=parse(text=TeX("$\\phi=0.2$")),
                            `0.5`=parse(text=TeX("$\\phi=0.5$")),
                            `0.8`=parse(text=TeX("$\\phi=0.8$"))
                          ))

results$k.studies_lab <- factor(results$k.studies, 
                                labels=c(
                                  `20`=parse(text=TeX("$N_{studies}=20$")),
                                  `50`=parse(text=TeX("$N_{studies}=50$"))
                                ))


########
# derive standard error of mu estimate from CI
results$mu_se <- (results$mu_ci_ub - results$mu_ci_lb) / (2 * qt(0.975, results$k.studies - 1))


########
# derive confidence interval with approx. df Satterwaithe method (3/((1/k.studies) + (1/k.species) + (1/k.species)))
# alternative: m-l-1 (m refers to the number of upper-level units and l refers to the number of contextual variables)
results <- results |>
  mutate(df.hm = (3/(1/k.studies + 1/k.species + 1/k.species))) 

# obtain confidence interval and coverage rates
results$mu_ci_ub_df.hm <- results$mu_est + (results$mu_se * qt(0.975, results$df.hm - 1))
results$mu_ci_lb_df.hm <- results$mu_est - (results$mu_se * qt(0.975, results$df.hm - 1))
results$mu_ci_width_df.hm <- results$mu_ci_ub_df.hm - results$mu_ci_lb_df.hm
results$mu_cov_df.hm <- (results$mu_ci_lb_df.hm <= results$mu) & (results$mu_ci_ub_df.hm >= results$mu)




# remove ML-VCV-0.5 model --- expect now to have 5,760,000 rows... (480k x 12 models)
results <- results |> 
  filter(!model_type %in% c("ML-VCV-0.5"))



########
## Converged rate % and average comp.time of models per condition 
conv.summary <- results |> 
  group_by(model_type, k.studies) |> 
  summarise(n_models = n(), comp_time = mean(comp.time)) |> 
  mutate(perc_models = (n_models*100) / 720000, 
         comp_time = round(comp_time, 2)) |> 
  select(-n_models) |> 
  arrange(k.studies) ##to order ascendinng by k.studies
#print(xtable(conv.summary, digits=c(0,0,2,2,2)), include.rownames=FALSE) ##save table for supporting information


## derive sample variance of estimates (for MCSE)
sample_var <- results |> 
  group_by(model_type, rho) |>
  summarise(mu_mean = mean(mu_est),
            u_mean = mean(sigma.u_est),
            s_mean = mean(sigma.s_est),
            n_mean = mean(sigma.n_est),
            p_mean = mean(sigma.p_est),
            mu_S2 = sum((mu_est - mean(mu_est))^2) / (n() - 1),
            u_S2 = sum((sigma.u_est - mean(sigma.u_est))^2) / (n() - 1),
            s_S2 = sum((sigma.s_est - mean(sigma.s_est))^2) / (n() - 1),
            n_S2 = sum((sigma.n_est - mean(sigma.n_est))^2) / (n() - 1),
            p_S2 = sum((sigma.p_est - mean(sigma.p_est))^2) / (n() - 1)) |> 
  ungroup()


## derive the bias Monte Carlo SE (per model, method and condition) for overall mean
mu_mcse <- sample_var |> 
  group_by(model_type, rho, mu_mean) |> 
  summarise(mu_mcse = round(sqrt(mu_S2/n()),5)) |> 
  arrange(rho)
#print(xtable(mu_mcse, digits=c(0,2,2,4,3)), include.rownames=FALSE) ##save table for supporting information





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
  group_by(model_type, rho,  rho_lab, sigma2.s, sigma2.u, sigma2.p, sigma2.n, k.studies) |> #here look at model_type as there is no CRVE methods
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
  facet_wrap(k.studies~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


cov_plot.hm <-
  ggplot(cov.hm, aes(x=factor(model_type), y=cov_prop, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




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



### Combine into Figure 1
figure1 <- wrap_plots(mu_est_plot.a, mu_mse_plot.a, cov_plot.a, ci_plot.a) +
  plot_annotation(tag_levels='A')

ggsave("output/study2/SFig/Sfigure5_PML_means.png", plot = figure1, width = 12, height = 11)





################################################################################
# ------------------------ b. CRVE methods ----------------------------
################################################################################


# derive coverage for this subset across all varying conditions
cov.b <- results |>
  group_by(model_type, CR_method, rho, sigma2.s, sigma2.u, k.studies) |> 
  summarise(cov_prop = mean(mu_cov, na.rm = TRUE))

## derive the coverage Monte Carlo SE (per model, method and condition)
cov_mcse <- cov.b |> 
  filter(CR_method=="none") |> 
  group_by(model_type, rho) |> 
  summarise(mean_cov = mean(cov_prop),  # Compute the mean coverage
            cov_mcse = round(sqrt((mean_cov * (1 - mean_cov)) / n()), 5)) |>  # Apply MCSE formula
  arrange(rho) 

#print(xtable(cov_mcse, digits=c(0,2,2,3,4)), include.rownames=FALSE) ##save table for supporting information



####################
# (1) Coverage 
####################

cov_plot.b <-
  cov.b |> 
  ggplot(aes(x=factor(rho), y=cov_prop, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col5mods,4, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="", x=TeX("$\\phi$"), y =TeX("Coverage rate of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  scale_y_continuous(limits=c(0.38, 0.98), breaks=seq(0.05, 0.95, 0.1))+
  facet_wrap(~factor(CR_method), ncol=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



####################
# (2) Confidence interval widths
####################

ci_plot.b <-
  results |> 
  ggplot(aes(x=factor(rho), y=mu_ci_width, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col5mods, name="model")+
  scale_fill_manual(values=alpha(col5mods, 0.4), name="model") +
  labs(title="", x=TeX("$\\phi$"), y=TeX("Confidence interval width of $\\hat{\\mu}$"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~factor(CR_method), ncol=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Combine into Figure 4
figure4 <- cov_plot.b + ci_plot.b +
  plot_annotation(tag_levels='A')

ggsave("output/study2/Figures/figure4.png", plot = figure4, width = 10, height = 12)





################################################################################
# ------------------ c. Estimate of random effect variances----------------------
################################################################################

### the FE model is removed from total bias and MSE measures as it is not representative (I didn't extract the Tau2 measure (residual variance) from the models


### 
res.c <- results |> 
  filter(CR_method %in% c("none")) |> 
  mutate(sigma2.u_bias = sigma.u_est - sigma2.u,
         sigma2.s_bias = sigma.s_est - sigma2.s,
         sigma2.n_bias = sigma.n_est - sigma2.n,
         sigma2.p_bias = sigma.p_est - sigma2.p,
         sigma2.total = sigma2.u + sigma2.s + sigma2.n + sigma2.p,
         sigma2.total_est = sigma.u_est + sigma.s_est + sigma.n_est + sigma.p_est,
         sigma2.total_bias = sigma2.total_est - sigma2.total,
         sigma2.total_mse = (sigma2.total_est - sigma2.total)^2)


# derive mean bias and MSE of variance components (if needed for effect size plots)
sigma.mean <- res.c |>
  group_by(model_type, CR_method, rho, sigma2.s, sigma2.u, k.studies) |> 
  summarise(sigma2.u_bias_mean = mean(sigma2.u_bias),
            sigma2.s_bias_mean = mean(sigma2.s_bias),
            sigma2.u_mse_mean = mean(sigma.u_mse),
            sigma2.s_mse_mean = mean(sigma.s_mse))


## derive Monte Carlo SE (per model, method and condition) for all variance components
s2_mcse <- sample_var |> 
  group_by(model_type, rho) |> 
  summarise(s2.u_mcse = round(sqrt(u_S2/n()),5),
            s2.s_mcse = round(sqrt(s_S2/n()),5),
            s2.n_mcse = round(sqrt(n_S2/n()),5),
            s2.p_mcse = round(sqrt(p_S2/n()),5)) |> 
  arrange(rho)
#print(xtable(s2_mcse, digits=c(0,2,2,3,3,3,3)), include.rownames=FALSE) ##save table for supporting information




####################
# (1) Variance estimates 
####################

## Estimate
sigma2.u_plot_est_0.05 <-
  res.c |> 
  filter(sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.u), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sigma2.u_plot_est_0.3 <-
  res.c |> 
  filter(sigma2.p == 0.3 & sigma2.n == 0.3 & sigma2.s == 0.3 & sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u = 0.3$"), x = "", y = TeX("$\\hat{\\sigma}^2_u$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.u), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## Among study (s)
sigma2.s_plot_est_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_s$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.s), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.s_plot_est_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.p == 0.3 & 
           sigma2.n == 0.3 & sigma2.s == 0.3 & sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s = 0.3$"),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_s$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.s), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




## Among species non-phylogeny (n)
sigma2.n_plot_est_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.n == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.n_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_n = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_n$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.n), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.n_plot_est_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE")  & sigma2.p == 0.3 & 
           sigma2.n == 0.3 & sigma2.s == 0.3 & sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_n = 0.3$"),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_n$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.n), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## Among species phylogeny (p)
sigma2.p_plot_est_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.p == 0.05 & sigma2.n == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.n_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_p = 0.05$"), x = "", y = TeX("$\\hat{\\sigma}^2_p$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.p), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.p_plot_est_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.p == 0.3 & 
           sigma2.n == 0.3 & sigma2.s == 0.3 & sigma2.u == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_p = 0.3$"),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_p$"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = sigma2.p), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



### Combine into Figure 5
figure5 <- sigma2.u_plot_est_0.3 / 
           sigma2.s_plot_est_0.3 / 
           sigma2.n_plot_est_0.3 / 
          sigma2.p_plot_est_0.3 +
  plot_annotation(tag_levels='A')

ggsave("output/study2/Figures/figure5.png", plot = figure5, width = 8, height = 12.8)


### supplementary SFigure
Sfigure6 <- sigma2.u_plot_est_0.05 / 
  sigma2.s_plot_est_0.05 / 
  sigma2.n_plot_est_0.05 / 
  sigma2.p_plot_est_0.05 +
  plot_annotation(tag_levels='A')

ggsave("output/study2/SFig/Sfigure6_sigma2.png", plot = Sfigure6, width = 8, height = 12.5)



####################
# (2) Variance bias
####################

### Within study (u)
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



## Among study (s)
sigma2.s_plot_bias_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.s_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
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
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX(""),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_s$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




## Non-phylogeny effect (n)
sigma2.n_plot_bias_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.n == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.n_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX(""),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_n$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.n_plot_bias_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.n == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.n_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col6mods, guide="none")+
  scale_fill_manual(values=alpha(col6mods, 0.4), guide="none") +
  labs(title = TeX(""),
       x = "",
       y = TeX("$\\hat{\\sigma}^2_n$ bias"))+
  geom_boxplot(width=0.4)+
  geom_hline(aes(yintercept = 0), colour="darkgray")+ 
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





####################
# (2) Variance MSE
####################

## Within study (u)
sigma2.u_plot_mse_0.05 <-
  res.c |> 
  filter(sigma2.u == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.u_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
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
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_u$=0.3"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_u$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## Among study (s)
sigma2.s_plot_mse_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.s == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
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
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s$=0.3"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_s$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




## Non-phylogeny effect (n)
sigma2.n_plot_mse_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.n == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_n$=0.05"), 
       x = "",
       y = TeX("$\\hat{\\sigma}^2_n$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.n_plot_mse_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.n == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_s$=0.3"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_s$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



## Phylogeny effect (p)
sigma2.p_plot_mse_0.05 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.p == 0.05) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_p$=0.05"), 
       x = "",
       y = TeX("$\\hat{\\sigma}^2_p$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


sigma2.p_plot_mse_0.3 <-
  res.c |> 
  filter(!model_type %in% c("FE", "RE") & sigma2.p == 0.3) |>  
  ggplot(aes(x=factor(model_type), y=sigma.s_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = TeX("$\\sigma^2_p$=0.3"), 
       x = "", 
       y = TeX("$\\hat{\\sigma}^2_p$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




### supplementary SFigure 7
Sfigure7 <- sigma2.u_plot_mse_0.05 / 
  sigma2.s_plot_mse_0.05 / 
  sigma2.n_plot_mse_0.05 / 
  sigma2.p_plot_mse_0.05 +
  plot_annotation(tag_levels='A')

ggsave("output/study2/SFig/Sfigure7_sigma2_mse.png", plot = Sfigure7, width = 8, height = 12.5)



### supplementary SFigure 8
Sfigure8 <- sigma2.u_plot_mse_0.3 / 
  sigma2.s_plot_mse_0.3 / 
  sigma2.n_plot_mse_0.3 / 
  sigma2.p_plot_mse_0.3 +
  plot_annotation(tag_levels='A')

ggsave("output/study2/SFig/Sfigure8_sigma2_mse.png", plot = Sfigure8, width = 8, height = 12.5)





####################
# (3) Total variance estimate
####################


#### All sigma.2 conditions
sigma2.T_plot_bias_all <-
  res.c |> 
  filter(!model_type %in% c("FE")) |>  
  ggplot(aes(x=factor(model_type), y=sigma2.total_bias, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
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
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title = "", x = "", y = TeX("$\\hat{\\sigma}^2_{total}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sigma2.T_plot_mse_all




### supplementary SFigure 9
Sfigure9 <- sigma2.T_plot_bias_all / 
             sigma2.T_plot_mse_all + 
        plot_annotation(tag_levels='A') 


ggsave("output/study2/SFig/Sfigure9_sigma2_total.png", plot = Sfigure9, width = 9, height = 9.5)



