################################################################################
# Plot results
################################################################################

load("~/PhD/2_MetaRVE/results/sims/collated_sim_pilot_results.RDATA") # load pilot results

results <- dat


###
# load librairies
library(pacman)
p_load(tidyverse, ggplot2, broom, here,
       cowplot, ggdist, patchwork,
       scales, latex2exp)

# create new variable for plots: model_type
results$model_type <- gsub("-CR1|-CR0", "", results$model)


# reorder factor levels for model_type and model
results$model <- factor(results$model, levels=c("PML", "PML-CR0", "PML-CR1",
                                                "PML-VCV-0.2", "PML-VCV-0.2-CR0", "PML-VCV-0.2-CR1",
                                                "PML-VCV-0.5", "PML-VCV-0.5-CR0", "PML-VCV-0.5-CR1",
                                                "PML-VCV-0.8", "PML-VCV-0.8-CR0", "PML-VCV-0.8-CR1"))


results$model_type <- factor(results$model_type, levels=c("PML", "PML-VCV-0.2", "PML-VCV-0.5", "PML-VCV-0.8"))


results$CR_method <- factor(results$CR_method, levels=c("none", "CR0", "CR1", "CR2"))



########
# derive coverage rate for model fixed coefficients
results <- results |>
  mutate(b1.1_cov =  b1.1_ci_lb < b1.1 &  b1.1_ci_ub > b1.1,
         b1.2_cov =  b1.2_ci_lb < b1.2 &  b1.2_ci_ub > b1.2,
         b2_cov =  b2_ci_lb < b2 &  b2_ci_ub > b2,
         b3_cov =  b3_ci_lb < b3 &  b3_ci_ub > b3)


# derive MSE for model fixed coefficients
results <- results |>
  mutate(b1.1_cov =  b1.1_ci_lb < b1.1 &  b1.1_ci_ub > b1.1,
         b1.2_cov =  b1.2_ci_lb < b1.2 &  b1.2_ci_ub > b1.2,
         b2_cov =  b2_ci_lb < b2 &  b2_ci_ub > b2,
         b3_cov =  b3_ci_lb < b3 &  b3_ci_ub > b3)


# derive confidence interval width for model coefficients
results <- results |>
  mutate(b1.1_mse = mean((b1.1_est - b1.1)^2),
         b1.2_mse = mean((b1.2_est - b1.2)^2),
         b2_mse = mean((b2_est - b2)^2),
         b3_mse = mean((b3_est - b3)^2)
)

# derive RMSE
results <- results |>
  group_by(model, scenario) |>
  mutate(b0_rmse = sqrt(mean((b0_est - b0)^2)))


# derive sample variance of b0 estimate
sample_var <- results |> 
  group_by(model, model_type, CR_method) |>
  summarise(S2 = sum((b0_est - mean(b0_est))^2) / (length(b0_est) - 1)) |> 
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



## data setup - filter out CRVE methods and only keep VCV-0.5 for plots
res.a <- results |> 
  filter(CR_method %in% c("none")) 


# derive coverage for this subset across all varying conditions
cov.a <- res.a |>
  group_by(model_type, rho, sigma2.s,
           sigma2.u, sigma2.p, sigma2.n,
           k.studies, k.species) |> 
  summarise(cov_prop_b0 = mean(b0_cov, na.rm = TRUE),
            cov_prop_b1.1 = mean(b1.1_cov, na.rm = TRUE),
            cov_prop_b1.2 = mean(b1.2_cov, na.rm = TRUE),
            cov_prop_b2 = mean(b2_cov, na.rm = TRUE),
            cov_prop_b3 = mean(b3_cov, na.rm = TRUE))




####################
# (1) Bias 
####################

b0_est_plot <-
  ggplot(res.a, aes(x=factor(model_type), y=b0_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="", x="", y = TeX("$\\hat{\\b0}$"))+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=b0, colour="darkgray", linewidth=0.5)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b1.1_est_plot <-
  ggplot(res.a, aes(x=factor(model_type), y=b1.1_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="", x="", y = TeX("$\\hat{\\b1.1}$"))+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=b1.1, colour="darkgray", linewidth=0.5)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b1.2_est_plot <-
  ggplot(res.a, aes(x=factor(model_type), y=b1.2_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="", x="", y = TeX("$\\hat{\\b1.2}$"))+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=b1.2, colour="darkgray", linewidth=0.5)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b2_est_plot <-
  ggplot(res.a, aes(x=factor(model_type), y=b2_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="", x="", y = TeX("$\\hat{\\b2}$"))+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=b2, colour="darkgray", linewidth=0.5)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



b3_est_plot <-
  ggplot(res.a, aes(x=factor(model_type), y=b3_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="", x="", y = TeX("$\\hat{\\b3}$"))+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=b3, colour="darkgray", linewidth=0.5)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





####################
# (2) MSE
####################


b0_mse_plot.a <-
  ggplot(res.a, aes(x=factor(model_type), y=b0_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title="",x="", y = TeX("$\\hat{\\b0}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2/b0_mse_all.png", plot = b0_mse_plot.a, width = 10, height = 6)




b1.1_mse_plot <-
  ggplot(res.a, aes(x=factor(model_type), y=b1.1_mse, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title="",x="", y = TeX("$\\hat{\\b1.1}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/b1.1_mse.png", plot = b1.1_mse_plot, width = 10, height = 6)




b1.2_mse_plot <-
  ggplot(res.a, aes(x=factor(model_type), y=b1.2_mse, color=model_type, fill=model_type)) + 
  #stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title="",x="", y = TeX("$\\hat{\\b1.2}$ MSE"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/b1.2_mse.png", plot = b1.2_mse_plot, width = 10, height = 6)




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

cov_plot.b0 <-
  ggplot(cov.a, aes(x=factor(model_type), y=cov_prop_b0, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\b0}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


cov_plot.b1.1 <-
  ggplot(cov.a, aes(x=factor(model_type), y=cov_prop_b1.1, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\beta_{1.1}}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b1.1.png", plot = cov_plot.b1.1, width = 10, height = 6)


cov_plot.b1.2 <-
  ggplot(cov.a, aes(x=factor(model_type), y=cov_prop_b1.2, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\beta_{1.2}}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b1.2.png", plot = cov_plot.b1.2, width = 10, height = 6)


cov_plot.b2 <-
  ggplot(cov.a, aes(x=factor(model_type), y=cov_prop_b2, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\beta_{2}}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b2.png", plot = cov_plot.b2, width = 10, height = 6)


cov_plot.b3 <-
  ggplot(cov.a, aes(x=factor(model_type), y=cov_prop_b3, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\beta_{3}}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(~rho, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b3.png", plot = cov_plot.b3, width = 10, height = 6)






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
ggsave("output/study2/Fig1/ci_width_mu.png", plot = ci_plot.a, width = 10, height = 6)



ci_plot.b1.1 <-
  ggplot(res.a, aes(x=factor(model_type), y=b1.1_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Confidence interval width of $\\hat{\\beta_{1.1}}$"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/ci_width_b1.1.png", plot = ci_plot.b1.1, width = 10, height = 6)



ci_plot.b1.2 <-
  ggplot(res.a, aes(x=factor(model_type), y=b1.2_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title="",x="", y = TeX("Confidence interval width of $\\hat{\\beta_{1.1}}$"))+
  geom_boxplot(width=0.4)+
  facet_wrap(~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/ci_width_b1.2.png", plot = ci_plot.b1.2, width = 10, height = 6)

