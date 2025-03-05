################################################################################
# Plot results
################################################################################

load("~/PhD/2_MetaRVE/results/sims/all_results_study2.sub.RDATA")
results <- combined_dat


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


# derive MSE for model fixed coefficients
results <- results |>
  mutate(b1_mse = (b1_est - b1)^2,
         b2_mse = (b2_est - b2)^2,
         b3_mse = (b3_est - b3)^2
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
                            `0.2`=parse(text=TeX("$\\phi=0.2$")),
                            `0.5`=parse(text=TeX("$\\phi=0.5$")),
                            `0.8`=parse(text=TeX("$\\phi=0.8$"))
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

# derive coverage for this subset across all varying conditions
cov <- results |>
  group_by(model_type, rho, sigma2.s,
           sigma2.u, sigma2.p, sigma2.n,
           k.studies, k.species, rho_lab,
           CR_method) |> 
  summarise(cov_prop_b0 = mean(b0_cov, na.rm = TRUE),
            cov_prop_b1 = mean(b1_cov, na.rm = TRUE),
            cov_prop_b2 = mean(b2_cov, na.rm = TRUE),
            cov_prop_b3 = mean(b3_cov, na.rm = TRUE))


####################
# b0 coeff 
####################

b0_est_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b0_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("$\\hat{\\beta_{0}}$"), x="", y = "")+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=results$b0[1], colour="darkgray", linewidth=0.5)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b0_mse_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b0_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title=TeX("$\\hat{\\beta_{0}}$ MSE"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/b0_mse_b.png", plot = b0_mse_plot, width = 10, height = 6)


cov_plot.b0 <- cov |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=cov_prop_b0, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  scale_y_continuous(limits=c(0.85, 1), breaks=seq(0.85, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{0}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/cov_b0_b.png", plot = cov_plot.b0, width = 10, height = 6)


ci_plot.b0 <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b0_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("Confidence interval width of $\\hat{\\beta_{0}}$"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/ci_width_b0_b.png", plot = ci_plot.b0, width = 10, height = 6)


### Combine into Figure S.b0
fig.beta.0 <- b0_est_plot + b0_mse_plot + cov_plot.b0 + ci_plot.b0 +
  plot_annotation(tag_levels='A')

ggsave("output/study2.sub/Figures/figure_beta_0.png", plot = fig.beta.0, width = 11, height = 11)





####################
# b1 coeff 
####################

b1_est_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b1_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("$\\hat{\\beta_{1}}$"), x="", y = "")+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=results$b1[1], colour="darkgray", linewidth=0.5)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b1_mse_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b1_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title=TeX("$\\hat{\\beta_{1}}$ MSE"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/b1_mse_b.png", plot = b1_mse_plot, width = 10, height = 6)


cov_plot.b1 <- cov |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=cov_prop_b1, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  scale_y_continuous(limits=c(0.85, 1), breaks=seq(0.85, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{1}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/cov_b1_b.png", plot = cov_plot.b1, width = 10, height = 6)


ci_plot.b1 <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b1_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("Confidence interval width of $\\hat{\\beta_{1}}$"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/ci_width_b1_b.png", plot = ci_plot.b1, width = 10, height = 6)


### Combine into Figure S.b1
fig.b1 <- b1_est_plot + b1_mse_plot + cov_plot.b1 + ci_plot.b1 +
  plot_annotation(tag_levels='A')

ggsave("output/study2.sub/Figures/figure_b1.png", plot = fig.b1, width = 11, height = 11)




####################
# b2 coeff 
####################

b2_est_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b2_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("$\\hat{\\beta_{2}}$"), x="", y = "")+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=results$b2[1], colour="darkgray", linewidth=0.5)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b2_mse_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b2_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title=TeX("$\\hat{\\beta_{2}}$ MSE"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/b2_mse_b.png", plot = b2_mse_plot, width = 10, height = 6)


cov_plot.b2 <- cov |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=cov_prop_b2, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  scale_y_continuous(limits=c(0.85, 1), breaks=seq(0.85, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{2}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/cov_b2_b.png", plot = cov_plot.b2, width = 10, height = 6)


ci_plot.b2 <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b2_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("Confidence interval width of $\\hat{\\beta_{2}}$"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/ci_width_b2_b.png", plot = ci_plot.b2, width = 10, height = 6)


### Combine into Figure S.b2
fig.b2 <- b2_est_plot + b2_mse_plot + cov_plot.b2 + ci_plot.b2 +
  plot_annotation(tag_levels='A')

ggsave("output/study2.sub/Figures/figure_b2.png", plot = fig.b2, width = 11, height = 11)



####################
# b3 coeff 
####################

b3_est_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b3_est, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("$\\hat{\\beta_{3}}$"), x="", y = "")+
  #scale_y_continuous(limit=c(-0.6, 1), breaks=round(seq(-0.6, 1, 0.2), 3)) +
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=results$b3[1], colour="darkgray", linewidth=0.5)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


b3_mse_plot <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b3_mse, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limit=c(0, 0.6), breaks=round(seq(0, 0.6, 0.2), 3)) +
  labs(title=TeX("$\\hat{\\beta_{3}}$ MSE"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/b3_mse_b.png", plot = b3_mse_plot, width = 10, height = 6)


cov_plot.b3 <- cov |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=cov_prop_b3, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  scale_y_continuous(limits=c(0.85, 1), breaks=seq(0.85, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{3}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/cov_b3_b.png", plot = cov_plot.b3, width = 10, height = 6)


ci_plot.b3 <- results |> 
  filter(CR_method=="none") |>
  ggplot(aes(x=factor(model_type), y=b3_ci_width, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  labs(title=TeX("Confidence interval width of $\\hat{\\beta_{3}}$"),x="", y = "")+
  geom_boxplot(width=0.4)+
  facet_grid(cols=vars(rho_lab),
             #rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/ci_width_b3_b.png", plot = ci_plot.b3, width = 10, height = 6)


### Combine into Figure S.b3
fig.b3 <- b3_est_plot + b3_mse_plot + cov_plot.b3 + ci_plot.b3 +
  plot_annotation(tag_levels='A')

ggsave("output/study2.sub/Figures/figure_b3.png", plot = fig.b3, width = 11, height = 11)







################################################################################
################################################################################

# Including CRVE methods


cov_plot.b0_crve <- cov |> 
  ggplot(aes(x=factor(model_type), y=cov_prop_b0, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limits=c(0.8, 1), breaks=seq(0.8, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{0}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b0_crve.png", plot = cov_plot.b0_crve, width = 9, height = 11)


cov_plot.b1_crve <- cov |> 
  ggplot(aes(x=factor(model_type), y=cov_prop_b1, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limits=c(0.85, 1), breaks=seq(0.85, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{1}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b1_crve.png", plot = cov_plot.b1_crve, width = 9, height = 11)



cov_plot.b2_crve <- cov |> 
  ggplot(aes(x=factor(model_type), y=cov_prop_b2, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limits=c(0.85, 1), breaks=seq(0.85, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{2}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b2_crve.png", plot = cov_plot.b2_crve, width = 9, height = 11)




cov_plot.b3_crve <- cov |> 
  ggplot(aes(x=factor(model_type), y=cov_prop_b3, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  scale_y_continuous(limits=c(0.85, 1), breaks=seq(0.85, 1, 0.05))+
  labs(title= TeX("Coverage rate of $\\hat{\\beta_{3}}$"),x="", y = TeX(""))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_grid(cols=vars(rho_lab),
             rows=vars(CR_method),
             labeller=label_parsed) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/study2.sub/cov_b3_crve.png", plot = cov_plot.b3_crve, width = 9, height = 11)







################################################################################

# Study size and coverage rate 
cov |> 
  filter(CR_method=="none") |>
ggplot(aes(x=factor(model_type), y=cov_prop_b1, color=model_type, fill=model_type)) + 
  stat_halfeye() +
  scale_color_manual(values=col5mods, guide="none")+
  scale_fill_manual(values=alpha(col5mods, 0.4), guide="none") +
  #scale_y_continuous(limits=c(0.92, 0.97), breaks=seq(0.92, 0.97, 0.05))+
  labs(title="",x="", y = TeX("Coverage rate of $\\hat{\\beta_{1}}$"))+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  facet_wrap(k.studies~rho_lab, labeller=label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("output/study2.sub/cov_b1_b.png", plot = cov_plot.b1, width = 10, height = 6)
