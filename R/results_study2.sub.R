################################################################################
# Plot results
################################################################################

load("~/PhD/2_MetaRVE/results/sims/all_results_study2.RDATA") ### change this to OSF or repo where i'll store large files


###
# load librairies
library(pacman)
p_load(tidyverse, ggplot2, broom, here,
       cowplot, ggdist, patchwork,
       scales, latex2exp)

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
########

# derive coverage rate for model fixed coefficients
results <- results |>
  mutate(b1.1_cov =  b1.1_ci_lb < b1.1 &  b1.1_ci_ub > b1.1,
         b1.2_cov =  b1.2_ci_lb < b1.2 &  b1.2_ci_ub > b1.2,
         b2_cov =  b2_ci_lb < b2 &  b2_ci_ub > b2,
         b3_cov =  b3_ci_lb < b3 &  b3_ci_ub > b3)


# derive confidence interval width for model coefficients
results <- results |>
  mutate(mu_ci_width = mu_ci_ub - mu_ci_lb,
         b1.1_ci_width = b1.1_ci_ub - b1.1_ci_lb,
         b1.2_ci_width = b1.2_ci_ub - b1.2_ci_lb,
         b2_ci_width = b2_ci_ub - b2_ci_lb,
         b3_ci_width = b3_ci_ub - b3_ci_lb)

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