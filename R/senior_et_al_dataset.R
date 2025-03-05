library(readr); library(ggplot2); library(patchwork)

# Senior, A. M., Grueber, C. E., Kamiya, T., Lagisz, M., O'Dwyer, K., Santos, E. S., & Nakagawa, S. (2016). Heterogeneity in ecological and evolutionary meta-analyses: its magnitude and implications. Ecology, 97(12), 3293–3299. https://doi.org/10.1002/ecy.1591 

# load meta-analysis data
# downloaded from: https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.1591&file=ecy1591-sup-0004-DataS1.zip
dat <- read_csv("data/Data_Package_Part_4_Meta-Analysis_data.csv")


# get histogram of number effect sizes per study
hist.k <- ggplot(dat, aes(K)) +
  geom_histogram(binwidth=20, fill="lightgrey", color="black") + 
  scale_x_continuous(breaks=seq(0, 800, 100)) +
  scale_y_continuous(breaks=seq(0, 27, 2)) +
  labs(title="", x="Number of effect sizes per study", y="Frequency") +
  theme_bw()

# summarise mean, median, mode, max, min of effect sizes
mean(dat$K)
median(dat$K)
as.numeric(names(sort(-table(dat$K))[1]))
max(dat$K)
min(dat$K)

# get histogram of number effect sizes per study
hist.sp <- ggplot(dat, aes(Species)) +
  geom_histogram(bins=15, fill="lightgrey", color="black", na.rm = TRUE) + 
  labs(title="", x="Number of species per study", y="Frequency") +
  theme_bw()



# get hist.plot 
hist.plot <- hist.k + hist.sp + plot_annotation(tag_levels='A')

# save plot
ggsave(hist.plot, file="output/senior_et_al_hist.png", width=10, height=6, dpi=300)
