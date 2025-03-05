# load libraries
library(lme4); library(nlme); library(parameters); library(metafor)
library(ape); library(MASS)

# load data (example simdat from study2.sub)
load("~/Projects/sim-meta-eco-evo/data/simdat_1.RDATA")

k <- length(unique(dat$id))
k.species <- length(unique(dat$species.id))
k.studies <- length(unique(dat$study))
dat$species.phylo <- dat$species.id

# simulate tree
tree <- rtree(k.species, tip.label=seq_len(k.species))
tree <- compute.brlen(tree, power = 1)
P <- vcv(tree, corr=TRUE)
P <- P[order(as.numeric(rownames(P))), order(as.numeric(rownames(P)))]


# fit rma.mv model
mod.rma <- rma.mv(yi, vi, random = list(~ 1 | species.id, ~ 1 | species.phylo, ~ 1 | study.id, ~ 1| id),
                  #mods = ~ factor(x1) + x2 + factor(x3),
                  R = list(species.phylo=P),
                  data = dat,
                  test = "t", dfs = "contain")
df.rma <- mod.rma$ddf[[1]]


# fit simple lmer model
mod <- lmer(yi ~ 1 + factor(x1) + x2 + factor(x3) +
              (1|study.id) + 
              (1|species.id),
            data = dat)

model_parameters(mod, df_method = "satterthwaite")
dof_kenward(mod)
dof_satterthwaite(mod)
dof_betwithin(mod)





### (1) Use df from lmer output on simplistic model (without phylogenetic random effect)
df.1 <- dof_satterthwaite(mod)[[1]]


### (2) Approximate df given advice from Gotelli and Ellison (2004) + Satterwaith (1946)
# using harmonic mean of each parameter to get overall df
df.2 <- 3/((1/k) + (1/k.species) + (1/k.species))
df.3 <- 4/((1/k.studies) + (1/k.species) + (1/k.species) + (1/k))

## Compare three df estimates
data.frame(model = c("rma", "lmer", "harmonic_mean"),
           df = c(df.rma, df.1, df.2))


## Non-parametric (cluster wild bootstrap & permutations, sampling) - design based 

### (3) Cluster wild bootstrap
