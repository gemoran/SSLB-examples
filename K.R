library(reshape2)
library(tidyverse)
library(xtable)

load("sim_study_1/results/sim_study_1.RData")

result_melt <- melt(result)
result_tibble <- as.tibble(result_melt)
dat <- spread(result_tibble, Var2, value)
nrep <- dim(result)[3]

dat <- filter(dat, !(Var1 %in% c("SSBiEM", "FABIA")))
labels <- unique(dat$Var1)
labels <- as.character(labels)

# sim_1

dat <- group_by(dat, Var1)
K_sim_1 <- summarize(dat, K = paste(format(round(mean(K), 1), nsmall = 1), " {\\scriptsize (", format(round(sd(K)/sqrt(nrep), 2), nsmall = 2), ")}", sep = ""))

# sim_2
load("sim_study_2/results/sim_study_2.RData")

result_melt <- melt(result)
result_tibble <- as.tibble(result_melt)
dat <- spread(result_tibble, Var2, value)
nrep <- dim(result)[3]

dat <- filter(dat, !(Var1 %in% c("SSBiEM", "FABIA")))
labels <- unique(dat$Var1)
labels <- as.character(labels)

dat <- group_by(dat, Var1)
K_sim_2 <- summarize(dat, K = paste(format(round(mean(K), 1), nsmall = 1), " {\\scriptsize (", format(round(sd(K)/sqrt(nrep), 2), nsmall = 2), ")}", sep = ""))

K_sim <- cbind(K_sim_1, K_sim_2$K)

print(xtable((K_sim)), include.rownames = F, sanitize.text.function = function(x) {x})

#---------------------

load("sim_study_3/results/sim_study_3.RData")

result_melt <- melt(result)
result_tibble <- as.tibble(result_melt)
dat <- spread(result_tibble, Var2, value)
nrep <- dim(result)[3]

dat <- filter(dat, !(Var1 %in% c("SSBiEM", "FABIA")))
labels <- unique(dat$Var1)
labels <- as.character(labels)

# sim_3

dat <- group_by(dat, Var1)
K_sim_3 <- summarize(dat, K = paste(format(round(mean(K), 1), nsmall = 1), " {\\scriptsize (", format(round(sd(K)/sqrt(nrep), 2), nsmall = 2), ")}", sep = ""))

# sim_4

load("sim_study_4/results/sim_study_4.RData")

result_melt <- melt(result)
result_tibble <- as.tibble(result_melt)
dat <- spread(result_tibble, Var2, value)
nrep <- dim(result)[3]

dat <- filter(dat, !(Var1 %in% c("SSBiEM", "FABIA")))
labels <- unique(dat$Var1)
labels <- as.character(labels)

dat <- group_by(dat, Var1)
K_sim_4 <- summarize(dat, K = paste(format(round(mean(K), 1), nsmall = 1), " {\\scriptsize (", format(round(sd(K)/sqrt(nrep), 2), nsmall = 2), ")}", sep = ""))

K_sim <- cbind(K_sim_3, K_sim_4$K)

print(xtable((K_sim)), include.rownames = F, sanitize.text.function = function(x) {x})


