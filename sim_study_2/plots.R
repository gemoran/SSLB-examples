# R script to create consensus and relevance/recovery plots

library(reshape2)
library(tidyverse)

load( "results/sim_study_2.RData")

result_melt <- melt(result)
result_tibble <- as_tibble(result_melt)

dat <- spread(result_tibble, Var2, value)

# consensus 
dat <- filter(dat, !(Var1 %in% c("FABIA_truth","Bicmix_truth", "FABIAS", "FABIAS_truth", "CCS")))
labels <- unique(dat$Var1)
labels <- as.character(labels)
labels[which(labels == "SSLB_PY")] <- "SSLB-PY"
labels[which(labels == "SSLB_IBP")] <- "SSLB-IBP"
labels[which(labels == "SSLB_BB")] <- "SSLB-BB"
labels[which(labels == "Bicmix")] <- "BicMix"


g <- ggplot(dat, aes(Var1, consensus))
g <- g + geom_boxplot(aes(fill = Var1))
g <- g + scale_y_continuous(limits = c(0, 1))
g <- g + labs(x = "", y = "Consensus")
g <- g + scale_x_discrete(labels = labels)
g <- g + theme(axis.text.x=element_text(color = "black", size=15, angle=0)) 
g <- g + theme(legend.position="none")
g <- g + theme(axis.ticks = element_blank())
g <- g + theme(text = element_text(size=20))
g

pdf(file = "figures/consensus.pdf", width = 6, height = 3.5, pointsize = 14)
g
dev.off()

# relevance and recovery

g <- ggplot(dat, aes(recovery, relevance))
g <- g + geom_point(aes(color = Var1, shape = Var1))
g <- g + scale_x_continuous(limits = c(0, 1))
g <- g + scale_y_continuous(limits = c(0, 1))
g <- g + labs(x = "Recovery", y = "Relevance")
g <- g + scale_color_discrete(labels = labels)
g <- g + scale_shape_discrete(labels = labels)
g <- g + theme(legend.title=element_blank()) 
g <- g + theme(text = element_text(size=20))
g


pdf(file = "figures/rr.pdf", width = 6, height = 4, pointsize = 14)
g
dev.off()



