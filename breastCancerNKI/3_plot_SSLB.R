# R script to plot SSLB biclusters

library(pals)

source('../SSLB_functions.R', echo=TRUE)

load("SSLB_result/SSLB_out.RData")
load("data/gene_info.RData")
load("data/sample_info.RData")

Y <- read.table("data/Y_raw.txt", stringsAsFactors = F)
Y <- as.matrix(Y)
N <- nrow(Y)
G <- ncol(Y)

gene_info[] <- lapply(gene_info, as.character)

er <- sample_info$er
er[er == 0] <- -1

er_order <- order(er)
er_factor <- factor(er)

X_SSLB <- out$X
K_SSLB <- out$K
B_SSLB <- out$B

# which biclusters contain ESR1
B_SSLB[which(gene_info$HUGO.gene.symbol == "ESR1"), ]

# HER2 markers
B_SSLB[which(gene_info$HUGO.gene.symbol == "ERBB2"), ]

# sign change to make ERRB2 positive in B for consistent colors
X_SSLB[, 2] <- -X_SSLB[, 2]
B_SSLB[, 2] <- -B_SSLB[, 2]

# order biclusters (columns)
X_SSLB <- X_SSLB[, union(c(3, 2), 1:K_SSLB)]
B_SSLB <- B_SSLB[, union(c(3, 2), 1:K_SSLB)]

# significant difference
p_vals <- numeric(K_SSLB)
for (k in 1:K_SSLB) {
  p_vals[k] <- wilcox.test(X_SSLB[, k] ~ er_factor)$p.value
}

er_k <- which.min(p_vals)
er_k < 0.01/K_SSLB

# PLOTS -----------------------------------------------------
# order patients within ER status
g1 <- which((X_SSLB[, 1] < 0 & er == -1) & X_SSLB[, 2] <= 0)
ord1 <- order(X_SSLB[g1, 1], decreasing = F)
ord1 <- g1[ord1]

g2 <- which((X_SSLB[, 1] < 0 & er == -1) & X_SSLB[, 2] > 0)
ord2 <- order(X_SSLB[g2, 1], decreasing = F)
ord2 <- g2[ord2]

g3 <- which((X_SSLB[, 1] == 0 & er == -1))
ord3 <- order(X_SSLB[g3, 2], decreasing = T)
ord3 <- g3[ord3]

g4 <- which(X_SSLB[, 1] < 0 & er == 1)
ord4 <- order(X_SSLB[g4, 2], decreasing = F)
ord4 <- g4[ord4]

g5 <- which(X_SSLB[, 1] >= 0 & er == 1)
ord5 <- order(X_SSLB[g5, 2], decreasing = T)
ord5 <- g5[ord5]

ord <- c(ord1, ord2, ord3, ord4, ord5)

# cap for viualization
X_SSLB[X_SSLB > 5] <- 5
X_SSLB[X_SSLB < -5] <- -5


# get non-zero genes
B_subset <- c(which(B_SSLB[, 1] > 0), which(B_SSLB[, 1] < 0))

Y_subset <- Y[ord, B_subset]

# set abs(Y) > rescale -> Y = rescale for plot
rescale <- 0.25
Y_subset[Y_subset > rescale] <- rescale
Y_subset[Y_subset < -rescale] <- -rescale

Y_plot <- plot_matrix(Y_subset, title = "Gene Expression Matrix", legend = F) +
  theme(text = element_text(size= 12)) 

er_plot <- plot_matrix(er[ord], title = "ER") +
  theme(text = element_text(size= 12))

# plot with first 10 biclusters
SSLB_plot <- plot_matrix(X_SSLB[ord, 1:10]/max(abs(X_SSLB)) * rescale, 
                         title = "SSLB Factor Matrix", legend = T) +
  theme(text = element_text(size= 12)) + theme(legend.position = "left")

png(file = "figures/SSLB_fig1_reorder.png", width = 10, height = 4, units = "in", res = 300)

grid.arrange(SSLB_plot, Y_plot, er_plot, nrow = 1, widths = c(10, 15, 1.25))

dev.off()


#---------------------------------------------------

# plot with all biclusters

SSLB_plot <- plot_matrix(X_SSLB[ord,]/max(abs(X_SSLB)) * rescale, 
                         title = "SSLB Factor Matrix") +
  theme(text = element_text(size= 12))

pdf(file = "figures/X_SSLB_all.pdf", width = 5, height = 4)

SSLB_plot

dev.off()
