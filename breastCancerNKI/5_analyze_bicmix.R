# R script to analyze/plot bicmix result

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pals)

load("data/gene_info.RData")
load("data/sample_info.RData")

option = c("raw", "norm")
i = 2

dir = paste("bicmix_result/", option[i], sep = "")

X_bicmix <- try(as.matrix(read.table(paste(dir, "/EX", sep = ""), header = F)), silent = TRUE)
if (inherits(X_bicmix, 'try-error')) {
  X_bicmix <- matrix(0, nrow = K_init, ncol = N)
}
X_bicmix <- t(X_bicmix)
B_bicmix <- try(as.matrix(read.table(paste(dir, "/LAM", sep = ""), header = F)), silent = TRUE)
if (inherits(B_bicmix, 'try-error')) {
  B_bicmix <- matrix(0, nrow = G, ncol = K_init)
}

if (any(X_bicmix != 0)) {
  
  X_norm <- apply(X_bicmix, 2, function(x) sum(abs(x)))
  B_norm <- apply(B_bicmix, 2, function(x) sum(abs(x)))
  
  d <- sqrt(X_norm/B_norm)
  X_bicmix <- t(1/d * t(X_bicmix))
  B_bicmix <- t(d * t(B_bicmix))
  
  X_bicmix[abs(X_bicmix) < 10^(-10)] <- 0
  
}

K_bicmix = ncol(X_bicmix)

gene_info[] <- lapply(gene_info, as.character)

er <- sample_info$er
er[er == 0] <- -1

er_order <- order(er)
er_factor <- factor(er)

er <- sample_info$er
er[er == 0] <- -1

er_order <- order(er)
er_factor <- factor(er)

# significant difference
p_vals <- numeric(K_bicmix)
for (k in 1:K_bicmix) {
  p_vals[k] <- wilcox.test(X_bicmix[, k] ~ er_factor)$p.value
}

er_k <- which.min(p_vals)
p_vals < 0.01/K_bicmix

X_bicmix = X_bicmix[, order(p_vals)]
B_bicmix = B_bicmix[, order(p_vals)]

# which biclusters contain ESR1
B_bicmix[which(gene_info$HUGO.gene.symbol == "ESR1"), ]

# HER2 markers
B_bicmix[which(gene_info$HUGO.gene.symbol == "ERBB2"), ]

# sign change to make ESR1 positive in B for consistent colors
X_bicmix[, 1] <- -X_bicmix[, 1]
B_bicmix[, 1] <- -B_bicmix[, 1]

#  Testing further gene markers (ZMZ14) --------------------------------
# luminal A (ER+, HER2-)
markers <- c("KRT8", "FOXA1", "XBP1", "GATA3", "ADH1B")
B_bicmix[which(gene_info$HUGO.gene.symbol %in% markers), ]

in_bic <- B_bicmix[which(gene_info$HUGO.gene.symbol %in% markers), 1]

found <- sum(in_bic != 0)
tot <- length(markers)
n_bic <- sum(B_bicmix[,1] != 0)

mat <- matrix(c(found, tot - found, n_bic - found, G - tot - (n_bic - found)), 2, 2)
fisher.test(mat)


# ENRICHMENT ---------------------------------

enrich_list <- vector("list", K_bicmix)
n_enrich <- numeric(K_bicmix)
for (i in 1:K_bicmix) {
  test <- gene_info$EntrezGene.ID[which(B_bicmix[, i] != 0)]
  enrich_list[[i]] <- enrichGO(test, ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05)
  n_enrich[i] <- length(enrich_list[[i]]$Description)
}

save(enrich_list, file = "bicmix_result/enriched_raw.RData")


er_plot <- plot_matrix(er[er_order], title = "ER") +
  theme(text = element_text(size= 12))

# plot with first 10 biclusters
bicmix_plot <- plot_matrix(X_bicmix[er_order, ], 
                         title = "BicMix Factor Matrix (K = 13)", legend = F) +
  theme(text = element_text(size= 12)) 

pdf(file = "figures/X_bicmix.pdf", width = 5, height = 4)

grid.arrange(er_plot, bicmix_plot, nrow = 1, widths = c(1.5, 15))

dev.off()
