load("data/gene_info.RData")
load("data/sample_info.RData")
load("ISA_result/out_ISA.RData")

gene_info[] <- lapply(gene_info, as.character)

er <- sample_info$er
er[er == 0] <- -1

er_order <- order(er)
er_factor <- factor(er)

#---------

X_SSBiEM = as.matrix(read.csv("matlab/X.txt", header = F))
B_SSBiEM = as.matrix(read.csv("matlab/B.txt", header = F))

X_SSBiEM = X_SSBiEM[, c(1:14)]
B_SSBiEM = B_SSBiEM[, c(1:14)]

pdf(file = "figures/X_SSBiEM.pdf", width = 5, height = 4)

plot_matrix(X_SSBiEM[er_order,],
            title = "SSBiEM Factor Matrix (K = 2)") +
  theme(text = element_text(size= 12))

dev.off()


X_ISA = out_ISA$rows
B_ISA = out_ISA$columns

K_ISA = ncol(X_ISA)

er_plot <- plot_matrix(er[er_order], title = "ER") +
  theme(text = element_text(size= 12))

isa_plot <- plot_matrix(X_ISA[er_order, ], 
                           title = "ISA Factor Matrix (K = 540)", legend = F) +
  theme(text = element_text(size= 12)) 


pdf(file = "figures/X_ISA.pdf", width = 10, height = 4)

grid.arrange(er_plot, isa_plot, nrow = 1, widths = c(1.5, 30))

dev.off()

p_vals <- numeric(ncol(X_ISA))
for (k in 1:ncol(X_ISA)) {
  p_vals[k] <- wilcox.test(X_ISA[, k] ~ er_factor)$p.value
}


enrich_list <- vector("list", K_ISA)
n_enrich <- numeric(K_ISA)
for (i in 1:K_ISA) {
  test <- gene_info$EntrezGene.ID[which(B_ISA[, i] != 0)]
  enrich_list[[i]] <- enrichGO(test, ont = "BP", OrgDb = org.Hs.eg.db, qvalueCutoff = 0.05)
  n_enrich[i] <- length(enrich_list[[i]]$Description)
}


save(enrich_list, file = "ISA_result/enrichment.RData")

