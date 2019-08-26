# R script to conduct gene enrichment analysis (with clusterProfiler) on genes found by SSLB 

source('../SSLB_functions.R', echo=TRUE)

library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

load("SSLB_result/SSLB_output.RData")

genes <- read.table(file = "data/gene_info.txt", stringsAsFactors = F)
genes <- as.matrix(genes)

X_SSLB <- out$X
B_SSLB <- out$B
K_SSLB <- out$K

reorder <- union(c(4, 10,
                   5, 27, 26, 19, 29, 24,
                   3, 16, 21, 33,
                   1, 7, 14, 6, 22, 23,
                   2, 28, 15, 20, 11, 40, 
                   32, 13, 39, 42, 41,
                   8, 34, 31, 30, 25,
                   9, 18, 12,
                   35, 43, 55, 63,
                   17, 38, 43), 1:out$K)

X_SSLB <- X_SSLB[, reorder]
B_SSLB <- B_SSLB[, reorder]

SSLB_hugo_list <- vector("list", K_SSLB)
SSLB_entrez_list <- vector("list", K_SSLB)
SSLB_enrichment <- vector("list", K_SSLB)

# convert HUGO gene ids to entrez for clusterprofiler
for (k in 1:K_SSLB) {
  hugo_symbol <- genes[which(B_SSLB[, k] != 0)]
  SSLB_hugo_list[[k]] <- hugo_symbol
  
  entrez_id <- mapIds(org.Mm.eg.db, hugo_symbol, 'ENTREZID', 'SYMBOL')
  SSLB_entrez_list[[k]] <- as.numeric(entrez_id)
}

SSLB_hugo <- sapply(SSLB_hugo_list, '[', seq(max(lengths(SSLB_hugo_list))))
SSLB_hugo[is.na(SSLB_hugo)] <- ""
write.table(SSLB_hugo,file="SSLB_result/SSLB_hugo.csv", quote=F,sep=",",row.names=F)

entrez_all <- mapIds(org.Mm.eg.db, as.vector(genes), 'ENTREZID', 'SYMBOL')

# Gene enrichment
for (k in 1:K_SSLB) {
  enriched <- enrichGO(SSLB_entrez_list[[k]], ont = "BP", OrgDb = org.Mm.eg.db, qvalueCutoff = 0.05, universe = entrez_all)
  SSLB_enrichment[[k]] <- enriched
}

save(SSLB_enrichment, file = "SSLB_result/GO_enrichment.RData")

n_enriched <- numeric(K_SSLB)
for (k in 1:K_SSLB) {
  n_enriched[k] <- length(SSLB_enrichment[[k]]$Description)
}

# number of enriched biclusters
sum(n_enriched != 0)/K_SSLB

# Plots --------------------------
k <- 44
enrich_simple <- simplify(SSLB_enrichment[[k]], cutoff = 0.9)

pdf(file = "figures/bic44_enrichment.pdf", width = 11, height = 7)

emapplot(enrich_simple, showCategory = 30)

dev.off()

pdf(file = "figures/bic44_enrichment.pdf", width = 10, height = 8)

emapplot(SSLB_enrichment[[k]], showCategory = 25)

dev.off()
