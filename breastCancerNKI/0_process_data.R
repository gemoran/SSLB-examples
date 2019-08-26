# R script to process data for breast cancer NKI dataset

library(breastCancerNKI)
library(impute)
library(preprocessCore)

# Load data from breastCancerNKI package
data(nki)
nki_exprs <- exprs(nki)
gene_info <- fData(nki)
sample_info <- pData(nki)

# Remove genes with > 10% missing
nki_missing <- apply(nki_exprs, 1, function(x) sum(is.na(x)) > 0.1 * ncol(nki_exprs))
nki_exprs <- nki_exprs[!nki_missing, ]

gene_info <- gene_info[!nki_missing, ]
  
# Impute remaining missing values using the impute package
Y <- impute.knn(nki_exprs)$data
rownames(Y) <- NULL
colnames(Y) <-  NULL
 
write.table(t(Y), file ="data/Y_raw.txt", row.names = FALSE, col.names = FALSE)
write.table(Y, file ="data/Y_raw_bicmix.txt", row.names = FALSE, col.names = FALSE)

save(gene_info, file = "data/gene_info.RData")
save(sample_info, file = "data/sample_info.RData")
