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
 
N = ncol(Y)
G = nrow(Y)

write.table(t(Y), file ="data/Y_raw.txt", row.names = FALSE, col.names = FALSE)
write.table(Y, file ="data/Y_raw_bicmix.txt", row.names = FALSE, col.names = FALSE)

save(gene_info, file = "data/gene_info.RData")
save(sample_info, file = "data/sample_info.RData")

# Normalize quantiles (average reference)
Y_normalized <- normalize.quantiles(Y)

# Normalize quantiles (gaussian reference)
Y_std_normal <- normalize.quantiles.use.target(Y, target = qnorm(seq(0.000001, 1-0.000001, length.out = G)))

Y_normalized <- t(Y_normalized)
Y_std_normal = t(Y_std_normal)

# Write data matrix to file
# Y is N x G matrix
# t(Y) is G x N matrix

write.table(t(Y_normalized), file = "data/Y_bicmix_quantile.txt", row.names = FALSE, col.names = FALSE)
write.table(Y_normalized, file ="data/Y_quantile.txt", row.names = FALSE, col.names = FALSE)

write.table(Y_std_normal, file ="data/Y_norm.txt", row.names = FALSE, col.names = FALSE)
write.table(t(Y_std_normal), file ="data/Y_bicmix_norm.txt", row.names = FALSE, col.names = FALSE)

