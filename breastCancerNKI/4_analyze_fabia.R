# R script to analyze/plot fabia result

load("data/gene_info.RData")
load("data/sample_info.RData")

load("fabia_result/fabia_out_50.RData")

X_FABIA <- t(out_FABIA@Z)
B_FABIA <- out_FABIA@L

FABIA_biclusters <- extractBic(out_FABIA)
B_FABIA_support <- matrix(0, nrow = G, ncol = K_FABIA)
X_FABIA_support <- matrix(0, nrow = N, ncol = K_FABIA)

for(k in 1:K_FABIA) {
  B_FABIA_support[FABIA_biclusters$numn[k,]$numng, k] <- 1
  X_FABIA_support[FABIA_biclusters$numn[k,]$numnp, k] <- 1
}

X_FABIA[X_FABIA_support == 0] <- 0
B_FABIA[B_FABIA_support == 0] <- 0

zeroes_B <- which(apply(B_FABIA, 2, function(x) all(x == 0)))
zeroes_X <- which(apply(X_FABIA, 2, function(x) all(x == 0)))
zeroes <- union(zeroes_B, zeroes_X)

if (length(zeroes) > 0) {
  X_FABIA <- as.matrix(X_FABIA[, -zeroes])
  B_FABIA <- as.matrix(B_FABIA[, -zeroes])
}

N <- nrow(X_FABIA)
G <- nrow(B_FABIA)

gene_info[] <- lapply(gene_info, as.character)

er <- sample_info$er
er[er == 0] <- -1

er_order <- order(er)
er_factor <- factor(er)

K_FABIA <- 50
p_vals_fabia <- numeric(K_FABIA)
for (k in 1:K_FABIA) {
  p_vals_fabia[k] <- wilcox.test(X_FABIA[, k] ~ er_factor)$p.value
}

er_cols <- which(p_vals_fabia < 0.01/K_FABIA)
er_cols <- er_cols[order(p_vals_fabia[er_cols], decreasing = F)]

er_plot <- plot_matrix(er[er_order], title = "ER") +
  theme(text = element_text(size= 12))

X_FABIA[X_FABIA > 5] <- 5
X_FABIA[X_FABIA < -5] <- -5


fabia_plot_50 <- plot_matrix(X_FABIA[er_order, union(er_cols, 1:K_FABIA)], 
                         title = "FABIA Factor Matrix (K = 50)") +
  theme(text = element_text(size= 12))

#-------------------------------------------------------------------------

load("fabia_result/fabia_out_10.RData")

K_FABIA <- 10

X_FABIA <- t(out_FABIA@Z)
B_FABIA <- out_FABIA@L

FABIA_biclusters <- extractBic(out_FABIA)
B_FABIA_support <- matrix(0, nrow = G, ncol = K_FABIA)
X_FABIA_support <- matrix(0, nrow = N, ncol = K_FABIA)

for(k in 1:K_FABIA) {
  B_FABIA_support[FABIA_biclusters$numn[k,]$numng, k] <- 1
  X_FABIA_support[FABIA_biclusters$numn[k,]$numnp, k] <- 1
}

X_FABIA[X_FABIA_support == 0] <- 0
B_FABIA[B_FABIA_support == 0] <- 0

zeroes_B <- which(apply(B_FABIA, 2, function(x) all(x == 0)))
zeroes_X <- which(apply(X_FABIA, 2, function(x) all(x == 0)))
zeroes <- union(zeroes_B, zeroes_X)

if (length(zeroes) > 0) {
  X_FABIA <- as.matrix(X_FABIA[, -zeroes])
  B_FABIA <- as.matrix(B_FABIA[, -zeroes])
}


p_vals_fabia <- numeric(K_FABIA)
for (k in 1:K_FABIA) {
  p_vals_fabia[k] <- wilcox.test(X_FABIA[, k] ~ er_factor)$p.value
}

er_plot <- plot_matrix(er[er_order], title = "ER") +
  theme(text = element_text(size= 12))


er_cols <- which(p_vals_fabia < 0.01/K_FABIA)
er_cols <- er_cols[order(p_vals_fabia[er_cols], decreasing = F)]


X_FABIA[X_FABIA > 5] <- 5
X_FABIA[X_FABIA < -5] <- -5


fabia_plot_10 <- plot_matrix(X_FABIA[er_order, union(c(er_cols), 1:K_FABIA)], 
                         title = "FABIA Factor Matrix (K = 10)") +
  theme(text = element_text(size= 12))

pdf(file = "figures/FABIA_output.pdf", width = 8, height = 4)

grid.arrange(er_plot, fabia_plot_10, fabia_plot_50, nrow = 1, widths = c(1.5,11,  20))

dev.off()



