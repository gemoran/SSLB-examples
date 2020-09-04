library(fabia)
library(SSLB)
library(mvtnorm)
library(isa2)
library(biclust)

# Requires these functions
source("../SSLB_functions.R")

set.seed(123456789)

N <- 300 # number of samples
G <- 1000 # number of features
K <- 15 # number of biclusters
K_init <- 30 # number of initial biclusters

r = 1


dir <- paste("data/rep", r, sep = "")
# system2("mkdir", dir)

# READ DATA
Y <- as.matrix(read.table(file = paste(dir, "/Y.txt", sep = "")))
Y <- t(Y)
X <- as.matrix(read.table(file = paste(dir, "/X.txt", sep = "")))
B <- as.matrix(read.table(file = paste(dir, "/B.txt", sep = "")))

# pdf("figures/X_true.pdf", width = 2, height = 4)
# 
# plot_matrix(X, title = "True X (K = 15)", xlab = "K", ylab = "N")
# 
# dev.off()
# 
# pdf("figures/B_true.pdf", width = 2, height = 4)
# 
# plot_matrix(B, title = "True B (K = 15)", xlab = "K", ylab = "N")
# 
# dev.off()

## SSLB 

lambda1 <- 1
lambda1_tilde <- 1
lambda0s <- c(1, 5, 10, 50, 100, 500, 1000, 10000, 100000, 1000000, 10000000)
lambda0_tildes <- c(1, rep(5, length(lambda0s)-1))

SSLB_output <- capture.output(out <- SSLB(Y, 
                                          K_init,
                                          lambda0s = lambda0s,
                                          lambda0_tildes = lambda0_tildes,
                                          lambda1 = lambda1,
                                          lambda1_tilde = lambda1_tilde, 
                                          alpha = 1, 
                                          d = 0,
                                          IBP = 0))


X_SSLB <- out$X
B_SSLB <- out$B
K_SSLB <- ncol(B_SSLB)

result_SSLB = analyze_bic(X_SSLB, B_SSLB, X, B)

X_SSLB = result_SSLB$X_found
B_SSLB = result_SSLB$B_found

X_SSLB = X_SSLB[, union(c(12, 1:9, 11), 1:K_SSLB)]
B_SSLB = B_SSLB[, union(c(12, 1:9, 11), 1:K_SSLB)]


pdf("figures/X_SSLB.pdf", width = 2, height = 4)

plot_matrix(X_SSLB, title = paste("X (SSLB, K = ", K_SSLB, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()


pdf("figures/B_SSLB.pdf", width = 2, height = 4)

plot_matrix(B_SSLB, title = paste("B (SSLB, K = ", K_SSLB, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

#--------------------------------------

X_bicmix <- try(as.matrix(read.table(paste(dir, "/result/EX", sep = ""), header = F)), silent = TRUE)
if (inherits(X_bicmix, 'try-error')) {
  X_bicmix <- matrix(0, nrow = K, ncol = N)
}
X_bicmix <- t(X_bicmix)
B_bicmix <- try(as.matrix(read.table(paste(dir, "/result/LAM", sep = ""), header = F)), silent = TRUE)
if (inherits(B_bicmix, 'try-error')) {
  B_bicmix <- matrix(0, nrow = G, ncol = K)
}

X_bicmix[abs(X_bicmix) < 10e-10] <- 0
B_bicmix[abs(B_bicmix) < 10e-10] <- 0

result_bicmix = analyze_bic(X_bicmix, B_bicmix, X, B)

X_bicmix[X_bicmix!=0]=1
B_bicmix[B_bicmix!=0]=1


K_bicmix = ncol(X_bicmix)

pdf("figures/X_bicmix.pdf", width = 2, height = 4)

plot_matrix(X_bicmix, title = paste("X (BicMix, K = ", K_bicmix, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

pdf("figures/B_bicmix.pdf", width = 2, height = 4)

plot_matrix(B_bicmix, title = paste("B (BicMix, K = ", K_bicmix, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

#-----------------------------------

FABIA_output <- capture.output(out_FABIA <- fabia(t(Y), K), type = "message")

X_FABIA <- t(out_FABIA@Z)
B_FABIA <- out_FABIA@L

K_FABIA = ncol(X_FABIA)

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

K_FABIA = ncol(X_FABIA)


result_fabia = analyze_bic(X_FABIA, B_FABIA, X, B)
X_FABIA = result_fabia$X_found
B_FABIA = result_fabia$B_found


pdf("figures/X_fabia.pdf", width = 2, height = 4)

plot_matrix(X_FABIA, title = paste("X (FABIA, K = ", K_FABIA, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

pdf("figures/B_fabia.pdf", width = 2, height = 4)

plot_matrix(B_FABIA, title = paste("B (FABIA, K = ", K_FABIA, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

#-----------


X_SSBiEM = as.matrix(read.csv(paste(dir, "/matlab/X.txt", sep = ""), header = F))
B_SSBiEM = as.matrix(read.csv(paste(dir, "/matlab/B.txt", sep = ""), header = F))
X_SSBiEM_support = as.matrix(read.csv(paste(dir, "/matlab/X_support.txt", sep = ""), header = F))
B_SSBiEM_support = as.matrix(read.csv(paste(dir, "/matlab/B_support.txt", sep = ""), header = F))

X_SSBiEM[X_SSBiEM_support == 0] = 0
B_SSBiEM[B_SSBiEM_support == 0] = 0

K_SSBiEM = ncol(X_SSBiEM)

result_ssbiem = analyze_bic(X_SSBiEM, B_SSBiEM, X, B)
X_SSBiEM = result_ssbiem$X_found
B_SSBiEM = result_ssbiem$B_found


pdf("figures/X_ssbiem_large.pdf", width = 2, height = 4)

plot_matrix(X_SSBiEM, title = paste("X (SSBiEM, K = ", K_SSBiEM, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

pdf("figures/B_ssbiem_large.pdf", width = 2, height = 4)

plot_matrix(B_SSBiEM, title = paste("B (SSBiEM, K = ", K_SSBiEM, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

#-------------------------------------------------

out = isa(Y)

X_ISA = out$rows
B_ISA = out$columns

K_ISA = ncol(X_ISA)

result_ISA = analyze_bic(X_ISA, B_ISA, X, B)
X_ISA = result_ISA$X_found
B_ISA = result_ISA$B_found

pdf("figures/X_ISA.pdf", width = 2, height = 4)

plot_matrix(X_ISA, title = paste("X (ISA, K = ", K_ISA, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

pdf("figures/B_ISA.pdf", width = 2, height = 4)

plot_matrix(B_ISA, title = paste("B (ISA, K = ", K_ISA, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

#--------------------


out = biclust(exp(Y), method = BCSpectral())

X_spectral = as.matrix(apply(out@RowxNumber, 2, as.numeric))
B_spectral = as.matrix(apply(out@NumberxCol, 1, as.numeric))

K_spectral = ncol(X_spectral)

result_spectral = analyze_bic(X_spectral, B_spectral, X, B)
X_spectral = result_spectral$X_found
B_spectral = result_spectral$B_found

pdf("figures/X_spectral.pdf", width = 2, height = 4)

plot_matrix(X_spectral, title = paste("X (Spectral, K = ", K_spectral, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

pdf("figures/B_spectral.pdf", width = 2, height = 4)

plot_matrix(B_spectral, title = paste("B (Spectral, K = ", K_spectral, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()
#-----------------


plaid_out = capture.output(out <- biclust(Y, method = BCPlaid(), max.layers = K))

X_plaid = as.matrix(apply(out@RowxNumber, 2, as.numeric))
B_plaid = as.matrix(apply(out@NumberxCol, 1, as.numeric))

K_plaid = ncol(X_plaid)

if (is.null(K_plaid)) {
  K_plaid = 0
}

result_plaid = analyze_bic(X_plaid, B_plaid, X, B)
X_plaid = result_plaid$X_found
B_plaid = result_plaid$B_found

pdf("figures/X_plaid.pdf", width = 2, height = 4)

plot_matrix(X_plaid, title = paste("X (Plaid, K = ", K_plaid, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

pdf("figures/B_plaid.pdf", width = 2, height = 4)

plot_matrix(B_plaid, title = paste("B (Plaid, K = ", K_plaid, ")", sep = ""), xlab = "K", ylab = "N")

dev.off()

