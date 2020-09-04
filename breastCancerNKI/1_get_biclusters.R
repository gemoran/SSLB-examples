# R script to find biclusters using SSLB, FABIA, BicMix, ISA, Plaid and Spectral

source('../SSLB_functions.R', echo=TRUE)

library(SSLB)
library(fabia)
library(isa2)
library(biclust)
library(preprocessCore)

bicmix_dir <- as.character(as.matrix(read.table("../bicmix_directory.txt", stringsAsFactors = F)))
paste("../", bicmix_dir, sep = "")

# -------------------------------------
# Read in data 
#--------------------------------------

option = c("raw", "norm")
i = 1

Y <- read.table(paste("data/Y_", option[i], ".txt", sep = ""), stringsAsFactors = F)
Y <- as.matrix(Y)
N <- nrow(Y)
G <- ncol(Y)

K_init <- 50

# Run SSLB ----------------------------

lambda1 <- 1
lambda1_tilde <- 1
lambda0s <- c(1, 5, 10, 50, 100, 500, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
lambda0_tildes <- lambda0s

time <- system.time(out <- SSLB(Y, K_init,
                                lambda0s = lambda0s,
                                lambda0_tildes = lambda0_tildes,
                                lambda1 = lambda1,
                                lambda1_tilde = lambda1_tilde, 
                                alpha = 1/N, 
                                a = 1/(G * K_init),
                                a_tilde = 1/(N * K_init),
                                d = 0,
                                IBP = 1, EPSILON = 0.01))

save(out, file = paste("SSLB_result/SSLB_out_", option[i], ".RData", sep = ""))

# Run FABIA ----------------------------------------------

# run first with number of biclusters K = 10 
K_FABIA <- 10

time <- system.time(out_FABIA_10 <- fabia(t(Y), K_FABIA))

save(out_FABIA, file = paste("fabia_result/fabia_", option[i], "_10.RData", sep = ""))

# run with number of biclusters from SSLB

K_FABIA = 30

time <- system.time(out_FABIA_K <- fabia(t(Y), K_FABIA))

save(out_FABIA, file = paste("fabia_result/fabia_", option[i], "_K.RData", sep = ""))

# Run Bicmix ---------------------------------------------
dir = "bicmix_result"
system2("mkdir", dir)

# Normalize quantiles for BicMix (Gaussian reference)
Y_std_normal <- normalize.quantiles.use.target(t(Y), target = qnorm(seq(0.000001, 1-0.000001, length.out = G)))
Y_std_normal = t(Y_std_normal)

write.table(t(Y_std_normal), file ="data/Y_bicmix_norm.txt", row.names = FALSE, col.names = FALSE)

for (k in 1:length(option)) {

  dir = paste(dir, "/", option[k], sep = "")
  system2("mkdir", dir)
  
  bicmix_output <- system2(bicmix_dir,
                           args = c("--nf", as.character(K_init),
                                    "--y", paste("data/Y_bicmix_", option[k], ".txt", sep = ""),
                                    "--out", dir, "--sep space"),
                           stdout = T)
}



# ISA ----------------------

out_ISA = isa(Y)

save(out_ISA, file = "ISA_result/out_ISA.RData")

X_ISA = out_ISA$rows
B_ISA = out_ISA$columns

# PLAID ------------------------

plaid_out = capture.output(out <- biclust(Y, method = BCPlaid(), max.layers = 100))

# Spectral --------------------------

out = biclust(exp(Y), method = BCSpectral())

