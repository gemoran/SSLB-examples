# R script to find biclusters using SSLB, FABIA, BicMix

source('../SSLB_functions.R', echo=TRUE)

library(SSLB)
library(fabia)
bicmix_dir <- "../../Code/BicMix/BicMix"

# -------------------------------------
# Read in data
#--------------------------------------

Y <- read.table("data/Y_raw.txt", stringsAsFactors = F)
Y <- as.matrix(Y)
N <- nrow(Y)
G <- ncol(Y)

K_init <- 50

# -------------------------------------
# Run SSLB
#--------------------------------------


lambda1 <- 1
lambda1_tilde <- 1
lambda0s <- c(1, 5, 10, 50, 100, 500, 1000, 10000, 100000, 1000000, 10000000, 100000000)
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
                                  IBP = 1, EPSILON = 0.05))
  

save(out, file = "SSLB_result/SSLB_out.RData")


# -------------------------------------
# Run FABIA
#--------------------------------------

# run first with number of biclusters K = 50 
K_FABIA <- 50

time <- system.time(out_FABIA_50 <- fabia(t(Y), K_FABIA))

save(out_FABIA_50, file = "fabia_result/fabia_out_50.RData")

# run with number of biclusters K = 10 

K_FABIA <- 10

time <- system.time(out_FABIA_10 <- fabia(t(Y), K_FABIA))


save(out_FABIA_10, file = "fabia_result/fabia_out_10.RData")


# -------------------------------------
# Run BicMix
#--------------------------------------


system2("mkdir", "bicmix_result")

  dir <- "bicmix_result"
  bicmix_output <- system2(bicmix_dir,
                           args = c("--nf", as.character(K_init),
                                    "--y", "data/Y_raw_bicmix.txt",
                                    "--out", dir, "--sep space"),
                           stdout = T)
  
  
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




