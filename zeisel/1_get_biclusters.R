# R script to fit biclusters for Zeisel data
library(preprocessCore)

source('../SSLB_functions.R', echo=TRUE)

library(SSLB)
library(fabia)
bicmix_dir <- "../../Code/BicMix/BicMix"

# -------------------------------------
# Read in data
#--------------------------------------

Y_raw <- read.table("data/Y_raw.txt", stringsAsFactors = F)

# normalize quantiles (average quantile normalization)
Y <- normalize.quantiles(t(Y_raw))
Y <- t(Y)

# save normalized gene expression for bicmix
write.table(t(Y), file = "data/Y_ave_bicmix.txt", row.names = F, col.names = F)

Y <- as.matrix(Y)
N <- nrow(Y)
G <- ncol(Y)

K_init <- 100

# -------------------------------------
# Run SSLB
#--------------------------------------

lambda1 <- 1
lambda1_tilde <- 1
lambda0s <- c(1, 5, 10, 50, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000)
lambda0_tildes <- lambda0s


time <- system.time(out <- SSLB(Y, K_init,
                                  lambda0s = lambda0s,
                                  lambda0_tildes = lambda0_tildes,
                                  lambda1 = lambda1,
                                  lambda1_tilde = lambda1_tilde,
                                  alpha = 1/(N),
                                  a = 1/(G * K_init),
                                  a_tilde = 1/(N * K_init),
                                  d = 0,
                                  IBP = 1, EPSILON = 0.05))



save(out, file = "SSLB_result/SSLB_output.RData")


# -------------------------------------
# Run FABIA
#--------------------------------------

time <- system.time(  out_FABIA <- fabia(t(Y), K_FABIA))

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

  save(out_FABIA, file = "fabia_result/fabia_output.RData")



# -------------------------------------
# Run BicMix
#--------------------------------------


system2("mkdir", "bicmix_result")
dir <- "bicmix_result"

bicmix_output <- system2(bicmix_dir,
                           args = c("--nf", as.character(K_init),
                                    "--y", "data/Y_ave_bicmix.txt",
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





