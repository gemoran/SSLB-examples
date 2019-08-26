# R script for Simulation 2 of SSLB paper

library(fabia)
library(SSLB)
library(mvtnorm)

# Requires these functions
source("../SSLB_functions.R")

# directory for C code for BicMix (downloaded from https://www.cs.princeton.edu/~bee/software/BicMix-Code-for-distribution.zip)
# CHANGE TO CORRECT DIRECTORY
bicmix_dir <- "../../Code/BicMix/BicMix"

GEN_DATA <- F # if true, generate new data. if false, read in data
RUN_SSLB <- T
RUN_BICMIX <- T
RUN_FABIA <- T

set.seed(123456789)

N <- 300 # number of samples
G <- 1000 # number of features
K <- 15 # number of biclusters
K_init <- 30 # number of initial biclusters

nrep <- 50 # number of replications

#-----------------------------------

methods <- c("SSLB_IBP", "SSLB_PY", "SSLB_BB", "Bicmix", "FABIA")
metrics <- c("consensus", "recovery", "relevance", "K")
n_methods <- length(methods)
n_metrics <- length(metrics)
result <- array(NA, dim = c(n_methods, n_metrics, nrep))
rownames(result) <- methods
colnames(result) <- metrics


for (r in 1:nrep) {
  
  dir <- paste("data/rep", r, sep = "")
 # system2("mkdir", dir)
  
  # generate data
  if (GEN_DATA) {
    get_data <- generate_dense_bic(n_f = N, n_l = G, n_bic = K, n_dense = 5, min_f = 5, max_f = 20,
                                   min_l = 10, max_l = 50, overlap_f = 5, overlap_l = 15,
                                   mean_f = 2, sd_f = 1, mean_l = 3, sd_l = 1, 
                                   sd_f_dense = 2, sd_l_dense = 2,
                                   sd_f_noise = 0.2, sd_l_noise = 0.2, sd_epsilon = 1)
    
    Y <- get_data$data
    X <- get_data$factors_bic
    B <- get_data$loadings_bic
    
    write.table(t(Y), file = paste(dir, "/Y.txt", sep = ""), row.names = FALSE, col.names = FALSE)
    write.table(X, file = paste(dir, "/X.txt", sep = ""), row.names = FALSE, col.names = FALSE)
    write.table(B, file = paste(dir, "/B.txt", sep = ""), row.names = FALSE, col.names = FALSE)
  } else {
    Y <- as.matrix(read.table(file = paste(dir, "/Y.txt", sep = "")))
    Y <- t(Y)
    X <- as.matrix(read.table(file = paste(dir, "/X.txt", sep = "")))
    B <- as.matrix(read.table(file = paste(dir, "/B.txt", sep = "")))
  }
  
  
  
  #---------------------------
  # SSLB
  #---------------------------
  if (RUN_SSLB) {
    lambda1 <- 1
    lambda1_tilde <- 1
    lambda0s <- c(1, 5, 10, 50, 100, 500, 1000, 10000, 100000, 1000000, 10000000)
    lambda0_tildes <- lambda0s
    
    
    SSLB_output <- capture.output(out <- SSLB(Y, 
                                              K_init,
                                              lambda0s = lambda0s,
                                              lambda0_tildes = lambda0_tildes,
                                              lambda1 = lambda1,
                                              lambda1_tilde = lambda1_tilde, 
                                              alpha = 1, 
                                              d = 0.5,
                                              IBP = 1, EPSILON = 0.01))
    
    X_col <- which(apply(out$X, 2, function(x) sum(x != 0)) < 0.5 * N)
    B_col <- which(apply(out$B, 2, function(x) sum(x != 0)) < 0.5 * G)
    keep <- intersect(X_col, B_col)
    X_SSLB <- out$X[, keep]
    B_SSLB <- out$B[, keep]
    K_SSLB <- ncol(B_SSLB)
    
    if(out$K > 0) {
      result_SSLB <- analyze_bic(X_SSLB, B_SSLB, X, B)
      
      result["SSLB_PY", "consensus", r] <- result_SSLB$consensus
      result["SSLB_PY", "recovery", r] <- result_SSLB$recovery
      result["SSLB_PY", "relevance", r] <- result_SSLB$relevance
      result["SSLB_PY", "K", r] <- K_SSLB
    } else {
      result["SSLB_PY", "K", r] <- K_SSLB
    }
    
    SSLB_output <- capture.output(out <- SSLB(Y, 
                                              K_init,
                                              lambda0s = lambda0s,
                                              lambda0_tildes = lambda0_tildes,
                                              lambda1 = lambda1,
                                              lambda1_tilde = lambda1_tilde, 
                                              alpha = 1, 
                                              d = 0,
                                              IBP = 1, EPSILON = 0.01))
    
    X_SSLB <- out$X
    B_SSLB <- out$B
    
    X_col <- which(apply(X_SSLB, 2, function(x) sum(x != 0)) < 0.5 * N)
    B_col <- which(apply(B_SSLB, 2, function(x) sum(x != 0)) < 0.5 * G)
    keep <- intersect(X_col, B_col)
    X_SSLB <- out$X[, keep]
    B_SSLB <- out$B[, keep]
    K_SSLB <- ncol(B_SSLB)
    
    
    if(out$K > 0) {
      result_SSLB <- analyze_bic(X_SSLB, B_SSLB, X, B)
      
      result["SSLB_IBP", "consensus", r] <- result_SSLB$consensus
      result["SSLB_IBP", "recovery", r] <- result_SSLB$recovery
      result["SSLB_IBP", "relevance", r] <- result_SSLB$relevance
      result["SSLB_IBP", "K", r] <- K_SSLB
    } else {
      result["SSLB_IBP", "K", r] <- K_SSLB
    }
    
    SSLB_output <- capture.output(out <- SSLB(Y, 
                                              K_init,
                                              lambda0s = lambda0s,
                                              lambda0_tildes = lambda0_tildes,
                                              lambda1 = lambda1,
                                              lambda1_tilde = lambda1_tilde, 
                                              alpha = 1, 
                                              d = 0,
                                              IBP = 0, EPSILON = 0.01))
    
    X_col <- which(apply(out$X, 2, function(x) sum(x != 0)) < 0.5 * N)
    B_col <- which(apply(out$B, 2, function(x) sum(x != 0)) < 0.5 * G)
    keep <- intersect(X_col, B_col)
    X_SSLB <- out$X[, keep]
    B_SSLB <- out$B[, keep]
    K_SSLB <- ncol(B_SSLB)
    
    
    if(out$K > 0) {
      result_SSLB <- analyze_bic(X_SSLB, B_SSLB, X, B)
      
      result["SSLB_BB", "consensus", r] <- result_SSLB$consensus
      result["SSLB_BB", "recovery", r] <- result_SSLB$recovery
      result["SSLB_BB", "relevance", r] <- result_SSLB$relevance
      result["SSLB_BB", "K", r] <- K_SSLB
    } else {
      result["SSLB_BB", "K", r] <- K_SSLB
    }
    
  }
  
  if (RUN_BICMIX) {
    #---------------------------
    # Bicmix
    #---------------------------
    
    system2("mkdir", paste(dir, "/result", sep = ""))
    bicmix_output <- system2(bicmix_dir, 
                             args = c("--nf", as.character(K_init),
                                      "--y", paste(dir, "/Y.txt", sep = ""),
                                      "--out", paste(dir, "/result", sep = ""), "--sep space"),
                             stdout = T)
    
    X_bicmix <- try(as.matrix(read.table(paste(dir, "/result/EX", sep = ""), header = F)), silent = TRUE)
    if (inherits(X_bicmix, 'try-error')) {
      X_bicmix <- matrix(0, nrow = K, ncol = N)
    }
    X_bicmix <- t(X_bicmix)
    B_bicmix <- try(as.matrix(read.table(paste(dir, "/result/LAM", sep = ""), header = F)), silent = TRUE)
    if (inherits(B_bicmix, 'try-error')) {
      B_bicmix <- matrix(0, nrow = G, ncol = K)
    }
    Z_bicmix <- try(as.matrix(read.table(paste(dir, "/result/Z", sep = ""), header = F)), silent = TRUE)
    O_bicmix <- try(as.matrix(read.table(paste(dir, "/result/O", sep = ""), header = F)), silent = TRUE)
    keep <- rbind(Z_bicmix[1,], O_bicmix[1,])
    keep <- apply(keep, 2, sum)
    keep <- which(keep == 2)
    
    X_bicmix <- X_bicmix[, keep]
    B_bicmix <- B_bicmix[, keep]
    
    X_bicmix[abs(X_bicmix) < 10e-10] <- 0
    B_bicmix[abs(B_bicmix) < 10e-10] <- 0
    
    zeroes_B <- which(apply(B_bicmix, 2, function(x) all(x == 0)))
    zeroes_X <- which(apply(X_bicmix, 2, function(x) all(x == 0)))
    zeroes <- union(zeroes_B, zeroes_X)
    
    if (length(zeroes) > 0) {
      X_bicmix <- as.matrix(X_bicmix[, -zeroes])
      B_bicmix <- as.matrix(B_bicmix[, -zeroes])
    }
    
    X_norm <- apply(X_bicmix, 2, function(x) sqrt(sum(x^2)))
    B_norm <- apply(B_bicmix, 2, function(x) sqrt(sum(x^2)))
    
    d <- sqrt(X_norm/B_norm)
    X_bicmix <- t(1/d * t(X_bicmix))
    B_bicmix <- t(d * t(B_bicmix))
    
    K_bicmix <- ncol(X_bicmix)
    
    if (K_bicmix > 0) {
      result_bicmix <- analyze_bic(X_bicmix, B_bicmix, X, B)
      
      result["Bicmix", "consensus", r] <- result_bicmix$consensus
      result["Bicmix", "recovery", r] <- result_bicmix$recovery
      result["Bicmix", "relevance", r] <- result_bicmix$relevance
      result["Bicmix", "K", r] <- K_bicmix
    } else {
      result["Bicmix", "K", r] <- K_bicmix
    }
    
  }
  #---------------------------
  # FABIA
  #---------------------------
  
  if (RUN_FABIA) {

    Y_FABIA <- Y
    K_svd <- K
    FABIA_output <- capture.output(out_FABIA <- fabia(t(Y_FABIA), K_svd), type = "message")
    
    FABIA_biclusters <- extractBic(out_FABIA)
    B_FABIA <- matrix(0, nrow = G, ncol = K)
    X_FABIA <- matrix(0, nrow = N, ncol = K)
    
    for(k in 1:K_svd) {
      B_FABIA[FABIA_biclusters$numn[k,]$numng, k] <- 1
      X_FABIA[FABIA_biclusters$numn[k,]$numnp, k] <- 1
    }
    
    zeroes_B <- which(apply(B_FABIA, 2, function(x) all(x == 0)))
    zeroes_X <- which(apply(X_FABIA, 2, function(x) all(x == 0)))
    zeroes <- union(zeroes_B, zeroes_X)
    
    if (length(zeroes) > 0) {
      X_FABIA <- as.matrix(X_FABIA[, -zeroes])
      B_FABIA <- as.matrix(B_FABIA[, -zeroes])
    }
    
    X_col <- which(apply(X_FABIA, 2, function(x) sum(x != 0)) < 0.5 * N)
    B_col <- which(apply(B_FABIA, 2, function(x) sum(x != 0)) < 0.5 * G)
    keep <- intersect(X_col, B_col)
    X_FABIA <- X_FABIA[, keep]
    B_FABIA <- B_FABIA[, keep]
    K_FABIA <- ncol(B_FABIA)
    

    
    K_FABIA <- ncol(X_FABIA)
    
    if (K_FABIA > 0) {
      result_FABIA <- analyze_bic(X_FABIA, B_FABIA, X, B)
      
      result["FABIA", "consensus", r] <- result_FABIA$consensus
      result["FABIA", "recovery", r] <- result_FABIA$recovery
      result["FABIA", "relevance", r] <- result_FABIA$relevance
      result["FABIA", "K", r] <- K_FABIA
      
    } else {
      result["FABIA", "K", r] <- K_FABIA
    }
    
  }

  print(r)
}



save(result, file = "results/sim_study_2.RData")
