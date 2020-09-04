# R script for Simulation 3 of SSLB paper

library(fabia)
library(SSLB)
library(mvtnorm)
library(isa2)
library(biclust)

# Requires these functions
source("../SSLB_functions.R")

# directory for C code for BicMix (downloaded from https://www.cs.princeton.edu/~bee/software/BicMix-Code-for-distribution.zip)
bicmix_dir <- as.character(as.matrix(read.table("../bicmix_directory.txt", stringsAsFactors = F)))
paste("../", bicmix_dir, sep = "")

GEN_DATA <- F # if true, generate new data. if false, read in data
RUN_SSLB <- T
RUN_BICMIX <- F
RUN_FABIA <- F
READ_SSBiEM = F
RUN_ISA = F
RUN_SPECTRAL = F
RUN_PLAID = F

set.seed(123456789)

N <- 300 # number of samples
G <- 1000 # number of features
K <- 15 # number of biclusters
K_init <- 30 # number of initial biclusters

nrep <- 50 # number of replications

#-----------------------------------

methods <- c("SSLB_IBP", "SSLB_PY", "SSLB_BB", "Bicmix", "FABIA", "SSBiEM",
             "ISA", "Spectral", "Plaid")
metrics <- c("consensus", "recovery", "relevance", "K")
n_methods <- length(methods)
n_metrics <- length(metrics)
result <- array(NA, dim = c(n_methods, n_metrics, nrep))
rownames(result) <- methods
colnames(result) <- metrics

# system2("mkdir", "data")
# system2("mkdir", "results")

for (r in 1:nrep) {
  
  dir <- paste("data/rep", r, sep = "")
  # system2("mkdir", dir)
  
  # generate data
  if (GEN_DATA) {
    get_data <- generate_sparse_bic_poisson(n_f = N, n_l = G, n_bic = K, min_f = 5, max_f = 20,
                                    min_l = 10, max_l = 50, overlap_f = 5, overlap_l = 15,
                                    mean_f = 2, sd_f = 1, mean_l = 3, sd_l = 1, 
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
    lambda0_tildes <- c(1, rep(5, length(lambda0s)-1))
    
    
    SSLB_output <- capture.output(out <- SSLB(Y, 
                                              K_init,
                                              lambda0s = lambda0s,
                                              lambda0_tildes = lambda0_tildes,
                                              lambda1 = lambda1,
                                              lambda1_tilde = lambda1_tilde, 
                                              alpha = 1, 
                                              d = 0.5,
                                              IBP = 1, EPSILON = 0.01))
    
    X_SSLB <- out$X
    B_SSLB <- out$B
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
    
    X_SSLB <- out$X
    B_SSLB <- out$B
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
  
  
  #---------------------------
  # Bicmix
  #---------------------------
  if (RUN_BICMIX) {
    
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
    
    FABIA_output <- capture.output(out_FABIA <- fabia(t(Y), K), type = "message")
    
    FABIA_biclusters <- extractBic(out_FABIA)
    B_FABIA <- matrix(0, nrow = G, ncol = K)
    X_FABIA <- matrix(0, nrow = N, ncol = K)
    
    for(k in 1:K) {
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
  
  #-------------------------------
  # SSBiEM
  #-------------------------------
  
  if (READ_SSBiEM) {
    X_SSBiEM = as.matrix(read.csv(paste(dir, "/matlab/X_support.txt", sep = ""), header = F))
    B_SSBiEM = as.matrix(read.csv(paste(dir, "/matlab/B_support.txt", sep = ""), header = F))
    
    K_SSBiEM = ncol(X_SSBiEM)
    
    if(K_SSBiEM > 0) {
      result_SSBiEM <- analyze_bic(X_SSBiEM, B_SSBiEM, X, B)
      
      result["SSBiEM", "consensus", r] <- result_SSBiEM$consensus
      result["SSBiEM", "recovery", r] <- result_SSBiEM$recovery
      result["SSBiEM", "relevance", r] <- result_SSBiEM$relevance
      result["SSBiEM", "K", r] <- K_SSBiEM
    } else {
      result["SSBiEM", "K", r] <- K_SSBiEM
    }
  }
  
  
  #-------------------------------
  # ISA
  #-------------------------------
  
  if (RUN_ISA) {
    
    out = isa(Y)
    
    X_ISA = out$rows
    B_ISA = out$columns
    
    K_ISA = ncol(X_ISA)
    
    if(K_ISA > 0) {
      result_ISA <- analyze_bic(X_ISA, B_ISA, X, B)
      
      result["ISA", "consensus", r] <- result_ISA$consensus
      result["ISA", "recovery", r] <- result_ISA$recovery
      result["ISA", "relevance", r] <- result_ISA$relevance
      result["ISA", "K", r] <- K_ISA
    } else {
      result["ISA", "K", r] <- K_ISA
    }
  }
  
  #-------------------------------
  # Spectral
  #-------------------------------
  
  if (RUN_SPECTRAL) {
    
    out = biclust(Y, method = BCSpectral(), normalization = "irrc")
    
    X_spectral = as.matrix(apply(out@RowxNumber, 2, as.numeric))
    B_spectral = as.matrix(apply(out@NumberxCol, 1, as.numeric))
    
    K_spectral = ncol(X_spectral)
    
    if(K_spectral > 0) {
      result_spectral <- analyze_bic(X_spectral, B_spectral, X, B)
      
      result["Spectral", "consensus", r] <- result_spectral$consensus
      result["Spectral", "recovery", r] <- result_spectral$recovery
      result["Spectral", "relevance", r] <- result_spectral$relevance
      result["Spectral", "K", r] <- K_spectral
    } else {
      result["Spectral", "K", r] <- K_spectral
    }
  }
  
  #-------------------------------
  # Plaid
  #-------------------------------
  
  if (RUN_PLAID) {
    
    plaid_out = capture.output(out <- biclust(Y, method = BCPlaid()))
    
    X_plaid = as.matrix(apply(out@RowxNumber, 2, as.numeric))
    B_plaid = as.matrix(apply(out@NumberxCol, 1, as.numeric))
    
    K_plaid = ncol(X_plaid)
    
    if (is.null(K_plaid)) {
      K_plaid = 0
    }
    
    if(K_plaid > 0) {
      result_plaid <- analyze_bic(X_plaid, X_plaid, X, B)
      
      result["Plaid", "consensus", r] <- result_plaid$consensus
      result["Plaid", "recovery", r] <- result_plaid$recovery
      result["Plaid", "relevance", r] <- result_plaid$relevance
      result["Plaid", "K", r] <- K_plaid
    } else {
      result["Plaid", "K", r] <- K_plaid
    }
  }
  
  
  
  print(r)
}


save(result, file = "results/sim_study_3.RData")





