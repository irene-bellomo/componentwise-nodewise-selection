#Implementation of Componentwise Nodewise Selection
# DGP A, nodes are connected up to distance 2
library(fda)
library(mvtnorm)
library(glmnet)
library(igraph)
library(pROC)

##################
# Initialization #
##################
set.seed(2025)
p <- 50
func.path <- ""
#The functions for the two DGPs are taken from 'High-dimensionalfunctional graphical model structure learning via neighborhood selection'
#approach by Boxin Zhao, Percy S. Zhai, Y. Samuel Wang, and Mladen Kolar.
save.path <- ""



# Global Parameter Settings
mu <- 15 # number of basis used to generate data
M <- 5 # number of basis used for neighborhood selection
n <- 100
tau <- t<-100 # number of observations
L <- 100 # number of lambdas
K <- 5 # number of folds of SCV
t.vec <- c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0) # vector of t, where threshold epsilon = t * lambda
len.t <- length(t.vec)
l_basis<-6

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
source(paste(func.path,"A.prelim.func.R", sep="/")) # For Model A generation
source(paste(func.path,"D.prelim.func.R", sep="/")) # For Model D generation
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))   # FGLasso



####################################
#     PART 1: DATA GENERATION      #
####################################

#    Generate Random Functions     
#      and Observation h_ijk       

# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# 0. Generate precision matrix and real adjacency matrix
set.seed(2025)
Theta <- cov.mat.model.A(p, mu) # p*mu by p*mu large square matrix
G.true <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(Theta[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
      G.true[i,j] <- 1
  }
}
set.seed(2025)
R <- 50
adj_list <- vector("list", R)
score_list <- vector("list", R) 

for(r in seq_len(R)){
  cat("Monte Carlo simulation:", r, "\n")
  delta <- rmvnorm(n, sigma = solve(Theta))
  
  obs.time <- seq(1/tau, 1, length.out = tau)
  b.mat.list <- lapply(seq_len(p), function(j) fda.fourier.mat(obs.time, mu))
  
  h <- array(0, c(n, p, tau))
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      idx <- ((j-1)*mu + 1):(j*mu)
      h[i,j,] <- b.mat.list[[j]] %*% delta[i, idx] + rnorm(tau, 0, 0.5)
    }
  }
  
  basis <- create.bspline.basis(c(0,1), nbasis = l_basis)
  coef_ar <- array(NA, dim = c(n, p, l_basis))
  
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      fdj <- smooth.basis(seq(0,1,length.out=tau), h[i,j,], basis)$fd
      coef_ar[i,j,] <- fdj$coefs
    }
  }
  
  coef_s <- coef_ar
  for(j in seq_len(p)){
    for(k in seq_len(l_basis)){
      coef_s[,j,k] <- scale(coef_s[,j,k])
    }
  }
  
  adj_est   <- matrix(0, p, p)
  score_mat <- matrix(0, p, p) 
  
  for(j in seq_len(p)){
    Yj <- coef_s[,j,]  
    Others <- setdiff(seq_len(p), j)
    Xj     <- do.call(cbind, lapply(Others, function(k) coef_s[,k,]))
    
    fit   <- cv.glmnet(Xj, Yj, family="mgaussian", alpha=1)
    coefs <- coef(fit, s="lambda.1se")
    
    idx0 <- 1
    for(k in Others){
      beta_mat <- sapply(coefs, function(mm) mm[idx0:(idx0+5),1])
      score    <- sqrt(sum(beta_mat^2))
      score_mat[j,k] <- score
      if(score > 0) adj_est[j,k] <- 1
      idx0 <- idx0 + l_basis
    }
  }
  
  adj_list[[r]]   <- (adj_est & t(adj_est)) * 1
  score_list[[r]] <- (score_mat + t(score_mat)) / 2
}

thresholds <- seq(0, max(unlist(score_list)), length.out = 100)
roc_curve <- function(score_list, G, thresholds){
  upper <- upper.tri(G)
  TPR <- FPR <- numeric(length(thresholds))
  for(i in seq_along(thresholds)){
    th <- thresholds[i]
    TP <- FP <- FN <- TN <- 0
    for(s in seq_along(score_list)){
      est <- (score_list[[s]] >= th)*1  
      est <- (est & t(est))*1
      TP   <- TP + sum(est[upper]==1 & G[upper]==1)
      FP   <- FP + sum(est[upper]==1 & G[upper]==0)
      FN   <- FN + sum(est[upper]==0 & G[upper]==1)
      TN   <- TN + sum(est[upper]==0 & G[upper]==0)
    }
    TPR[i] <- TP/(TP+FN)
    FPR[i] <- FP/(FP+TN)
  }
  list(FPR=FPR, TPR=TPR)
}


roc_my <- roc_curve(score_list, G.true, thresholds)
mc.results<-"" # directory where the files with the Monte Carlo simulation results 
# for the FGLasso and Neighborhood Selection methods are located for dgp A (from  Boxin Zhao et al.)
setwd(mc.results) 

fpr_mat <- matrix(NA, nrow=100, ncol=50)
tpr_mat <- matrix(NA, nrow=100, ncol=50)
# Initialize matrices

fpr_fglasso <- tpr_fglasso <- matrix(NA, nrow=100, ncol=50)
fpr_gX <- tpr_gX <- matrix(NA, nrow=100, ncol=50)
fpr_gY <- tpr_gY <- matrix(NA, nrow=100, ncol=50)

# Load each file and save the ROC results of the three methods

for(i in 1:50){
  fname <- paste0("ROC.A.050.RunInd", i, ".RData")
  load(fname)
  
  fpr_fglasso[, i] <- roc[, "FPR.fglasso"]
  tpr_fglasso[, i] <- roc[, "TPR.fglasso"]
  
  fpr_gX[, i] <- roc[, "FPR.and.gX"]
  tpr_gX[, i] <- roc[, "TPR.and.gX"]
  
  fpr_gY[, i] <- roc[, "FPR.and.gY"]
  tpr_gY[, i] <- roc[, "TPR.and.gY"]
}

fpr_fglasso_mean <- rowMeans(fpr_fglasso)
tpr_fglasso_mean <- rowMeans(tpr_fglasso)

fpr_gX_mean <- rowMeans(fpr_gX)
tpr_gX_mean <- rowMeans(tpr_gX)

fpr_gY_mean <- rowMeans(fpr_gY)
tpr_gY_mean <- rowMeans(tpr_gY)
plot(roc_my$FPR, roc_my$TPR, type="l", col="purple", lwd=3,
     xlab="False Positive Rate", ylab="True Positive Rate", xlim=c(0,1), ylim=c(0,1))

lines(fpr_fglasso_mean, tpr_fglasso_mean, col="#E69F00", lwd=3)
lines(fpr_gX_mean, tpr_gX_mean, col="#333333", lwd=3)

legend("bottomright",
       legend = c("CNS", "FGLasso", "FPCA - Zhao"),
       col = c("purple", "#E69F00", "#333333"),
       lty = c(1, 1, 1),
       lwd = 2)

auc_trap <- function(fpr, tpr){
  ord <- order(fpr)
  x   <- fpr[ord]
  y   <- tpr[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
auc_mio <- auc_trap(roc_my$FPR, roc_my$TPR)
auc_fglasso <- auc_trap(fpr_fglasso_mean, tpr_fglasso_mean)
auc_gX <- auc_trap(fpr_gX_mean, tpr_gX_mean)

cat("AUC – CNS:",     round(auc_mio, 3), "\n")
cat("AUC – FGLasso:",         round(auc_fglasso, 3), "\n")
cat("AUC – Zhao:",   round(auc_gX, 3), "\n")


get_metrics <- function(est, G) {
  upper_idx <- upper.tri(G)
  
  TP <- sum(est[upper_idx] == 1 & G[upper_idx] == 1)
  FP <- sum(est[upper_idx] == 1 & G[upper_idx] == 0)
  FN <- sum(est[upper_idx] == 0 & G[upper_idx] == 1)
  TN <- sum(est[upper_idx] == 0 & G[upper_idx] == 0)
  
  TDR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)  # True Discovery Rate
  FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)  # False Discovery Rate
  ACC <- (TP + TN) / (TP + FP + FN + TN)            # Accuracy
  
  return(c(TDR = TDR, FDR = FDR, ACC = ACC))
}

results <- t(sapply(adj_list, get_metrics, G = G.true))

round(apply(results, 2, mean, na.rm = TRUE), 3)



#########################################
#Elastic net (alpha=0.5)
set.seed(2025)
R <- 50
adj_list_e <- vector("list", R)
score_list_e <- vector("list", R) 

for(r in seq_len(R)){
  cat("Monte Carlo simulation:", r, "\n")
  delta <- rmvnorm(n, sigma = solve(Theta))
  
  obs.time <- seq(1/tau, 1, length.out = tau)
  b.mat.list <- lapply(seq_len(p), function(j) fda.fourier.mat(obs.time, mu))
  
  h <- array(0, c(n, p, tau))
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      idx <- ((j-1)*mu + 1):(j*mu)
      h[i,j,] <- b.mat.list[[j]] %*% delta[i, idx] + rnorm(tau, 0, 0.5)
    }
  }
  
  basis <- create.bspline.basis(c(0,1), nbasis = l_basis)
  coef_ar <- array(NA, dim = c(n, p, l_basis))
  
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      fdj <- smooth.basis(seq(0,1,length.out=tau), h[i,j,], basis)$fd
      coef_ar[i,j,] <- fdj$coefs
    }
  }
  
  coef_s <- coef_ar
  for(j in seq_len(p)){
    for(k in seq_len(l_basis)){
      coef_s[,j,k] <- scale(coef_s[,j,k])
    }
  }
  
  adj_est   <- matrix(0, p, p)
  score_mat <- matrix(0, p, p) 
  
  for(j in seq_len(p)){
    Yj <- coef_s[,j,]  
    Others <- setdiff(seq_len(p), j)
    Xj     <- do.call(cbind, lapply(Others, function(k) coef_s[,k,]))
    
    fit   <- cv.glmnet(Xj, Yj, family="mgaussian", alpha=0.5)
    coefs <- coef(fit, s="lambda.1se")
    
    idx0 <- 1
    for(k in Others){
      beta_mat <- sapply(coefs, function(mm) mm[idx0:(idx0+5),1])
      score    <- sqrt(sum(beta_mat^2))
      score_mat[j,k] <- score
      if(score > 0) adj_est[j,k] <- 1
      idx0 <- idx0 + l_basis
    }
  }
  
  adj_list_e[[r]]   <- (adj_est & t(adj_est)) * 1
  score_list_e[[r]] <- (score_mat + t(score_mat)) / 2
}

thresholds <- seq(0, max(unlist(score_list_e)), length.out = 100)

roc_my_e <- roc_curve(score_list_e, G.true, thresholds)
plot(roc_my$FPR, roc_my$TPR,
     type="l", col="purple", lwd=3,
     xlab="False Positive Rate", ylab="True Positive Rate",
     xlim=c(0,1), ylim=c(0,1))

lines(roc_my_e$FPR, roc_my_e$TPR, type="l", col="darkgreen", lwd=3)
lines(fpr_fglasso_mean, tpr_fglasso_mean, col="#E69F00", lwd=3)
lines(fpr_gX_mean, tpr_gX_mean, col="#333333", lwd=3)

legend("bottomright",
       legend = c("CNS - Lasso", "FGLasso", "FPCA - Zhao", "CNS - Elastic Net"),
       col = c("purple", "#E69F00", "#333333", "darkgreen"),
       lty = c(1, 1, 1,1),
       lwd = 2)

auc_trap <- function(fpr, tpr){
  ord <- order(fpr)
  x   <- fpr[ord]
  y   <- tpr[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
auc_mio <- auc_trap(roc_my$FPR, roc_my$TPR)
auc_mio_e <- auc_trap(roc_my_e$FPR, roc_my_e$TPR)
auc_fglasso <- auc_trap(fpr_fglasso_mean, tpr_fglasso_mean)
auc_gX <- auc_trap(fpr_gX_mean, tpr_gX_mean)

cat("AUC – CNS Lasso:",     round(auc_mio, 3), "\n")
cat("AUC – CNS Elastic Net:",     round(auc_mio_e, 3), "\n")
cat("AUC – FGLasso:",         round(auc_fglasso, 3), "\n")
cat("AUC – Zhao:",   round(auc_gX, 3), "\n")


get_metrics <- function(est, G) {
  upper_idx <- upper.tri(G)
  
  TP <- sum(est[upper_idx] == 1 & G[upper_idx] == 1)
  FP <- sum(est[upper_idx] == 1 & G[upper_idx] == 0)
  FN <- sum(est[upper_idx] == 0 & G[upper_idx] == 1)
  TN <- sum(est[upper_idx] == 0 & G[upper_idx] == 0)
  
  TDR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)  # True Discovery Rate
  FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)  # False Discovery Rate
  ACC <- (TP + TN) / (TP + FP + FN + TN)            # Accuracy
  
  return(c(TDR = TDR, FDR = FDR, ACC = ACC))
}

results <- t(sapply(adj_list_e, get_metrics, G = G.true))

round(apply(results, 2, mean, na.rm = TRUE), 3)





#########################################
# Adalasso (with ridge coefficients)
set.seed(2025)
R         <- 50
l_basis   <- 6   # number of basis spline
adj_list_a   <- vector("list", R)
score_list_a <- vector("list", R)

for(r in seq_len(R)){
  cat("Monte Carlo simulation:", r, "\n")
  
  delta <-rmvnorm(n, sigma = solve(Theta))
  
  obs.time  <- seq(1/tau, 1, length.out = tau)
  b.mat.list <- lapply(seq_len(p), function(j) fda.fourier.mat(obs.time, mu))
  h <- array(0, c(n, p, tau))
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      idx <- ((j-1)*mu + 1):(j*mu)
      h[i,j,] <- b.mat.list[[j]] %*% delta[i, idx] + rnorm(tau, 0, 0.5)
    }
  }
  
  basis   <- create.bspline.basis(c(0,1), nbasis = l_basis)
  coef_ar <- array(NA, dim = c(n, p, l_basis))
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      fdj <- smooth.basis(seq(0,1,length.out = tau), h[i,j,], basis)$fd
      coef_ar[i,j,] <- fdj$coefs
    }
  }
  coef_s <- coef_ar
  for(j in seq_len(p)){
    for(k in seq_len(l_basis)){
      coef_s[,j,k] <- scale(coef_s[,j,k])
    }
  }
  
  score_mat_a <- matrix(0, p, p)
  adj_est_a   <- matrix(0, p, p)
  
  for(j in seq_len(p)){
    Yj     <- coef_s[, j, ]       
    Others <- setdiff(seq_len(p), j)
    Xj     <- do.call(cbind, lapply(Others, function(k) coef_s[, k, ]))
    
    fit_ridge    <- cv.glmnet(Xj, Yj, family="mgaussian", alpha=0)
    ridge_coefs  <- coef(fit_ridge, s = "lambda.min")
    
    penalty_factors <- numeric((p-1) * l_basis)
    idx_pf <- 1
    for(k in Others){
      beta_k_ridge <- sapply(ridge_coefs,
                             function(mat) mat[idx_pf:(idx_pf + l_basis - 1), 1])
      group_norms  <- sqrt(colSums(beta_k_ridge^2))
      penalty_factors[idx_pf:(idx_pf + l_basis - 1)] <- 1 / (group_norms + 1e-4)
      idx_pf <- idx_pf + l_basis
    }
    
    fit_adalasso <- cv.glmnet(Xj, Yj,
                              family="mgaussian",
                              alpha=1,
                              penalty.factor=penalty_factors)
    coef_list <- coef(fit_adalasso, s = "lambda.1se")
    
    idx0 <- 1
    for(k in Others){
      beta_mat <- sapply(coef_list,
                         function(mm) mm[idx0:(idx0 + l_basis - 1), 1])
      score    <- sqrt(sum(beta_mat^2))
      score_mat_a[j, k] <- score
      if(score > 0) adj_est_a[j, k] <- 1
      idx0 <- idx0 + l_basis
    }
  }
  
  adj_list_a[[r]]   <- (adj_est_a & t(adj_est_a)) * 1
  score_list_a[[r]] <- (score_mat_a + t(score_mat_a)) / 2
}


thresholds <- seq(0, max(unlist(score_list_a)), length.out = 100)
roc_my_a   <- roc_curve(score_list_a, G.true, thresholds)


plot(roc_my_a$FPR, roc_my_a$TPR,
     type="l", col="#E69F00", lwd=3,
     xlab="False Positive Rate", ylab="True Positive Rate",
     xlim=c(0,1), ylim=c(0,1))

lines(roc_my_e$FPR, roc_my_e$TPR, type="l", col="darkgreen", lwd=3)
lines(roc_my$FPR, roc_my_a$TPR, type="l", col="blue", lwd=3, lty=2)
lines(fpr_fglasso_mean, tpr_fglasso_mean, col="purple", lwd=3)
lines(fpr_gX_mean, tpr_gX_mean, col="#333333", lwd=3, lty=2)

legend("bottomright",
       legend = c("CNS - Lasso", "FGLasso", "NS - Zhao", "CNS - Elastic Net","CNS - AdaLasso"),
       col = c("blue", "purple", "#333333", "darkgreen", "#E69F00"),
       lty = c(2,1, 2, 1,1),
       lwd = 2)

auc_mio <- auc_trap(roc_my$FPR, roc_my$TPR)
auc_mio_e <- auc_trap(roc_my_e$FPR, roc_my_e$TPR)
auc_mio_a <- auc_trap(roc_my_a$FPR, roc_my_e$TPR)
auc_fglasso <- auc_trap(fpr_fglasso_mean, tpr_fglasso_mean)
auc_gX <- auc_trap(fpr_gX_mean, tpr_gX_mean)

cat("AUC – CNS Lasso:",     round(auc_mio, 3), "\n")
cat("AUC – CNS Elastic Net:",     round(auc_mio_e, 3), "\n")
cat("AUC – CNS AdaLasso:",     round(auc_mio_a, 3), "\n")
cat("AUC – FGLasso:",         round(auc_fglasso, 3), "\n")
cat("AUC – Zhao:",   round(auc_gX, 3), "\n")


results <- t(sapply(adj_list_e, get_metrics, G = G.true))

round(apply(results, 2, mean, na.rm = TRUE), 3)

results <- t(sapply(adj_list_a, get_metrics, G = G.true))
round(apply(results, 2, mean, na.rm = TRUE), 3)

setwd(save.path)
save(adj_list, adj_list_e, adj_list_a,score_list,score_list_e,score_list_a, roc_my, roc_my_e,roc_my_a, file = "roc dgp a.RData")



####################################################
#Implementation of Componentwise Nodewise Selection
# DGP B, nodes are randomly connected 
set.seed(2025)

##################
# Initialization #
##################
set.seed(2025)
p <- 50
func.path <- ""
#The functions for the two DGPs are taken from 'High-dimensionalfunctional graphical model structure learning via neighborhood selection'
#approach by Boxin Zhao, Percy S. Zhai, Y. Samuel Wang, and Mladen Kolar.
save.path <- ""


# Global Parameter Settings
mu <- 15 # number of basis used to generate data
M <- 5 # number of basis used for neighborhood selection
n <- 100
tau <- t<-100 # number of observations
L <- 100 # number of lambdas
K <- 5 # number of folds of SCV
t.vec <- c(0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0) # vector of t, where threshold epsilon = t * lambda
len.t <- length(t.vec)
l_basis<-6

# Read R files
source(paste(func.path,"ADMM.new.R", sep="/"))      # ADMM Function
source(paste(func.path,"prec.rec.R", sep="/"))      # Precision & Recall
source(paste(func.path,"auc.R", sep="/"))           # Computing AUC from ROC
source(paste(func.path,"FPCA.score.R", sep="/"))    # Computing FPCA scores
source(paste(func.path,"bases.func.R", sep="/"))    # Basis functions for FPCA
source(paste(func.path,"A.prelim.func.R", sep="/")) # For Model A generation
source(paste(func.path,"D.prelim.func.R", sep="/")) # For Model D generation
source(paste(func.path,"ProxAlg_FGM.R", sep="/"))   # FGLasso

####################################
#     PART 1: DATA GENERATION      #
####################################

#    Generate Random Functions     
#      and Observation h_ijk       

# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# 0. Generate precision matrix and real adjacency matrix
set.seed(2025)
Theta <- cov.mat.model.D(p, mu) # p*mu by p*mu large square matrix
G.true <- matrix(0, p, p) # p by p adjacency matrix
for(i in 1:p){
  for(j in 1:p){
    if(sum(abs(Theta[((i-1)*mu+1):(i*mu), ((j-1)*mu+1):(j*mu)])) > 0)
      G.true[i,j] <- 1
  }
}
set.seed(2025)
R <- 50
adj_list_b <- vector("list", R)
score_list_b <- vector("list", R) 

for(r in seq_len(R)){
  cat("Monte Carlo simulation:", r, "\n")
  delta <- rmvnorm(n, sigma = solve(Theta))
  
  obs.time <- seq(1/tau, 1, length.out = tau)
  b.mat.list <- lapply(seq_len(p), function(j) fda.fourier.mat(obs.time, mu))
  
  h <- array(0, c(n, p, tau))
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      idx <- ((j-1)*mu + 1):(j*mu)
      h[i,j,] <- b.mat.list[[j]] %*% delta[i, idx] + rnorm(tau, 0, 0.5)
    }
  }
  
  basis <- create.bspline.basis(c(0,1), nbasis = l_basis)
  coef_ar <- array(NA, dim = c(n, p, l_basis))
  
  for(i in seq_len(n)){
    for(j in seq_len(p)){
      fdj <- smooth.basis(seq(0,1,length.out=tau), h[i,j,], basis)$fd
      coef_ar[i,j,] <- fdj$coefs
    }
  }
  
  coef_s <- coef_ar
  for(j in seq_len(p)){
    for(k in seq_len(l_basis)){
      coef_s[,j,k] <- scale(coef_s[,j,k])
    }
  }
  
  adj_est   <- matrix(0, p, p)
  score_mat <- matrix(0, p, p) 
  
  for(j in seq_len(p)){
    Yj <- coef_s[,j,]  
    Others <- setdiff(seq_len(p), j)
    Xj     <- do.call(cbind, lapply(Others, function(k) coef_s[,k,]))
    
    fit   <- cv.glmnet(Xj, Yj, family="mgaussian", alpha=1)
    coefs <- coef(fit, s="lambda.1se")
    
    idx0 <- 1
    for(k in Others){
      beta_mat <- sapply(coefs, function(mm) mm[idx0:(idx0+5),1])
      score    <- sqrt(sum(beta_mat^2))
      score_mat[j,k] <- score
      if(score > 0) adj_est[j,k] <- 1
      idx0 <- idx0 + l_basis
    }
  }
  
  adj_list_b[[r]]   <- (adj_est & t(adj_est)) * 1
  score_list_b[[r]] <- (score_mat + t(score_mat)) / 2
}

thresholds <- seq(0, max(unlist(score_list)), length.out = 100)
roc_curve <- function(score_list, G, thresholds){
  upper <- upper.tri(G)
  TPR <- FPR <- numeric(length(thresholds))
  for(i in seq_along(thresholds)){
    th <- thresholds[i]
    TP <- FP <- FN <- TN <- 0
    for(s in seq_along(score_list)){
      est <- (score_list[[s]] >= th)*1  
      est <- (est & t(est))*1
      TP   <- TP + sum(est[upper]==1 & G[upper]==1)
      FP   <- FP + sum(est[upper]==1 & G[upper]==0)
      FN   <- FN + sum(est[upper]==0 & G[upper]==1)
      TN   <- TN + sum(est[upper]==0 & G[upper]==0)
    }
    TPR[i] <- TP/(TP+FN)
    FPR[i] <- FP/(FP+TN)
  }
  list(FPR=FPR, TPR=TPR)
}


roc_my_b <- roc_curve(score_list_b, G.true, thresholds)
setwd(mc.results) # directory where the files with the Monte Carlo simulation results 
# for the FGLasso and Neighborhood Selection methods are located
fpr_mat <- matrix(NA, nrow=100, ncol=50)
tpr_mat <- matrix(NA, nrow=100, ncol=50)

fpr_fglasso <- tpr_fglasso <- matrix(NA, nrow=100, ncol=50)
fpr_gX <- tpr_gX <- matrix(NA, nrow=100, ncol=50)
fpr_gY <- tpr_gY <- matrix(NA, nrow=100, ncol=50)

for(i in 1:50){
  fname <- paste0("ROC.D.050.RunInd", i, ".RData")
  load(fname) 
  
  fpr_fglasso[, i] <- roc[, "FPR.fglasso"]
  tpr_fglasso[, i] <- roc[, "TPR.fglasso"]
  
  fpr_gX[, i] <- roc[, "FPR.and.gX"]
  tpr_gX[, i] <- roc[, "TPR.and.gX"]
  
  fpr_gY[, i] <- roc[, "FPR.and.gY"]
  tpr_gY[, i] <- roc[, "TPR.and.gY"]
}

fpr_fglasso_mean <- rowMeans(fpr_fglasso)
tpr_fglasso_mean <- rowMeans(tpr_fglasso)

fpr_gX_mean <- rowMeans(fpr_gX)
tpr_gX_mean <- rowMeans(tpr_gX)

fpr_gY_mean <- rowMeans(fpr_gY)
tpr_gY_mean <- rowMeans(tpr_gY)
plot(roc_my_b$FPR, roc_my$TPR, type="l", col="blue", lwd=3,
     xlab="False Positive Rate", ylab="True Positive Rate", xlim=c(0,1), ylim=c(0,1))

lines(fpr_fglasso_mean, tpr_fglasso_mean, col="purple", lwd=3)
lines(fpr_gX_mean, tpr_gX_mean, col="#333333", lwd=3)

legend("bottomright",
       legend = c("CNS - Lasso", "FGLasso", "NS - Zhao"),
       col = c("blue", "purple", "#333333"),
       lty = c(1, 1, 1),
       lwd = 2)

auc_trap <- function(fpr, tpr){
  ord <- order(fpr)
  x   <- fpr[ord]
  y   <- tpr[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
auc_mio_b <- auc_trap(roc_my_b$FPR, roc_my_b$TPR)
auc_fglasso <- auc_trap(fpr_fglasso_mean, tpr_fglasso_mean)
auc_gX <- auc_trap(fpr_gX_mean, tpr_gX_mean)

cat("AUC – CNS:",     round(auc_mio_b, 3), "\n")
cat("AUC – FGLasso:",         round(auc_fglasso, 3), "\n")
cat("AUC – Zhao:",   round(auc_gX, 3), "\n")


get_metrics <- function(est, G) {
  upper_idx <- upper.tri(G)
  
  TP <- sum(est[upper_idx] == 1 & G[upper_idx] == 1)
  FP <- sum(est[upper_idx] == 1 & G[upper_idx] == 0)
  FN <- sum(est[upper_idx] == 0 & G[upper_idx] == 1)
  TN <- sum(est[upper_idx] == 0 & G[upper_idx] == 0)
  
  TDR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)  # True Discovery Rate
  FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)  # False Discovery Rate
  ACC <- (TP + TN) / (TP + FP + FN + TN)            # Accuracy
  
  return(c(TDR = TDR, FDR = FDR, ACC = ACC))
}

results <- t(sapply(adj_list_b, get_metrics, G = G.true))

round(apply(results, 2, mean, na.rm = TRUE), 3)


setwd(save.path)
save(adj_list, adj_list_b, score_list,score_list_b, roc_my,roc_my_b, file = "roc dgp b.RData")

     