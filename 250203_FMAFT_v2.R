###############################
# 1. EM algorithm function without covariate

em_mixture_aft_without_covariate <- function(logT, K = 2, max_iter = 200, tol = 1e-3) {
  n <- length(logT)
  
  # Initialize parameters
  sigma <- c(2, 5)         
  pi <- rep(1/K, K)        
  beta0 <- 2              
  gamma <- matrix(0, nrow = n, ncol = K)   
  
  # Design matrix
  X_design <- matrix(1, n, 1)
  
  for (iter in 1:max_iter) {
    
    ## E-step:
    for (k in 1:K) {
      gamma[, k] <- pi[k] * dnorm(logT, mean = beta0, sd = sigma[k])
    }
    gamma <- gamma / rowSums(gamma)
  
    pi_old    <- pi
    sigma_old <- sigma
    beta0_old <- beta0
    
    ## M-step:
    # (1) 
    pi <- colMeans(gamma)
    
    # (2) Intercept beta0 업데이트 (Weighted Least Squares)
    # 각 관측치에 대한 가중치: w_i = sum_{k=1}^K gamma[i,k] / sigma[k]^2
    sigma_sq <- sigma^2
    w <- rowSums(sweep(gamma, 2, sigma_sq, FUN = "/"))
    # 행렬 X가 1벡터니까까, weighted LS 해는:
    beta0 <- sum(w * logT) / sum(w)
    
    # (3) Update sigma
    residuals <- logT - beta0
    for (k in 1:K) {
      sigma[k] <- sqrt(sum(gamma[, k] * residuals^2) / sum(gamma[, k]))
    }
    
    # Convergence
    if (max(abs(pi - pi_old),
            abs(sigma - sigma_old),
            abs(beta0 - beta0_old)) < tol) {
      cat("Converged at iteration:", iter, "\n")
      return(list(pi = pi, sigma = sigma, beta0 = beta0, iter = iter))
    }
  }
  cat("Not converged in", max_iter, "iterations.\n")
  return(list(pi = pi, sigma = sigma, beta0 = beta0, iter = max_iter))
}

###############################
# Simulation & Result 
###############################

n_sim <- 100
pi_results    <- matrix(0, nrow = n_sim, ncol = 2)
sigma_results <- matrix(0, nrow = n_sim, ncol = 2)
beta0_results <- numeric(n_sim)
conv_iter     <- numeric(n_sim)

# 시작 시간 기록
start <- Sys.time()

for (i in 1:n_sim) {
  set.seed(i)
  n <- 200
  
  # True parameters
  pi_true    <- c(0.6, 0.4)
  sigma_true <- c(1, 2)
  beta0_true <- 1
  
  # 데이터 생성 (공변량 없음)
  # 각 관측치에 대해 성분 할당
  component <- sample(1:2, n, replace = TRUE, prob = pi_true)
  # 오차: 각 관측치에 대해 해당 성분의 sigma 적용
  epsilon <- rnorm(n, mean = 0, sd = sigma_true[component])
  
  # logT 생성: 모형에 따라 (beta0_true + epsilon)
  logT <- beta0_true + epsilon
  
  # EM algorithm
  result <- em_mixture_aft_without_covariate(logT, K = 2, max_iter = 200, tol = 1e-3)
  
  # Reslt
  pi_results[i, ]    <- result$pi
  sigma_results[i, ] <- result$sigma
  beta0_results[i]   <- result$beta0
  conv_iter[i]       <- result$iter
}

end <- Sys.time()
cat("Simulation time:", end - start, "\n")

# 결과정리리
summary_df <- data.frame(
  Parameter = c("pi1", "pi2", "sigma1", "sigma2", "beta0", "Convergence Iterations"),
  Estimate  = c(mean(pi_results[, 1]), mean(pi_results[, 2]),
                mean(sigma_results[, 1]), mean(sigma_results[, 2]),
                mean(beta0_results), mean(conv_iter)),
  True      = c(pi_true[1], pi_true[2], sigma_true[1], sigma_true[2], beta0_true, NA)
)
print(summary_df)

###############################
# Boxplots for beta0 and sigma estimates
###############################

# Boxplot for beta0 estimates
boxplot(beta0_results, 
        main = "Boxplot of beta0 Estimates", 
        ylab = "beta0", 
        col = "lightblue")

# Boxplot for sigma estimates (sigma1 and sigma2)
sigma_df <- data.frame(sigma1 = sigma_results[,1], sigma2 = sigma_results[,2])
boxplot(sigma_df, 
        main = "Boxplot of Sigma Estimates", 
        ylab = "Sigma Value", 
        names = c("sigma1", "sigma2"),
        col = c("lightcoral", "lightgoldenrod"))



#-------------------------------------------------------------#

# 2. EM algorithm function with covariate(em_mixture_aft_with_covariates)

em_mixture_aft_with_covariate <- function(logT, X, K = 2, max_iter = 200, tol = 1e-3) {
  n <- length(logT)
  
  # Initialize
  sigma <- c(2, 5)        
  pi <- rep(1/K, K)       
  beta0 <- 2               
  beta1 <- 2               
  gamma <- matrix(0, nrow = n, ncol = K)   
  
  # X : Design matrix
  X_design <- cbind(1, X)
  
  for (iter in 1:max_iter) {
    
    ## E-step: conditional probability
    for (k in 1:K) {
      gamma[, k] <- pi[k] * dnorm(logT, mean = beta0 + beta1 * X, sd = sigma[k])
    }
    gamma <- gamma / rowSums(gamma) 
    
    
    pi_old    <- pi
    sigma_old <- sigma
    beta0_old <- beta0
    beta1_old <- beta1
    
    ## M-step:
    # (1) Update pi
    pi <- colMeans(gamma)
    
    # (2) Update beta (Weighted LSE)
    
    # 최적화 문제:  minimize sum_i [ (sum_k gamma_ik/sigma_k^2) * (logT_i - (beta0+beta1 X_i))^2 ]
    # 각 관측치 i의 가중치는: w_i = sum_{k=1}^K gamma[i,k] / sigma[k]^2
    sigma_sq <- sigma^2
    
    # sweep()를 이용해서서 각 열별로 sigma_sq를 나누고, 행별 합 구함함
    w <- rowSums(sweep(gamma, 2, sigma_sq, FUN = "/"))
    W <- diag(w)  # Weight (diagonal) matrix
    beta_update <- solve(t(X_design) %*% W %*% X_design, t(X_design) %*% W %*% logT)
    beta0 <- beta_update[1]
    beta1 <- beta_update[2]
    
    # (3) Update sigma
    residuals <- logT - (beta0 + beta1 * X)
    for (k in 1:K) {
      sigma[k] <- sqrt(sum(gamma[, k] * residuals^2) / sum(gamma[, k]))
    }
    
    # Check Convergence
    if (max(abs(pi - pi_old),
            abs(sigma - sigma_old),
            abs(beta0 - beta0_old),
            abs(beta1 - beta1_old)) < tol) {
      cat("Converged at iteration:", iter, "\n")
      return(list(pi = pi, sigma = sigma, beta_0 = beta0, beta_1 = beta1, iter = iter))
    }
  }
  cat("Not converged in", max_iter, "iterations.\n")
  return(list(pi = pi, sigma = sigma, beta_0 = beta0, beta_1 = beta1, iter = max_iter))
}

###############################

# Result vector 
n_sim <- 100
pi_results    <- matrix(0, nrow = n_sim, ncol = 2)
sigma_results <- matrix(0, nrow = n_sim, ncol = 2)
beta0_results <- numeric(n_sim)
beta1_results <- numeric(n_sim)
conv_iter     <- numeric(n_sim)  

# 시작 시간 기록
start <- Sys.time()

# Simulation iteration
for (i in 1:n_sim) {
  set.seed(i)
  n <- 1000
  
  # True parameter
  pi_true    <- c(0.6, 0.4)
  sigma_true <- c(1, 2)
  beta0_true <- 1
  beta1_true <- 1
  
  # X Covariate
  X <- rbinom(n, 1, prob = 0.5)
  
  # Component (k= 1 or 2)
  component <- sample(1:2, n, replace = TRUE, prob = pi_true)
  
  # epsilon ~ N(0, sigma_k^2)
  epsilon <- rnorm(n, mean = 0, sd = sigma_true[component])
  
  # log(T)
  logT <- beta0_true + beta1_true * X + epsilon
  
  # EM algorithm
  result <- em_mixture_aft_with_covariate(logT, X, K = 2, max_iter = 200, tol = 1e-3)
  
  # Result
  pi_results[i, ]    <- result$pi
  sigma_results[i, ] <- result$sigma
  beta0_results[i]   <- result$beta_0
  beta1_results[i]   <- result$beta_1
  conv_iter[i]       <- result$iter
}

end <- Sys.time()

# Simulation time
cat("Simulation time:", end - start, "\n")


# Summary
summary_df <- data.frame(
  Parameter = c("pi1", "pi2", "sigma1", "sigma2", "beta0", "beta1", "Convergence Iterations"),
  Estimate  = c(mean(pi_results[, 1]), mean(pi_results[, 2]), 
                mean(sigma_results[, 1]), mean(sigma_results[, 2]), 
                mean(beta0_results), mean(beta1_results), mean(conv_iter)),
  True      = c(pi_true[1], pi_true[2], sigma_true[1], sigma_true[2], beta0_true, beta1_true, ".")
)

print(summary_df)

###############################
# Boxplot for Beta and Sigma estimates
###############################

# Boxplot for beta estimates (beta0 and beta1)
beta_df <- data.frame(beta0 = beta0_results, beta1 = beta1_results)
boxplot(beta_df, 
        main = "Boxplot of Beta Estimates", 
        ylab = "Beta Value",
        names = c("beta0", "beta1"),
        col = c("lightblue", "lightgreen"))

# Boxplot for sigma estimates (sigma1 and sigma2)
sigma_df <- data.frame(sigma1 = sigma_results[,1], sigma2 = sigma_results[,2])
boxplot(sigma_df, 
        main = "Boxplot of Sigma Estimates", 
        ylab = "Sigma Value", 
        names = c("sigma1", "sigma2"),
        col = c("lightcoral", "lightgoldenrod"))

#############################################################################

# 3. EM algorithm function with covariates(em_mixture_aft_with_covariates)

em_mixture_aft_with_covariates <- function(logT, X, K = 2, max_iter = 200, tol = 1e-3) {
  # X: matrix of covariates (each column one covariate)
  n <- length(logT)
  
  # Initialize parameters
  sigma <- c(2, 5)           # 초기 표준편차 (각 성분)
  pi <- rep(1/K, K)          # 혼합비 초기값
  # 회귀계수: intercept + (beta1, beta2)
  beta <- rep(2, ncol(cbind(1, X)))  # 여기서는 길이 3 (beta0, beta1, beta2)
  gamma <- matrix(0, nrow = n, ncol = K)  # 책임도 행렬
  
  # Design matrix: intercept + covariates X1, X2
  X_design <- cbind(1, X)
  
  for (iter in 1:max_iter) {
    
    ## E-step: 
    for (k in 1:K) {
      gamma[, k] <- pi[k] * dnorm(logT, mean = as.vector(X_design %*% beta), sd = sigma[k])
    }
    gamma <- gamma / rowSums(gamma)  
    
    pi_old <- pi
    sigma_old <- sigma
    beta_old <- beta
    
    ## M-step:
    # (1) Update pi
    pi <- colMeans(gamma)
    
    # (2) Update beta (Weighted Least Squares)
    # 각 관측치 i에 대한 가중치: w_i = sum_{k=1}^K gamma[i,k] / sigma[k]^2
    sigma_sq <- sigma^2
    w <- rowSums(sweep(gamma, 2, sigma_sq, FUN = "/"))
    W <- diag(w)  # 대각 가중치 행렬
    beta <- solve(t(X_design) %*% W %*% X_design, t(X_design) %*% W %*% logT)
    
    # (3) Update sigma
    residuals <- logT - as.vector(X_design %*% beta)
    for (k in 1:K) {
      sigma[k] <- sqrt(sum(gamma[, k] * residuals^2) / sum(gamma[, k]))
    }
    
    # Check convergence
    if (max(abs(pi - pi_old),
            abs(sigma - sigma_old),
            abs(beta - beta_old)) < tol) {
      cat("Converged at iteration:", iter, "\n")
      return(list(pi = pi, sigma = sigma, beta = beta, iter = iter))
    }
  }
  cat("Not converged in", max_iter, "iterations.\n")
  return(list(pi = pi, sigma = sigma, beta = beta, iter = max_iter))
}

##################

# Result vector
n_sim <- 100
pi_results    <- matrix(0, nrow = n_sim, ncol = 2)
sigma_results <- matrix(0, nrow = n_sim, ncol = 2)
beta0_results <- numeric(n_sim)
beta1_results <- numeric(n_sim)
beta2_results <- numeric(n_sim)
conv_iter     <- numeric(n_sim)

# Simulation time
start <- Sys.time()

for (i in 1:n_sim) {
  set.seed(i)
  n <- 200
  
  # True parameters
  pi_true    <- c(0.6, 0.4)
  sigma_true <- c(1, 2)
  beta0_true <- 1
  beta1_true <- 0
  beta2_true <- 1
  
  
  X1 <- rbinom(n, 1, prob = 0.5)
  X2 <- rnorm(n, 0, 1)
  
  X <- cbind(X1, X2)
  
  # Component
  component <- sample(1:2, n, replace = TRUE, prob = pi_true)
  
  # Epsilon
  epsilon <- rnorm(n, mean = 0, sd = sigma_true[component])
  
  # log(T)
  logT <- beta0_true + beta1_true * X1 + beta2_true * X2 + epsilon
  
  # EM algorithm
  result <- em_mixture_aft_with_covariates(logT, X, K = 2, max_iter = 200, tol = 1e-3)
  
  # Result
  pi_results[i, ]    <- result$pi
  sigma_results[i, ] <- result$sigma
  beta0_results[i]   <- result$beta[1]
  beta1_results[i]   <- result$beta[2]
  beta2_results[i]   <- result$beta[3]
  conv_iter[i]       <- result$iter
}

# 종료 시간 기록
end <- Sys.time()
cat("Simulation time:", end - start, "\n")

# 평균 추정 결과와 실제 값 비교
summary_df <- data.frame(
  Parameter = c("pi1", "pi2", "sigma1", "sigma2", "beta0", "beta1", "beta2", "Convergence Iterations"),
  Estimate  = c(mean(pi_results[, 1]), mean(pi_results[, 2]),
                mean(sigma_results[, 1]), mean(sigma_results[, 2]),
                mean(beta0_results), mean(beta1_results), mean(beta2_results),
                mean(conv_iter)),
  True      = c(pi_true[1], pi_true[2], sigma_true[1], sigma_true[2],
                beta0_true, beta1_true, beta2_true, NA)
)

print(summary_df)

#-----------------------------------------------------------------------##############################

# Boxplot for beta estimates (beta0, beta1, beta2)
beta_df <- data.frame(beta0 = beta0_results, beta1 = beta1_results, beta2 = beta2_results)
boxplot(beta_df, 
        main = "Boxplot of Beta Estimates", 
        ylab = "Beta Value", 
        names = c("beta0", "beta1", "beta2"),
        col = c("lightblue", "lightgreen", "lightpink"))

# Boxplot for sigma estimates (sigma1 and sigma2)
sigma_df <- data.frame(sigma1 = sigma_results[,1], sigma2 = sigma_results[,2])
boxplot(sigma_df, 
        main = "Boxplot of Sigma Estimates", 
        ylab = "Sigma Value", 
        names = c("sigma1", "sigma2"),
        col = c("lightcoral", "lightgoldenrod"))
