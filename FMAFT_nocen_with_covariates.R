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
