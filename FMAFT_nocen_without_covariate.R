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
