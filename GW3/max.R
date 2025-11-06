# 1.

library(splines)

setwd("C:/Users/Max McCourt/OneDrive/Documents/Master's Year/Semester 1/ESP/Group Works/GW3") ## comment out of submitted
data <- read.table("engcov.txt")



evaluate_x_xtilde_s <- function(data, k = 80, edur = 3.151, sdur = 0.469) {
  
  # create S
  S <- crossprod(diff(diag(k), diff=2))
  
  # create knots
  lower_bound <- min(data$julian) - 30
  upper_bound <- max(data$julian)
  f_seq <- seq(lower_bound, upper_bound, by=1)

  # sequence of k + 4 evenly spaced knots, where the middle k - 2 knots cover lower_bound to upper_bound
  middle_seq_length <- k-2
  middle_seq <- seq(lower_bound, upper_bound, length.out = middle_seq_length)
  
  # Compute step size for handling remaining knots
  step <- middle_seq[2] - middle_seq[1]
  
  # our k+4 evenly spaced knots, where 
  knots <- c(
    seq(to = min(middle_seq) - step, by = step, length.out = 3),  # 3 knots below
    middle_seq,                                                   # middle sequence
    seq(from = max(middle_seq) + step, by = step, length.out = 3) # 3 knots above
  )
  
  # create X_tilde
  X_tilde <- splineDesign(knots, f_seq)
  
  # create probability function for days from infection to death (pd)
  d <- 1:k
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)
  
  n <- length(data$date)
  m <- length(pd)
  r <- nrow(X_tilde)
  # Build weight matrix W such that W[i, t] = pd[30 + i - t] if in range
  W <- matrix(0, n, r)
  for (i in 1:n) {
    # create indexes of [30 + i - 1] to [30 + i - m]
    idxs <- (30 + i - (1:m))
    # subset to only valid indexes (0<=valid<=k)
    valid <- which(idxs > 0 & idxs <= r)
    # take i-th row of W at valid indexes to be pd of valid index
    W[i, idxs[valid]] <- pd[1:length(valid)]
  }
  # Matrix-multiply weight matrix W and X_tilde
  X <- W %*% X_tilde
  
  return(list(X_tilde=X_tilde,S=S,X=X))
}


part1 <- evaluate_x_xtilde_s(data)


# jackson's notes: 
#t <- data$julian
#lags <- 30:0
#sequence <- as.vector(sapply(t, function(x) x - lags)) 
#is more logical way of doing f_seq, but for this method we must remove duplicates (need all unique values t can take)

# max's notes:
# updated function so defining knots is cleaner
# used weight matrix to turn double loop into single loop for finding X (checked its the same dw)
# spelled tilde right lol

# 2.



pnll <- function(data, gamma, X, S, lambda) {
  y <- data$deaths
  B <- exp(gamma)
  # define penalty term
  penalty <- 1/2*lambda*t(B) %*% S %*% B
  # define mu
  mu <- X %*% B
  # define poisson log likelihood (use lfactorial to avoid numeric instability)
  ll_pois <- sum(y*log(mu) - mu - lfactorial(y))
  # define penalised negative log likelihood
  pnll <- -ll_pois + penalty
  return(pnll)
}



pnll_grad <- function(data, gamma, X, S, lambda){
  y <- data$deaths
  B <- exp(gamma)
  # derivative of the penalty term
  penalty_grad <- lambda*diag(B) %*% S %*% B
  # define mu
  mu <- X %*% B
  # derivative of the poisson log likelihood (check further)
  ll_pois_grad <- t(X) %*% (mu - y)
  # define derivative of penalised negative log likelihood
  pnll_grad <- - ll_pois_grad + penalty_grad
  return(pnll_grad)
}



# example:
gamma <- rpois(80, 0.2)
lambda <- 0.5

# testing penalised neg. log likelihood and its gradient
pnll(data, gamma, X, S, lambda)
pnll_grad(data, gamma, X, S, lambda)



# checking finite differencing is not working. need to come back to.
finite_diff <- function(f, gamma, eps = 1e-6) {
  k <- length(gamma)
  grad_approx <- numeric(k)
  
  for (i in 1:k) {
    gamma_eps <- gamma
    gamma_eps[i] <- gamma_eps[i] + eps
    grad_approx[i] <- (f(gamma_eps) - f(gamma)) / eps
  }
  
  return(grad_approx)
}


pnll_wrapper <- function(gamma_vec) {
  pnll(data, gamma_vec, X, S, lambda)
}


gamma0 <- rep(0, ncol(X))  # starting point


grad_numeric <- finite_diff(pnll_wrapper, gamma0)
grad_analytic <- pnll_grad(data, gamma0, X, S, lambda)


max_diff <- max(abs(grad_numeric - grad_analytic))
print(max_diff)
