# Jackson Cramer (s2274544), Maximillian McCourt (s2145762), Natalia Montalvo Cabornero (s1969053)


library(splines)

data <- read.table("engcov.txt")


# 1.

evaluate_x_xtilda_s <- function(data, k = 80, edur = 3.151, sdur = 0.469) {
  # create S
  S <- crossprod(diff(diag(k), diff=2))
  # create knots
  f_seq <- seq(min(data$julian) - 30, max(data$julian), by=1)
  lower_bound <- f_seq[1]
  upper_bound <- f_seq[length(f_seq)]
  # sequence of k + 4 evenly spaced knots, where the middle k-2 knots cover lower_bound to upper_bound
  num_knots <- k + 4
  middle_knots <- k - 2
  middle_seq <- seq(lower_bound, upper_bound, length.out = middle_knots)
  # Compute step size for handling remaining knots
  step <- middle_seq[2] - middle_seq[1]
  # Combine middle sequence with 2 knots before and 2 after
  knots <- c(
    middle_seq[1] - 3*step,
    middle_seq[1] - 2*step, 
    middle_seq[1] - step, 
    middle_seq, 
    middle_seq[middle_knots] + step, 
    middle_seq[middle_knots] + 2*step,
    middle_seq[middle_knots] + 3*step
  )
  # create X_tilda
  X_tilda <- splineDesign(knots, f_seq)
  # create X
  d <- 1:k
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)
  n <- length(data$date)
  X <- matrix(0, n, k)
  for (i in 1:n) {
    Xi <- rep(0, k)
    for (j in 1:min(29+i,k)) {
      idx <- 30+i-j
      current_x_tilda <- X_tilda[idx,]
      current_pd <- pd[j]
      Xi <- Xi + current_x_tilda*current_pd
      }
    X[i, ] <- Xi
  }
  return(list(X_tilda=X_tilda,S=S,X=X))
}


part1 <- evaluate_x_xtilda_s(data)

# note: 
#t <- data$julian
#lags <- 30:0
#sequence <- as.vector(sapply(t, function(x) x - lags)) 
#is more logical way of doing f_seq, but for this method we must remove duplicates (need all unique values t can take)

# 2.

penalised_neg_log_likelihood <- function(gamma, X, S, lambda) {
  y <- data$deaths
  B <- exp(gamma)
  penalty <- 1/2*lambda*t(B) %*% S %*% B
  mu <- X %*% B
  log_likelihood_pois <- sum(y*log(mu) - mu - log(factorial(y)))
  penalized_neg_ll = -log_likelihood_pois + penalty
  return(penalized_neg_ll)
}








