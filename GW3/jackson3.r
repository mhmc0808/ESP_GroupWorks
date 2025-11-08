
# 1.

library(splines)
library(ggplot2)

# setwd("C:/Users/Max McCourt/OneDrive/Documents/Master's Year/Semester 1/ESP/Group Works/GW3") ## comment out of submitted
data <- read.table("engcov.txt")

# 1

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

X <- part1$X
X_tilde <- part1$X_tilde
S <- part1$S


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



pnll <- function(y, gamma, X, S, lambda) {
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



pnll_grad <- function(y, gamma, X, S, lambda){
  B <- exp(gamma)
  # derivative of the penalty term
  penalty_grad <- lambda*diag(B) %*% S %*% B
  # define mu
  mu <- X %*% B
  # create F Matrix
  F_value <- diag(as.vector(y/mu - 1)) %*% X %*% diag(B) ## !! F is a special character. Dont use it as a variable name
  # sum over all rows to obtain vector of partial derivatives for each gamma_i
  ll_pois_grad <- colSums(F_value)
  # define derivative of penalised negative log likelihood
  pnll_grad <- -ll_pois_grad + penalty_grad
  return(pnll_grad)
}


# example:
gamma <- rpois(80, 0.2)
lambda <- 0.5

y <- data$deaths
# testing penalised neg. log likelihood and its gradient
pnll(y, gamma, X, S, lambda)
pnll_grad(y, gamma, X, S, lambda)


gamma <- rpois(80, 0.2)
lambda <- 0.5
y <- data$deaths
th0 <- gamma  # test paramater vector
fd <- numeric(length(th0)) # approx grad vec
nll0 <- pnll(y = y, gamma = th0, X = X, S = S, lambda = lambda) # pnll at th0
eps <- 1e-7 ## finite difference interval
for (i in 1:length(th0)) { # loop over parameters
  th1 <- th0
  th1[i] <- th1[i] + eps
  nll1 <- pnll(y = y, gamma = th1, X = X, S = S, lambda = lambda)
  fd[i] <- (nll1 - nll0) / eps
}

grad_test <- pnll_grad(y = y, gamma = th0, X = X, S = S, lambda = lambda)
cbind(fd, grad_test) # compare to check pnll_grad is code correctly

max(abs(grad_test-fd)) # find max difference - very small difference, so the code looks good!



# 3.


# fit the model using BFGS
fit <- optim(par = gamma,
             fn = pnll,
             gr = pnll_grad,
             y = y, X = X, S = S, 
             lambda = (5*10^-5), 
             method = "BFGS")

# compute fitted deaths
gamma_hat <- fit$par
# compute B = exp(gamma_hat)
B_hat <- exp(gamma_hat)
# compute fitted deaths (Poisson means)
fitted_deaths <- as.vector(X %*% B_hat)
true_deaths <- data$deaths
day_of_2020 <- data$julian

# combine into data frame
plot_df <- data.frame(
  day = day_of_2020,
  true_deaths = true_deaths,
  fitted_deaths = fitted_deaths)

ggplot(plot_df, aes(x = day)) +
  geom_line(aes(y = true_deaths, color = "Observed"), size = .7) +
  geom_line(aes(y = fitted_deaths, color = "Fitted"), size = .7) +
  labs(x = "Day of 2020",
       y = "Number of Deaths",
       title = "Observed vs Fitted Deaths Over Time",
       color = "Legend") +
  theme_minimal()


# 4.

# function to define log likelihood for use in BIC criterion calculation
ll <- function(y, beta_hat, X, S, lambda) {
  # define mu
  mu <- X %*% beta_hat
  # define poisson log likelihood (use lfactorial to avoid numeric instability)
  ll_pois <- sum(y*log(mu) - mu - lfactorial(y))
  return(ll_pois)
}


log_lambda_seq <- seq(-13,-7,length=50)
lambda_seq <- exp(log_lambda_seq) # in R log() is natural log, so take exp() instead of 10^
n <- length(y)

BIC <- rep(NA, length(lambda_seq))

for (i in 1:length(lambda_seq)){
  # take lambda for this iteration
  lambda <- lambda_seq[i]
  # define optim fit
  fit <- optim(
    par = gamma,        # predicting gamma 
    fn = pnll,          # using pnll as function
    gr = pnll_grad,     # using pnll_grad as gradient of function
    method = "BFGS",
    y=y,                # inputs for pnll and pnll_grad
    X = X,              # "
    S = S,              # "
    lambda = lambda    # "
    # control = list(maxit = 1000) # not sure why you are using this line here !!
  )
  # define beta_par
  beta_par <- exp(fit$par)
  # define mu_hat
  mu_hat <- X %*% beta_par
  # find W, H0 and H_lambda
  W <- diag(as.vector(y/(mu_hat^2)))
  H0 <- t(X) %*% W %*% X
  H_lambda <- H0 + lambda*S
  
  # find effective degrees of freedom of the model
  H_invH0 <- solve(H_lambda, H0) 
  EDF <- sum(diag(H_invH0))
  # n <- length(y) # length(y) doesnt change so define outside of loop !!
  
  # BIC criterion
  BIC[i] <- -2*ll(y, beta_par, X, S, lambda) + log(n)*EDF
}

# minimised BIC criterion
BIC_min <- min(BIC)
# lambda corresponding to minimised BIC criterion
lambda_bic <- lambda_seq[which(BIC==min(BIC))]
lambda_bic

# 5

n_bootstrap <- 200
# intialize empty matrix to store bootstrap replicates
# each row represents a different time point
# each column represents a different bootstrap sample
# each entry represents f(t) at a given time for a given bootstrap
f_hat_boot <- matrix(NA, nrow = nrow(X_tilde), ncol = n_bootstrap)

# update previous pnll function to include weights
pnll_w <- function(y, gamma, X, S, lambda, w) {
  B <- exp(gamma)
  # define penalty term
  penalty <- 1/2*lambda*t(B) %*% S %*% B
  # define mu
  mu <- X %*% B
  # define poisson log likelihood (use lfactorial to avoid numeric instability)
  ll_pois <- sum(w * (y*log(mu) - mu - lfactorial(y))) # weight the log likelihood
  # define penalised negative log likelihood
  pnll <- -ll_pois + penalty
  return(pnll)
}

# update previous pnll_grad function to include weights
pnll_grad_w <- function(y, gamma, X, S, lambda, w){
  B <- exp(gamma)
  # derivative of the penalty term
  penalty_grad <- lambda*diag(B) %*% S %*% B
  # define mu
  mu <- X %*% B
  # create F Matrix
  F_value <- diag(as.vector(w * (y/mu - 1))) %*% X %*% diag(B) # weight the log likelihood portion of the gradient
  # sum over all rows to obtain vector of partial derivatives for each gamma_i
  ll_pois_grad <- colSums(F_value)
  # define derivative of penalised negative log likelihood
  pnll_grad <- -ll_pois_grad + penalty_grad
  return(pnll_grad)
}

# perform the bootstrapping
for (b in 1:n_bootstrap) {
  # generate bootstrap weights 
  wb <- tabulate(sample(n, replace = TRUE), n)
  # fit penalized model with weights
  fit_b <- optim(
    par = gamma, 
    fn = pnll_w, 
    gr = pnll_grad_w,
    y = y, X = X, S = S, lambda = lambda_bic, w = wb,
    method = "BFGS"
  )
  
  # extract fitted coefficients
  gamma_b <- fit_b$par
  B_b <- exp(gamma_b)
  
  # compute f_hat for this bootstrap and add it to f_hat matrix (each column of matrix will represent a f_hat bootstrap sample)
  f_hat_boot[, b] <- X_tilde %*% B_b
}

# 6

f_hat_mean <- rowMeans(f_hat_boot) # obtain average f_hat across replicates for each time period
f_hat_lower <- apply(f_hat_boot, 1, quantile, probs = 0.025) # obtain 2.5% quantity of f_hat_boot
f_hat_upper <- apply(f_hat_boot, 1, quantile, probs = 0.975) # obtain 97.5% quantile of f_hat_boot
# need to define f_seq outside the function from part 1
lower_bound <- min(data$julian) - 30; upper_bound <- max(data$julian); f_seq <- seq(lower_bound, upper_bound, by=1)

# Combine bootstrap results into a data frame
plot_df_6 <- data.frame(
  day = f_seq, # days corresponding to X_tilde
  f_mean = f_hat_mean,
  f_lower = f_hat_lower,
  f_upper = f_hat_upper
)

ggplot() +
  # True deaths
  geom_line(aes(x = day_of_2020, y = true_deaths, color = "Observed Deaths"), size = 0.8) +
  # Fitted deaths
  geom_line(aes(x = day_of_2020, y = fitted_deaths, color = "Fitted Deaths"), size = 0.8) +
  # Daily infection rate with 95% CI
  geom_ribbon(aes(x = f_seq, ymin = f_hat_lower, ymax = f_hat_upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(x = f_seq, y = f_hat_mean, color = "Daily Infections"), size = 0.8) +
  # Secondary axis: rescale infections to match the main axis numerically
  scale_y_continuous(
    name = "Number of Deaths",
    sec.axis = sec_axis(~ ., name = "Daily New Infections")
  ) +
  labs(
    x = "Day of 2020",
    color = "Legend",
    title = "Observed & Fitted Deaths with Daily Infections (95% CI)"
  ) +
  theme_minimal()







