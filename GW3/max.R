                                                                                     #                                                                      #
                                                                                     #                                                                      #
                                                                                    #                                                                        #
                                                                                   #                                                                          #
                                                                                  #                                                                            #
                                                                                  #                                                                            #
                                                                                 #                                                                              #
                                                                                 #                                                                              #
                                                                                 #                                                                              #
                                                                                 #                                                                              #
                                                                                  #                                      #                                     #
                                                                                   #                  #                  #                  #                 #
                                                                                    #                ###                # #                ###               #                                        
                                                                                      ##              #               ##   ##               #             ##                                                
                                                                                        ###                        ###       ###                       ###                                                   
                                                                                           #######         ########             ########        #######                                                                                
                                                                                                   ########                             ########                                                                   
                                                                                                                                                                                                    












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
  # F Matrix
  F <- diag(as.vector(y/mu - 1)) %*% X %*% diag(B)
  # derivative of the poisson log likelihood is apparently colsums of F
  ll_pois_grad <- colSums(F)
  # define derivative of penalised negative log likelihood
  pnll_grad <- - ll_pois_grad + penalty_grad
  return(pnll_grad)
}



# example:
gamma <- rpois(80, 0.2)
lambda <- 0.5

y <- data$deaths
# testing penalised neg. log likelihood and its gradient
pnll(y, gamma, X, S, lambda)
pnll_grad(y, gamma, X, S, lambda)



# Check finite difference

finite_diff <- function(f, gamma, eps = 1e-6) {
  # adds small perturbation to each parameter of vector gamma, and stores the 
  # difference between f(gamma+eps_i) and f(gamma) over the perturbation size,
  # as in first principles of differentiation. The function returns all gradient
  # approximation changes of all parameters of gamma.
  
  # initialise empty numeric vector to store approximated gradient values
  k <- length(gamma)
  grad_approx <- numeric(k)
  
  # for each gamma
  for (i in 1:k) {
    # take copy of gamma vector
    gamma_eps <- gamma
    # add small perturbation epsilon to the i-th parameter of gamma vector
    gamma_eps[i] <- gamma_eps[i] + eps
    # approximate derivative for i-th parameter of gamma using first principles
    grad_approx[i] <- (f(gamma_eps) - f(gamma)) / eps
  }
  
  # return gradient approximations of all dimensions of gamma
  return(grad_approx)
}

# wraps function so we can use it in our finite_diff function
pnll_wrapper <- function(gamma_vec) {
  pnll(y, gamma_vec, X, S, lambda)
}

# initialise gamma of all 0 entries
gamma0 <- rep(0, ncol(X)) 

# finds the gradient approximations of gamma0 (for all dimensions) in pnll function
grad_numeric <- finite_diff(pnll_wrapper, gamma0)
# finds the gradient of pnll using pnll_grad function
grad_analytic <- pnll_grad(y, gamma0, X, S, lambda)

# takes max difference in dimensions of pnll between gradient approximations
# and pnll_grad function
max_diff <- max(abs(grad_numeric - grad_analytic))
print(max_diff)


# Max's notes: 3 functions covered. 
# - pnll(data, gamma, X, S, lambda) finds PNLL (penalised neg. log likelihood) 
# - pnll_grad(data, gamma, X, S, lambda) finds gradient of PNLL 
# - finite_diff creates a function to approximate gradients of a given function
#     (in our case pnll) by introducing small perturbations of gamma in each
#     dimension.
# SUCCESS! The max_diff < 0.001, seems to approximate gradient well.


# 3.

# may come back to



# 4.


# function to define log likelihood (based on pnll function in Q2)
loglik <- function(y, beta_hat, X, S, lambda) {
  # define penalty term
  penalty <- 1/2*lambda*t(beta_hat) %*% S %*% beta_hat
  # define mu
  mu <- X %*% beta_hat
  # define poisson log likelihood (use lfactorial to avoid numeric instability)
  ll_pois <- sum(y*log(mu) - mu - lfactorial(y))
  return(ll_pois)
}



log_lambda_seq <- seq(-13,-7,length=50)
lambda_seq <- 10^(log_lambda_seq)


BIC <- rep(NA, length(lambda_seq))

# iterate through lambda sequence
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
    lambda = lambda,    # "
    control = list(maxit = 1000)
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
  n <- length(y)
  
  # BIC criterion
  BIC[i] <- -2*loglik(y, beta_par, X, S, lambda) + log(n)*EDF
}
 
# minimised BIC criterion
BIC_min <- min(BIC)
# lambda corresponding to minimised BIC criterion
lambda_seq[which(BIC==min(BIC))]


# Max's notes: - Loop to find lambda corresponding to minimised BIC criterion
# I believe loop logic is correct and produces correct results
# However, it is inefficient. It takes ~50s to run. 
# In order to make this more efficient, we will need to find more efficient
# ways to perform matrix multiplication for finding EDF and also look back
# on loglik, pnll and pnll_grad functions and focus on optimisation.
