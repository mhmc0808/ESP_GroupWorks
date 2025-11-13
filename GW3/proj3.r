# Jackson Cramer (s2274544), Maximillian McCourt (s2145762), Natalia Montalvo Cabornero (s1969053)

# Our group worked collaboratively on the development and optimisation of all exercises.
# Although tasks were divided among us, we maintained a balanced workload through 
# regular code reviews, debugging sessions, and problem-solving meetings.
# Jackson focused ---
# Max focused on ---
# Natalia focused on ---

# Link to github repo: https://github.com/mhmc0808/ESP_GroupWorks.git
# Note about github repo: Our group used the same repo for the other projects.
# The code for project 3 can be found in the GW3 folder, and is titled proj3.r

######### GENERAL DESCRIPTION #########

# This project implements a smooth deconvolution to estimate the daily COVID-19 
# infections from observed death data. The model reconstructs the infection 
# trajectory by deconvolving death counts with a delay distribution from infection
# to death. 

#setwd("C:/Users/Max McCourt/OneDrive/Documents/Master's Year/Semester 1/ESP/Group Works/GW3")

# Import necessary libraries and read data
library(splines)
library(ggplot2)
data <- read.table("engcov.txt")

######### MODEL SETUP ##########

# evaluate_matrices() function constructs the following matrices and components 
# of the deconvolution model; 
#   X_tilde - B-spline basis matrix, 
#   X - Convolution matrix mapping infections to deaths via delay distribution, 
#   S - penalty matrix for smoothness regulatisation.
evaluate_matrices <- function(n, k = 80, edur = 3.151, sdur = 0.469, f_seq) {
  
  # Inputs:
  #   n - number of observed death days 
  #   k - number of B-spline basis functions
  #   edur - log(mean) time from infection to death
  #   sdur - log(standard deviation) time from infection to death
  #   f_seq - sequence of infection days
  
  # Penalty matrix for second-order differences
  S <- crossprod(diff(diag(k), diff=2))
  
  # Construct sequence of k+4 knots for B-splines
  # First, construct k-2 equally-spaced knots from lower to upper bound of f_seq
  middle_seq <- seq(f_seq[1], f_seq[length(f_seq)], length.out = k-2)
  # Then define step size of our knots
  step <- middle_seq[2] - middle_seq[1]
  # Use step size to extend sequence of knots
  knots <- c(
    seq(to = min(middle_seq) - step, by = step, length.out = 3), # 3 knots below
    middle_seq,                                                  # middle knots
    seq(from = max(middle_seq) + step, by = step, length.out = 3)# 3 knots above
  )
  
  # B-spline basis matrix for infection function f(t)
  X_tilde <- splineDesign(knots, f_seq)
  
  # Delay distribution from infection to death (pd)
  d <- 1:k
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)      # normalise to sum to 1
  print(length(pd))
  m <- length(pd)         # number of delay days  
  r <- nrow(X_tilde)      # number of infection days
  
  # Build convolution weight matrix W
  # W[i, t] = pd[30 + i - t] if in valid range, else 0
  W <- matrix(0, n, r)
  for (i in 1:n) {
    idxs <- (30 + i - (1:m))              # indices accounting for 30-day lag
    valid <- which(idxs > 0 & idxs <= r)  # restrict indices to valid range
    # Assign pd to valid infection days for i-th row
    W[i, idxs[valid]] <- pd[1:length(valid)] 
  }
  
  # Convolution matrix mapping infections to expected deaths
  X <- W %*% X_tilde
  
  return(list(X_tilde=X_tilde,S=S,X=X))
}

## --- Generate model matrices --- ##

# Define number of observed death days in our dataset
n <- length(data$date) 

# Lower and upper bounds of our days of infection
lb <- min(data$julian) - 30 # initialise 30 days before the first observed death
ub <- max(data$julian)
# Sequence of days of infection from defined lower to upper bound
f_seq <- seq(lb, ub, by=1)

# Construct model matrices 
model_matrices <- evaluate_matrices(n, f_seq=f_seq)

# Extract individual matrices
X <- model_matrices$X
X_tilde <- model_matrices$X_tilde
S <- model_matrices$S


######### PENALISED NEGATIVE LOG-LIKELIHOOD FUNCTIONS ##########

## --- Penalised Negative Log-Likelihood Function --- ##

# pnll() function finds the penalised negative log-likelihood for fitting a 
# Poisson deconvolution model, adding a smoothness penalty to prevent 
# overfitting. We transform gamma vector to exponential scale (beta) to enforce 
# positive coefficients for the B-spline basis functions.
pnll <- function(y, gamma, X, S, lambda) {
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter
  
  B <- exp(gamma)               # transformation of coefficients
  mu <- X %*% B                 # expected deaths
  
  # Smoothing penalty term
  penalty <- 1/2*lambda * sum(B * (S %*% B))
  
  # Poisson log-likelihood 
  # (note log(y!) value dropped as it will not affect beta optimisation)
  ll_pois <- sum(y*log(mu) - mu)
  
  # Penalised negative log likelihood
  pnll <- -ll_pois + penalty
  return(pnll)
}


## --- Gradient of the Penalised Negative Log-Likelihood Function --- ##

# pnll_grad() function computes the gradient of the penalised negative 
# log-likelihood defined in our previous function pnll().
pnll_grad <- function(y, gamma, X, S, lambda){
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter
  
  B <- exp(gamma)                     # transform to original scale
  mu <- X %*% B                       # expected deaths
  
  # Smoothness penalty gradient
  penalty_grad <- lambda * (S %*% B) * B  
  
  # Gradient of Poisson log-likelihood
  F_value <- diag(as.vector(y/mu - 1)) %*% X %*% diag(B)            
  ll_pois_grad <- colSums(F_value)
  
  # Gradient of penalised negative log-likelihood
  pnll_grad <- - ll_pois_grad + penalty_grad
  
  return(pnll_grad)
}


## --- Test functions using finite differences --- ##

# Extract data for y
y <- data$nhs
# Define arbitrary lambda and k
lambda <- 5e-5
k <- 80
# Initialise gamma as zero-vector
gamma <- rep(0,k)

# Initialise parameter vector as zero-vector
th0 <- gamma         
# Initialise vector to store finite-difference gradient estimates
fd <- numeric(length(th0))     

# Evaluate penalised negative log-likelihood at th0
pnll0 <- pnll(y = y, gamma = th0, X = X, S = S, lambda = lambda) 

# Small perturbation step for finite difference calculations
eps <- 1e-7                   

# Loop over parameters of th0
for (i in 1:length(th0)){ 
  # Define th1 as th0 but with perturbation of eps for i-th parameter 
  th1 <- th0; th1[i] <- th1[i] + eps
  # Compute pnll for perturbed th1
  pnll1 <- pnll(y = y, gamma = th1, X = X, S = S, lambda = lambda)
  # Finite-difference approximation with respect to parameter i 
  fd[i] <- (pnll1 - pnll0) / eps
}

# Find analytic gradient at th0 
grad_test <- pnll_grad(y = y, gamma = th0, X = X, S = S, lambda = lambda)

# Maximum absolute difference (across all parameters) between analytic gradient 
# and finite difference gradient approximation 
max_diff <- max(abs(grad_test-fd)) 
print(paste("Maximum gradient difference:", round(max_diff,8)))


######### INITIAL MODEL FITTING ##########

# Fit model using BFGS method with initial lambda for sanity check.
# This provides a preliminary fit using pnll() and pnll_grad() function to 
# verify the model is working correctly before proceeding with smoothing 
# parameter selection (using pre-defined y, X, S variables).
fit <- optim(par = gamma, fn = pnll, gr = pnll_grad,
             y = y, X = X, S = S, lambda = lambda, method = "BFGS", 
             control=list(maxit=1000, reltol=1e-8)) # ensure convergence

# Extract estimated parameters after optimisation
gamma_hat <- fit$par
B_hat <- exp(gamma_hat) # transform to original beta scale

# Compute fitted deaths and fitted infections using estimated beta vector
fitted_deaths <- as.vector(X %*% B_hat)
fitted_infections <- as.vector(X_tilde %*% B_hat)

# Compare with actual data of daily deaths
true_deaths <- data$nhs       
day_of_2020 <- data$julian    

# Comparison plot: observed vs fitted deaths
initial_plot <- ggplot() +
  # True deaths
  geom_point(aes(x = day_of_2020, y = true_deaths, color = "Observed Deaths")) +
  # Fitted deaths and infections
  geom_line(aes(x = day_of_2020, y = fitted_deaths, color = "Fitted Deaths"), linewidth = 0.8) +
  geom_line(aes(x = f_seq, y = fitted_infections, color = "Daily Infections"), linewidth = 0.8) +
  labs( # label axes and title
    x = "Day of 2020",
    y = "Counts",
    color = "Legend",
    title = "Observed & Fitted Deaths with Daily Infections (Sanity Check)"
  ) +
  scale_color_manual( # manually assign colors
    values = c(
      "Observed Deaths" = "black",
      "Fitted Deaths" = "red",
      "Daily Infections" = "blue"
    )
  ) +
  theme_minimal()

initial_plot # display initial plot


######### SMOOTHING PARAMETER SELECTION ##########

# number of observations in dataset
n <- length(y)                  

## --- Poisson Log-Likelihood Function --- ##

# ll() function computes the log-likelihood for fitting our Poisson 
# deconvolution model, given our observed data and estimated beta vector. 
ll <- function(y, beta_hat, X) {
  
  # Inputs:
  #   y - vector of observed death counts
  #   beta_hat - vector of estimated B-spline coefficients
  #   X - convolution design matrix
  
  # Expected deaths under current parameters
  mu <- X %*% beta_hat
  
  # Poisson log likelihood 
  ll_pois <- sum(y*log(mu) - mu - lgamma(y+1))
  return(ll_pois)
}

## --- Grid search for optimal smoothing parameter based on minimised BIC --- ##

# Define grid of smoothing parameter values (transformed to original scale)
lambda_seq <- exp(seq(-13,-7,length=50))

# Initialise vector to store BIC values
BIC <- numeric(length(lambda_seq))
# Initialise empty matrix to store fitted gamma parameters
gammas <- matrix(NA, nrow = length(lambda_seq), ncol=80)

# Grid search over lambda values to find optimal smoothing parameter
for (i in seq_along(lambda_seq)){                                             

  lambda <- lambda_seq[i] # current smoothing parameter
  
  # Fit penalised Poisson model with current lambda using BFGS optimisation
  # using previously fitted gamma_hat parameters
  fit <- optim(par = gamma_hat, fn = pnll, gr = pnll_grad,
               y = y, X = X, S = S, lambda = lambda, method = "BFGS", 
               control=list(maxit=1000, reltol=1e-8)) # ensure convergence 
  
  # Extract estimated parameters and convert to original scale
  gamma_i <- fit$par
  gammas[i,] <- gamma_i
  beta_par <- exp(gamma_i)
  # Find fitted mu parameter
  mu_hat <- as.vector(X %*% beta_par)
  
  # Compute weight vector for Hessian calculation
  w <- as.vector(y/(mu_hat^2)) 
  
  # Compute Hessian matrices:
  H0 <- crossprod(X * sqrt(w)) # Log-Likelihood Hessian, without penalty
  H_lambda <- H0 + lambda * S # Penalised Hessian, with smoothing parameter
  
  # Compute Effective Degrees of Freedom (EDF)
  R <- chol(H_lambda)
  H_invH0 <- backsolve(R, backsolve(R, H0, transpose = TRUE))
  EDF <- sum(diag(H_invH0))
  
  # Bayesian Information Criterion (BIC) calculation
  BIC[i] <- -2*ll(y, beta_par, X) + log(n)*EDF
}

# Identify optimal lambda and corresponding fitted gamma that minimise BIC
BIC_min <- min(BIC)
lambda_bic <- lambda_seq[which(BIC==BIC_min)]
gamma_bic <- gammas[which(BIC==BIC_min),]


######### BOOTSTRAP UNCERTAINTY QUANTIFICATION #########
# For bootstrap resamples, each observation contributes proportionally to the
# weights of the pnll and its gradient, so we need to modify both

## --- Weighted Penalised Negative Log-Likelihood Function --- ##

# pnll_w() function finds the weighted version of previously defined pnll() 
# function for non-parametric bootstrapping. 
pnll_w <- function(y, gamma, X, S, lambda, w) {
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis functions
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter
  #   w - vector of bootstrap weights 
  
  B <- exp(gamma)               # transform to original scale
  mu <- X %*% B                 # expected deaths
  
  # Smoothing penalty term
  penalty <- 1/2*lambda* sum(B*(S %*% B))   
  
  # Weighted Poisson log-likelihood
  # (note log(y!) value dropped as it will not affect beta optimisation)
  ll_pois <- sum(w * (y*log(mu) - mu))
  
  # Penalised negative log-likelihood with weights
  pnll <- -ll_pois + penalty
  return(pnll)
}

## --- Gradient of Weighted Penalised Negative Log-Likelihood Function --- ##

# pnll_grad_w() function finds the weighted version of previously defined 
# pnll_grad() function for non-parametric bootstrapping.
pnll_grad_w <- function(y, gamma, X, S, lambda, w){
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis functions
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter
  #   w - vector of bootstrap weights 
  
  B <- exp(gamma)               # transform to original scale
  mu <- X %*% B                 # expected deaths
  
  # Smoothing penalty gradient 
  penalty_grad <- lambda * (S %*% B) * B
  
  # Gradient of the weighted log-likelihood
  F_value <- diag(as.vector(w * (y/mu - 1))) %*% X %*% diag(B) 
  ll_pois_grad <- colSums(F_value)

  # Gradient of penalised negative log-likelihood
  pnll_grad <- -ll_pois_grad + penalty_grad
  return(pnll_grad)
}

## --- Bootstrap resampling and refitting --- ##

# Number of bootstrap replicates for uncertainty estimation
n_bootstrap <- 200

# Initialise empty matrix to store bootstrap replicates
#   Each column represents a different bootstrap sample
#   Each row represents f(t) value for day t
f_hat_boot <- matrix(0, nrow = nrow(X_tilde), ncol = n_bootstrap)

# Perform bootstrapping 
for (b in 1:n_bootstrap) {
  
  # Generate non-parametric bootstrap weights 
  wb <- tabulate(sample(n, replace = TRUE), n)
  
  # Fit penalised model using bootstrap weights via BFGS method, using optimised 
  # lambda and gamma parameters (that previously produced minimised BIC).
  fit_b <- optim(par = gamma_bic, fn = pnll_w, gr = pnll_grad_w,       
                 y = y, X = X, S = S, lambda = lambda_bic, method = "BFGS", 
                 w = wb,                                # bootstrap weights
                 control=list(maxit=1000, reltol=1e-8)) # ensure convergence              
  
  # Extract fitted parameters and transform to original scale
  B_b <- exp(fit_b$par) 
  
  # Compute estimated daily infection for this bootstrap sample
  f_hat_boot[, b] <- X_tilde %*% B_b
}


######### FINAL VISUALISATION AND RESULTS #########

## --- Summary statistics from the bootstrap distribution --- ##

# 95% Confidence intervals from bootstrap samples (2.5% and 97.5% quantiles)
f_hat_int <- apply(f_hat_boot, 1, quantile, probs = c(0.025, 0.975)) 

# Transform gamma vector with minimised BIC to original scale
B_bic_hat <- exp(gamma_bic)   

# Re-calculate fitted deaths and fitted infections using optimised beta vector
fitted_deaths <- as.vector(X %*% B_bic_hat)                                         
fitted_infections_optimised <- as.vector(X_tilde %*% B_bic_hat)

## --- Final plot --- ##

# Plots observed daily death data alongside fitted deaths and fitted 
# infection curve from model with optimised beta vector
final_plot <- ggplot() +
  # True deaths
  geom_point(aes(x = day_of_2020, y = true_deaths, color = "Observed Deaths"), size = 0.8) +
  # Fitted deaths
  geom_line(aes(x = day_of_2020, y = fitted_deaths, color = "Fitted Deaths"), linewidth = 0.8) +
  # Daily infection rate with 95% CI
  geom_ribbon(aes(x = f_seq, ymin = f_hat_int[1,], ymax = f_hat_int[2,]), fill = "blue", alpha = 0.2) +
  geom_line(aes(x = f_seq, y = fitted_infections_optimised, color = "Daily Infections"), linewidth = 0.8) +
  labs( # label axes and title
    x = "Day of 2020",
    y = "Counts",
    color = "Legend",
    title = "Observed & Fitted Deaths with Daily Infections (w/ 95% CI)"
  ) +
  scale_color_manual( # manually assign colors
    values = c(
      "Observed Deaths" = "black",
      "Fitted Deaths" = "red",
      "Daily Infections" = "blue"
    )
  ) +
  theme_minimal()

final_plot  # display final plot


