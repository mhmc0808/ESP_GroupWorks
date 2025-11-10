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

setwd("C:/Users/natmo/OneDrive/Escritorio/MSc/SEM1/Extended Statistical Programming/Assignments/GW3")
library(splines)
library(ggplot2)
data <- read.table("engcov.txt")


set.seed(3)  # remove

######### MODEL SET UP ##########

# Define the matrices and components of the deconvolution model, the function 
# evaluate_matrices() creates: B-spline basis matrix, Convolution matrix (X)  
# that maps infections to deaths via delay distribution and Penalty matrix (S) 
# for smoothing. 

evaluate_matrices <- function(data, k = 80, edur = 3.151, sdur = 0.469, lower_bound, 
                              upper_bound, f_seq) {
  
  # Inputs;
  #   k - number of B-spline basis functions
  #   edur - log(mean) time from infection to death
  #   sdur - log(standard deviation) time from infection to death
  
  # Create penalty matrix S
  S <- crossprod(diff(diag(k), diff=2))
  
  # Create the k+4 knots for B-splines
  # Sequence of k+4 evenly spaced knots
  middle_seq_length <- k-2
  middle_seq <- seq(lower_bound, upper_bound, length.out = middle_seq_length)
  
  # Step size for handling remaining knots
  step <- middle_seq[2] - middle_seq[1]
  
  # Knots sequence
  knots <- c(
    seq(to = min(middle_seq) - step, by = step, length.out = 3),  # 3 knots below
    middle_seq,                                                   # middle sequence
    seq(from = max(middle_seq) + step, by = step, length.out = 3) # 3 knots above
  )
  
  # Create X_tilde: B-spline basis matrix for infection function f(t)
  X_tilde <- splineDesign(knots, f_seq)
  
  # Delay distribution from infection to death (pd)
  d <- 1:k
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)      # normalise the sum to 1
  
  n <- length(data$date)  # number of death days
  m <- length(pd)         # maximum of days to consider    
  r <- nrow(X_tilde)      # number of infection days
  
  # Build weight matrix W for convolution
  # W[i, t] = pd[30 + i - t] if in valid range, 0 otherwise
  W <- matrix(0, n, r)
  
  for (i in 1:n) {
    idxs <- (30 + i - (1:m))              # indices for 30 day lag
    valid <- which(idxs > 0 & idxs <= r)  # keep indices in valid range
    # take i-th row of W at valid indexes to be pd of valid index
    W[i, idxs[valid]] <- pd[1:length(valid)]
  }
  
  # Convolution matrix X maps infections to expected deaths
  X <- W %*% X_tilde
  
  return(list(X_tilde=X_tilde,S=S,X=X))
}

## --- Generate model matrices --- ##

# Days of infection, extending the sequence 30 days before the first observed death
lower_bound <- min(data$julian) - 30
upper_bound <- max(data$julian)
f_seq <- seq(lower_bound, upper_bound, by=1)

model_matrices <- evaluate_matrices(data,lower_bound=lower_bound, 
                                    upper_bound=upper_bound, f_seq=f_seq)
X <- model_matrices$X
X_tilde <- model_matrices$X_tilde
S <- model_matrices$S



######### PENALISED NEGATIVE LOG-LIKELIHOOD FUNCTIONS ##########

# We are going to fit the model using a penalised Poisson likelihood with a 
# positive constraint

## --- Penalised negative log-likelihood --- ##
pnll <- function(y, gamma, X, S, lambda) {
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter

  B <- exp(gamma)               # transformation of coefficients
  mu <- pmax(X %*% B)           # expected deaths
  
  # Smoothing penalty term
  penalty <- 1/2*lambda * sum(B * (S %*% B))     # OPTIMISED VERSION
  # penalty <- 1/2*lambda*t(B) %*% S %*% B
  
  # Poisson log likelihood
  ll_pois <- sum(y*log(mu) - mu - lfactorial(y))
  
  # Penalised negative log likelihood
  pnll <- -ll_pois + penalty
  return(pnll)
}


## --- Gradient of penalized negative log-likelihood --- ##
pnll_grad <- function(y, gamma, X, S, lambda){
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter
  
  B <- exp(gamma)
  mu <- pmax(X %*% B)                 # expected deaths
  
  # Gradient of log-likelihood: dl/dgamma = colsums of diag(y_i/mu_i - 1)*X*diag(B)
  F_value <- diag(as.vector(y/mu - 1)) %*% X %*% diag(B)            
  ll_pois_grad <- colSums(F_value)
  
  # Smoothness penalty gradient
  penalty_grad <- lambda * (S %*% B) * B    # OPTIMISED
  
  # Derivative of penalised negative log likelihood
  pnll_grad <- - ll_pois_grad + penalty_grad
  
  return(pnll_grad)
}


## --- Test functions --- ##
y <- data$deaths
lambda <- 5*10^-5
k <- 80
gamma <- rep(0,k)
pnll(y, gamma, X, S, lambda)
pnll_grad(y, gamma, X, S, lambda)


## --- Gradient validation uing finite differences--- ##

th0 <- gamma                   # test paramater vector
fd <- numeric(length(th0))     # vector to store the approximation

# Penalised negative log-likelihood at th0
nll0 <- pnll(y = y, gamma = th0, X = X, S = S, lambda = lambda) 

eps <- 1e-7                    # perturbation size

# Loop over each parameter
for (i in 1:length(th0)) { 
  th1 <- th0
  th1[i] <- th1[i] + eps
  # Pnll at that perturbed parameter
  nll1 <- pnll(y = y, gamma = th1, X = X, S = S, lambda = lambda)
  # Approximate the partial derivative with the finite difference
  fd[i] <- (nll1 - nll0) / eps
}

# Find analytic gradient at th0
grad_test <- pnll_grad(y = y, gamma = th0, X = X, S = S, lambda = lambda)

# Finite-difference vs analytic gradient
cbind(fd, grad_test) 

# Maximum absolute difference 
max_diff <- max(abs(grad_test-fd)) 
print(paste("Maximum gradient difference:", round(max_diff,8)))
#  0.00027709



######### INITIAL MODEL FITTING ##########

# Fit model using BFGS with initial lambda for sanity check.
# This provides a preliminary fit to verify the model is working correctly
# before proceeding with smoothing parameter selection.


fit <- optim(par = gamma,         # initial parameter values
             fn = pnll,           # penalized negative log-likelihood
             gr = pnll_grad,      # gradient function for efficient optimization
             y = y,               # observed death counts
             X = X,               # convolution matrix
             S = S,               # penalty matrix
             lambda = lambda,  # initial smoothing parameter
             method = "BFGS")     # Quasi-Newton optimization method with gradient

# Extract estimated parameters after optimization
gamma_hat <- fit$par
B_hat <- exp(gamma_hat)

# Compute fitted deaths and fitted infections
fitted_deaths <- as.vector(X %*% B_hat)
fitted_infections <- as.vector(X_tilde %*% B_hat)

# Compare with actual data
true_deaths <- data$deaths        # actual deaths from data set
day_of_2020 <- data$julian

# Create data frame for plotting
plot_df <- data.frame(
  day = day_of_2020,
  true_deaths = true_deaths,
  fitted_deaths = fitted_deaths)

# Comparison plot: observed vs fitted deaths
initial_plot <- ggplot() +
  # True deaths
  geom_line(aes(x = day_of_2020, y = true_deaths, color = "Observed Deaths"), size = 0.8) +
  # Fitted deaths
  geom_line(aes(x = day_of_2020, y = fitted_deaths, color = "Fitted Deaths"), size = 0.8) +
  geom_line(aes(x = f_seq, y = fitted_infections, color = "Daily Infections"), size = 0.8) +
  # Secondary axis
  #scale_y_continuous(
  #  name = "Number of Deaths",
  #  sec.axis = sec_axis(~ ., name = "Daily New Infections")
  #) +
  labs(
    x = "Day of 2020",
    color = "Legend",
    title = "Observed & Fitted Deaths with Daily Infections (95% CI)"
  ) +
  theme_minimal()

initial_plot

######### SMOOTHING PARAMETER SELECTION ##########

n <- length(y)                       # number of observations
lfy <- lfactorial(y)

# Function to compute the Poisson log-likelihood for BIC calculation
ll <- function(y, beta_hat, X, lfly) {
  
  # Inputs:
  #   y - vector of observed death counts
  #   beta_hat - vector of estimated B-spline coefficients
  #   X - convolution design matrix
 
  # Expected deaths under current parameters
  mu <- X %*% beta_hat
  
  # Poisson log likelihood 
  ll_pois <- sum(y*log(mu) - mu - lfy)
  return(ll_pois)
}

# Define sequence of lambda values to test
log_lambda_seq <- seq(-13,-7,length=50)
lambda_seq <- exp(log_lambda_seq)    # convert to actual lambda values

BIC <- numeric(length(lambda_seq))   # initialize BIC storage vector

Xt <- t(X)  # transpose of X to save some time in the loop

# Grid search over lambda values to find optimal smoothing parameter
for (i in seq(lambda_seq)){
  
  # Current lambda value
  lambda <- lambda_seq[i]
  
  # Fit model with current lambda using BFGS optimization
  fit <- optim(par = gamma,# gamma_hat???           # initial parameter values
               fn = pnll,             # penalized negative log-likelihood
               gr = pnll_grad,        # gradient function
               y = y,                 # observed death counts
               X = X,                 # convolution matrix
               S = S,                 # penalty matrix
               lambda = lambda,       # initial smoothing parameter
               method = "BFGS")       # Optimization method
            
  # Extract estimated parameters
  beta_par <- exp(fit$par)
  mu_hat <- as.vector(X %*% beta_par)
  
  # Weight vector
  w <- as.vector(y/(mu_hat^2)) # n vector 
  
  # Compute Hessian matrices:
  # H0: Hessian of log-likelihood (without penalty)
  H0 <- Xt %*% (X * w)   # X^T diag(w) X
  
  #  H_lambda: Full Hessian (log-likelihood + penalty)
  H_lambda <- H0 + lambda * S
  
  ## Effective Degrees of Freedom (EDF) calculation
  R <- chol(H_lambda)
  H_invH0 <- backsolve(R, backsolve(R, H0, transpose = TRUE))
  EDF <- sum(diag(H_invH0))
  
  # Bayesian Information Criterion (BIC) calculation
  BIC[i] <- -2*ll(y, beta_par, X, lfy) + log(n)*EDF
}

# Find optimal lambda that minimizes BIC
BIC_min <- min(BIC)
lambda_bic <- lambda_seq[which.min(BIC)]
lambda_bic



######### BOOTSTRAP UNCERTAINTY QUANTIFICATION #########
# For bootstrap resamples, each observation contributes proportionally to the
# weights of the Pnll and its gradient, so we need to modify both


# Weighted version of penalized negative log-likelihood for bootstrap
pnll_w <- function(y, gamma, X, S, lambda, w) {
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter
  #   w - vector of bootstrap weights 
  
  B <- exp(gamma)               # transformation of coefficients
  mu <- pmax(X %*% B)           # expected deaths, avoiding zero
  
  # Smoothing penalty term
  penalty <- 1/2*lambda* sum(B*(S %*% B))   # OPTIMISED VERSION

  # Weighted Poisson log-likelihood
  ll_pois <- sum(w * (y*log(mu) - mu - lfactorial(y))) 
  
  # Penalized negative log-likelihood with weights
  pnll <- -ll_pois + penalty
  return(pnll)
}

# Weighted version of gradient function for bootstrap
pnll_grad_w <- function(y, gamma, X, S, lambda, w){
  
  # Inputs:
  #   y - vector of observed deaths
  #   gamma - vector of transformed parameters: beta = exp(gamma) to ensure 
  #           positive coefficients for the B-spline basis
  #   X - convolution matrix
  #   S - penalty matrix
  #   lambda - smoothing parameter
  #   w - vector of bootstrap weights 
  
  B <- exp(gamma)               # transformation of coefficients
  mu <- pmax(X %*% B)           # expected deaths

  # Gradient of the weighted log-likelihood
  F_value <- diag(as.vector(w * (y/mu - 1))) %*% X %*% diag(B) 
  ll_pois_grad <- colSums(F_value)
  
  # Gradient of penalty
  penalty_grad <- lambda * (S %*% B) * B   # OPTIMISED VERSION
  
  # Penalised negative log likelihood
  pnll_grad <- -ll_pois_grad + penalty_grad
  return(pnll_grad)
}

# -- Bootstrap resampling and refitting


# Number of bootstrap replicates for uncertainty estimation
n_bootstrap <- 200

# Initialise empty matrix to store bootstrap replicates
#   Each column represents a different bootstrap sample
#   Each row represents f(t) at a specific time point across all bootstrap samples
f_hat_boot <- matrix(0, nrow = nrow(X_tilde), ncol = n_bootstrap)


for (b in 1:n_bootstrap) {
  
  # Generate bootstrap weights 
  wb <- tabulate(sample(n, replace = TRUE), n)
  
  # Fit penalized model using bootstrap weights
  fit_b <- optim(par = gamma,          # initial parameter values
               fn = pnll_w,            # weighted penalized negative log-likelihood
               gr = pnll_grad_w,       # gradient function for efficient optimization
               y = y,                  # observed death counts
               X = X,                  # convolution matrix
               S = S,                  # penalty matrix
               lambda = lambda_bic,    # initial smoothing parameter
               method = "BFGS",        # Quasi-Newton optimization method with gradient
               w = wb)                 # weights for this iteration

  # Extract fitted parameters from bootstrap sample
  gamma_b <- fit_b$par
  B_b <- exp(gamma_b)
  
  # Compute infection rate estimate for this bootstrap sample
  f_hat_boot[, b] <- X_tilde %*% B_b
  
  gamma <- gamma_b 
}


######### FINAL VISUALIZATION AND RESULTS #########

## --- Summary statistics from the bootstrap distribution --- ##
# Mean infection curve across bootstrap samples
#f_hat_mean <- rowMeans(f_hat_boot) 
f_hat_lower <- apply(f_hat_boot, 1, quantile, probs = 0.025) # 2.5% quantile
f_hat_upper <- apply(f_hat_boot, 1, quantile, probs = 0.975) # 97.5% quantile

# Re-fit infections with lambda_bic
fit <- optim(par = gamma, 
             fn = pnll,
             gr = pnll_grad,
             y = y, X = X, S = S, 
             lambda = lambda_bic, 
             method = "BFGS")
gamma_hat <- fit$par
B_hat <- exp(gamma_hat)
fitted_infections_optimized <- as.vector(X_tilde %*% B_hat)

# Create data frame for infection curve with confidence intervals
plot_df_6 <- data.frame(
  day = f_seq,               # Day of year for infection estimates
  f_mean = f_hat_mean,       # Mean across bootstraps
  f_lower = f_hat_lower,     # Lower bound of 95% confidence interval
  f_upper = f_hat_upper     # Upper bound of 95% confidence interval
)

# Final plot showing:
#   Observed death data
#   Model with fitted deaths
#   Estimated infection curve

final_plot <- ggplot() +
  # True deaths
  geom_line(aes(x = day_of_2020, y = true_deaths, color = "Observed Deaths"), size = 0.8) +
  # Fitted deaths
  geom_line(aes(x = day_of_2020, y = fitted_deaths, color = "Fitted Deaths"), size = 0.8) +
  # Daily infection rate with 95% CI
  geom_ribbon(aes(x = f_seq, ymin = f_hat_lower, ymax = f_hat_upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(x = f_seq, y = fitted_infections_optimized, color = "Daily Infections"), size = 0.8) +
  # Secondary axis: rescale infections to match the main axis numerically
  #scale_y_continuous(
  #  name = "Number of Deaths",
  #  sec.axis = sec_axis(~ ., name = "Daily New Infections")
  #) +
  labs(
    x = "Day of 2020",
    color = "Legend",
    title = "Observed & Fitted Deaths with Daily Infections (95% CI)",
    subtitle = "Infection curve shows estimated daily new infections that 
    led to observed deaths"
  ) +
  theme_minimal()

final_plot

# Interpretation:
# A we can see in the graph the close match between the observed and fitted deaths 
# demonstrates that our model successfully finds the mortality parameter in our data, 
# providing confidence in the model to reconstruct infection dinamics.
#
# The blue infection curve shows the estimated daily new infections that resulted in
# the observed deaths. We can clearly see the infection waves, showing when transmission
# was increasing or decreasing, like a big increase around day 75 and a big 
# decrease around day 85.
#
# The blue shaded region represents the 95% confidence interval around the infection
# estimates, derived from 200 bootstrap replicates. Wider intervals indicate greater
# uncertainty, particularly at the boundaries where data was limited.
# 
# The approximately 30-day horizontal shift that we assumed when creating the model
# between infection peaks and corresponding death peaks visually confirms the 
# lognormal delay distribution used in the model.
#
# The overall visualization demonstrates successful deconvolution, we have effectively
# recovered the unobserved infection trajectory from the observed death data using
# the known delay distribution and smoothness assumptions.


# Natalia's notes:
# - added a description for the project
# - added comments and descriptions for all the functions
# - added interpretation of final plot
# - changed some of the matrices multiplications for optimization
# - changed EFD calculation to chol decomp for a faster approach maybe? -- not sure if actually makes a difference
# - added some sanity checks so some our values cannot be zero to avoid crashes
# - we need to think/ask tutors how are we presenting the report because
#   actually the first 3 exercises are just steps for alter building the model with
#   the weights and the right lambda
# - changed last plot --> error when running jackson's plots
# - need to check if we can optimize lambda loop and bootstrap loop 
# My code runs in ~ 50 sec
