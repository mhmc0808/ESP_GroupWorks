# Jackson Cramer (s2274544), Maximillian McCourt (s2145762), Natalia Montalvo Cabornero (s1969053)

# Exercise 1
n <- 1000
h_max <- 5

# works as expected, truncating may be inefficient
h <- rep(1:n, sample(1:h_max, n, replace=TRUE))[1:n] |> sample()
h


# Exercise 2

beta <- runif(n, 0, 1) # assign random sociability parameters
nc <- 15 # average no. contacts per person


# WE NEED TO ADJUST THIS FUNCTION TO MAKE IT SINGLE-LOOP

# Looped version of get.net
get.net_loop <- function(beta,h,nc=15){
    # initialise empty contact list with n entries
    contacts <- vector("list", n)
    # Loop through each pair of individuals
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            # if they are not in the same household
            if (h[i] != h[j]){
                # find probability of contact between person i and j
                prob_beta_ij <- (nc*beta[i]*beta[j])/(mean(beta)**2 * (n-1))
                # if a contact occurs, add each person to the others contact list
                if (sample(c(0,1), size=1, prob=c(1-prob_beta_ij, prob_beta_ij)) == 1){
                    contacts[[i]] <- c(contacts[[i]], j)
                    contacts[[j]] <- c(contacts[[j]], i)
                }
            }
        }
    }
    return(contacts)
}

alink <- get.net_loop(beta,h,nc)
print(alink)


# One looped version of get.net
get.net_loop <- function(beta,h,nc=15){
  n_beta <- length(beta) # !! notice how I am not referring to global variable n (like was done before)
  # initialise empty contact list with n entries
  contacts <- vector("list", n_beta)
  pairs <- combn(n_beta, 2)  # generates all possible combinations of (i, j) pairs
  # Loop through each pair of combinations
  for (k in 1:ncol(pairs)){
    i <- pairs[1, k]
    j <- pairs[2, k]
    # if they are not in the same household
    if (h[i] != h[j]){
      prob_beta_ij <- (nc * beta[i] * beta[j]) / (mean(beta)^2 * (n_beta - 1))
      if (sample(c(0,1), size=1, prob=c(1 - prob_beta_ij, prob_beta_ij)) == 1){
        contacts[[i]] <- c(contacts[[i]], j)
        contacts[[j]] <- c(contacts[[j]], i)
      }
    }
  }
  
  return(contacts)
}




# can also use rbinom(1,1, prob=prob_beta_ij) instead of sampling

# Exercise 3 - NATALIA'S

# SEIR model 
nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  ## SEIR stochastic simulation model.
  ## beta = transmission parameter
  ## h = vector indicating which household each person belongs to
  ## alink = list defining regular contacts of each person
  ## gamma = daily prob E -> I
  ## delta = daily prob I -> R
  ## nc = average number of contacts per person
  ## nt = number of days
  ## pinf = proportion of initial population to randomly start in state I
  
  x <- rep(0,n) ## initialize to susceptible state (= 0)
  
  ni <- max(1, round(pinf*n))   # at least one infected
  initial_infected <- sample(1:n, ni)  # randomly choosing by prob the first I
  x[initial_infected] <- 2 
  
  S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state
  ## initialize
  S[1] <- n-ni
  I[1] <- ni 
  
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    
    infected_id <- which(x==2)
    #####  S -> E
    ## x[x==0&u<beta*I[i-1]/n] <- 1 ## S -> E with prob beta*I[i-1]/n
    }
    
    S[i] <- sum(x==0)
    E[i] <- sum(x==1)
    I[i] <- sum(x==2)
    R[i] <- sum(x==3)

  return(list(S=S,E=E,I=I,R=R,beta=beta))
} 




# Ex 3 - MAX'S

# Function to simulate infection spread

nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
    # Our SEIR, t data
    # Using S=0, E=1, I=2, R=3 structure
    n <- length(beta)
    # defining alphas
    alpha_h <- alpha[1]
    alpha_c <- alpha[2]
    alpha_r <- alpha[3]
    # total population
    x <- rep(0, n)
    # initialise infected individuals at random
    x[sample(1:n, round(n*pinf))] <- 2
    # initialise SEIR counts
    S <- E <- I <- R <- rep(0,nt)
    S[1] <- sum(x == 0)
    I[1] <- sum(x == 2)
    # loop over days
    for (t in 2:nt){
        u <- runif(n)
        x[x==2 & u<delta] <- 3 # infected to recovered
        x[x==1 & u<gamma] <- 2 # exposed to infected
        # susceptible to exposed
        current_I <- which(x==2)
        for (i in current_I){
            # random mixing daily probability
            r_prob <- (alpha_r*nc*beta[i]*beta)/(mean(beta)**2 * (n - 1))
            x[x==0 & u<r_prob] <- 1
            # household mixing daily probability 
            # if x is susceptible, and in same household as i, and random number less than alpha_h probability
            household <- which(h==h[i] & x==0)
            x[household[u[household] < alpha_h]] <- 1
            # contacts mixing daily probability
            # if x is susceptible, and in contact list of i, and random number less than alpha_c probability
            contacts_i <- alink[[i]]
            x[contacts_i[x[contacts_i] == 0 & u[contacts_i] < alpha_c]] <- 1
            # count of SEIR population
        }
        S[t] <- sum(x==0); E[t] <- sum(x==1); I[t] <- sum(x==2); R[t] <- sum(x==3)
    }
    return(list(S=S,E=E,I=I,R=R,t=1:nt))
}


system.time(epi <- nseir(beta,h,alink))
# function operates correctly, runs in about 8.5 seconds with n=10,000.



# 4.

plot_dynamics = function(pop_states){
  # plot number of people in each group over time
  plot(pop_states$S,ylim=c(0,max(epi$S)),xlab="day", ylab="N", main="Simulated Population States over Time", cex.main=1) # S black
  points(pop_states$E, col=4)
  points(pop_states$I, col=2)
  points(pop_states$R, col=3)
}

plot_dynamics(epi)

