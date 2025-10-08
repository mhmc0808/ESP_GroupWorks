

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



contacts <- get.net_loop(beta,h,nc)
print(contacts)


# can also use rbinom(1,1, prob=prob_beta_ij) instead of sampling




# SEIR model 
nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  ## SEIR stochastic simulation model.
  ## n = population size; ni = initially infective; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## bmu = mean beta; bsc = var(beta) = bmu * bsc
  
  x <- rep(0,n) ## initialize to susceptible state
  beta <- rgamma(n,shape=bmu/bsc,scale=bsc) ## individual infection rates
  x[1:ni] <- 2 ## create some infectives
  S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state
  S[1] <- n-ni;I[1] <- ni ## initialize
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    x[x==0&u<beta*I[i-1]] <- 1 ## S -> E with prob beta*I[i-1]
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
} 



