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




# Vectorised version of get.net_vec
get.net_vec <- function(beta,h,nc=15){
    # lets collect all unique combinations of people i and j
    pairs <- t(combn(n, 2))
    # find all combinations where person i lives in the same household as person j
    h_matching <- h[pairs[,1]] == h[pairs[,2]]
    # we then remove all combinations where person i lives in the same household as person j
    pairs <- pairs[!h_matching,, drop=FALSE]
    # all probabilities of contact between person i and j
    probs_ij <- (nc*beta[pairs[,1]]*beta[pairs[,2]])/(mean(beta)**2 * (n-1))
    # now use probs_ij to generate all contacts between people 1 to n
    links <- rbinom(length(probs_ij), 1, probs_ij) == 1
    # we keep only the pairs that have a link from links
    pairs <- pairs[links,, drop=FALSE]
    # initialise empty contact list with n entries
    contacts <- vector("list", n)
    
    # Populate contact list with links for each person from 1 to n
    for (k in seq_len(nrow(pairs))) {
        contacts[[pairs[k, 1]]] <- c(contacts[[pairs[k, 1]]], pairs[k, 2])
        contacts[[pairs[k, 2]]] <- c(contacts[[pairs[k, 2]]], pairs[k, 1])
    }
    return(contacts)
}






# can also use rbinom(1,1, prob=prob_beta_ij) instead of sampling

# SEIR Model

nseir_old <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
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
  x[sample(n, max(1,round(n*pinf)))] <- 2 # use max() to ensure at least 1 person starts as infected
  # initialise SEIR counts
  S <- E <- I <- R <- rep(0,nt)
  S[1] <- sum(x == 0)
  I[1] <- sum(x == 2)
  # loop over days
  for (t in 2:nt){
    u <- runif(n)
    x[x==2 & u<delta] <- 3 # infected to recovered
    x[x==1 & u<gamma] <- 2 # exposed to infected
    susceptible <- which(x==0) # initialize those susceptible to infection
    # susceptible to exposed
    current_I <- which(x==2)
    for (i in current_I){
      u_hc <- runif(n)
      # if x is susceptible, and in same household as i, and random number less than alpha_h probability
      household <- which(h==h[i] & x==0)
      x[household[u_hc[household] < alpha_h]] <- 1
      # contacts mixing daily probability
      # if x is susceptible, and in contact list of i, and random number less than alpha_c probability
      contacts_i <- alink[[i]]
      x[contacts_i[x[contacts_i] == 0 & u_hc[contacts_i] < alpha_c]] <- 1
      # random mixing daily probability
      u_rm <- runif(n)
      r_prob <- alpha_r*nc*beta[i]*beta[susceptible]/(mean(beta)^2 * (n - 1))
      infected <- susceptible[u_rm[susceptible] < r_prob]
      x[infected] <- 1
      # count of SEIR population
    }
    S[t] <- sum(x==0)
    E[t] <- sum(x==1)
    I[t] <- sum(x==2)
    R[t] <- sum(x==3)
  }
  return(list(S=S,E=E,I=I,R=R,t=1:nt))
}







system.time(epi <- nseir(beta,h,alink))
# function operates correctly, runs in about 8.5 seconds with n=10,000.


# Exercise 4

plot_dynamics = function(pop_states, title=""){
  # plot number of people in each group over time
  ymax <- max(c(pop_states$S, pop_states$E, pop_states$I, pop_states$R))
  plot(pop_states$S,ylim=c(0,ymax), type="l", lwd=3, xlab="day", ylab="N", main=title) # S black
  lines(pop_states$E, col=4, lwd=3) # E blue
  lines(pop_states$I, col=2, lwd=3) # I red
  lines(pop_states$R, col=3, lwd=3) # R green
  legend(x=length(pop_states$R)/2, y=ymax/1.5, # length(pop_states$t)
         legend=c("Susceptible", "Exposed", "Infected", "Recovered"),
         col=c("black", "blue", "red", "green"),
         lty=NA, lwd=2,
         pch=16,
         cex=0.7, bty="n")
}

# Exercise 5

# default parameters
def_params = nseir(beta,h,alink)

# remove household and regular network structure
random_mixing = nseir(beta, h, alink, alpha=c(0,0,.04))

# beta vector set to contain average of previous veta for every element
avg_beta = 1:length(beta)
avg_beta[1:length(avg_beta)] = mean(beta)
constant_beta = nseir(avg_beta,h,alink)

# constant beta and random mixing
random_mix_constant_beta = nseir(avg_beta,h,alink, alpha=c(0,0,.04))

# plot all four scenarios side by side
par(mfcol=c(2,2), mar=c(4,4,2,1))
titles <- c("Default Parameters", "Random Mixing", "Constant Beta", "Random Mix & Constant Beta")
plots <- list(def_params, random_mixing, constant_beta, random_mix_constant_beta)
for (i in seq_along(plots)) {
  plot_dynamics(plots[[i]], titles[i])
}



