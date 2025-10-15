# Jackson Cramer (s2274544), Maximillian McCourt (s2145762), Natalia Montalvo Cabornero (s1969053)

# Exercise 1
n <- 1000
h_max <- 5

h <- rep(1:n, sample(1:h_max, n, replace=TRUE))[1:n] |> sample()

# Exercise 2

# One looped version of get.net
get.net <- function(beta,h,nc=15){
  # establish population using length of beta
  n <- length(beta) 
  # initialise empty contact list with n entries
  contacts <- vector("list", n)
  # record probability coefficient for efficiency
  coeff <- nc/(mean(beta)**2 * (n - 1))
  for (i in 1:(n-1)){
    # establish all people in future possible links that are not household members to person i
    non_h <- h[(i+1):n]!=h[i] # ??!! how does this handle when there are non-household members earlier than i
    # initialise all possible links for person i
    poss_links <- c((i+1):n)[non_h]
    # extract corresponding betas
    beta_poss_links <- beta[poss_links]
    # finds probabilities of contact between person i and possible future links
    prob_betas <- coeff*beta[i]*beta_poss_links
    # generate links based on probabilities
    link_idx <- rbinom(length(prob_betas), 1, prob_betas) == 1
    links <- poss_links[link_idx]
    contacts[[i]] <- c(contacts[[i]], links)
    for (j in links) {
      contacts[[j]] <- c(contacts[[j]], i)
    }
  }
  return(contacts)
}

#system.time(alink <- get.net(beta,h,nc))


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



#system.time(epi <- nseir(beta,h,alink))
# function operates correctly, runs in about 8.5 seconds with n=10,000.


# Exercise 4

plot_dynamics = function(pop_states, title=""){
  # plot number of people in each group over time
  ymax <- max(c(pop_states$S, pop_states$E, pop_states$I, pop_states$R))
  plot(pop_states$S,ylim=c(0,ymax), type="l", lwd=3, xlab="Time (days)", ylab="Population", main=title) # S black
  lines(pop_states$E, col=4, lwd=3) # E blue
  lines(pop_states$I, col=2, lwd=3) # I red
  lines(pop_states$R, col=3, lwd=3) # R green
  legend(x=length(pop_states$R)/2, y=ymax/1.5,
         legend=c("Susceptible", "Exposed", "Infected", "Recovered"),
         col=c("black", "blue", "red", "green"),
         lty=NA, lwd=2,
         pch=16,
         cex=0.7, bty="n")
}


# Exercise 5

beta <- runif(n, 0, 1) # assign random sociability parameters
alink <- get.net(beta,h,nc)

# default parameters
def_params = nseir(beta,h,alink)

# remove household and regular network structure
random_mixing = nseir(beta, h, alink, alpha=c(0,0,.04))


# beta vector set to contain average of previous beta for every element
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


# Comment on plots


