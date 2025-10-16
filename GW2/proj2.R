# Jackson Cramer (s2274544), Maximillian McCourt (s2145762), Natalia Montalvo Cabornero (s1969053)

# Our group worked collaboratively on the development and optimisation of all exercises.
# Although tasks were divided among us, we maintained a balanced workload through 
# regular code reviews, debugging sessions, and problem-solving meetings.
# Jackson focused for the plotting functions and household creation. 
# Max focused on implementing the nseir function.
# Natalia focusedto the get.net function and commenting.


######### GENERAL DESCRIPTION #########

# This project extends the SEIR epidemic model by incorporating household and social 
# network structures. Individuals are grouped into households and connected through 
# probabilistic contact networks based on sociability parameters. The model simulates 
# disease transmission through household, network, and random interactions to study 
# how social structure and individual variability affect epidemic dynamics.



#########  SETTING UP  ##########

# In this section, we define the population structure and create social networks.
# The goal is to represent both household connections and external contacts
# in order to model how social structure influences the spread of infectious diseases.


## --- Distribution into Households --- ##

n <- 5000  # Total population size
h_max <- 5 # Maximum number of people per household

# Each individual is randomly assigned to a household of size between 1 and h_max.
# This creates a heterogeneous household structure.
h <- rep(1:n, sample(1:h_max, n, replace=TRUE))[1:n] |> sample()


## --- Creation of networks --- ## 

# get.net() constructs a social contact network among individuals.
# Each person i can form contacts with others who are not in the same household.
# beta = individual sociability
# h = household structure
# nc = expected number of contacts per person
get.net <- function(beta,h,nc=15){
  
  # n <- length(beta) ## do we need to keep this?? Always going to be = og n
  
  # Initialize an empty contact list for all individuals
  contacts <- vector("list", n)
  
  # Compute a coefficient to normalize the expected number of contacts
  coeff <- nc/(mean(beta)**2 * (n - 1))
  
  # Loop over all individuals to create network links
  for (i in 1:(n-1)){
    
    # Identify all possible contacts who are not in the same household
    non_h <- h[(i+1):n]!=h[i] 
    poss_links <- c((i+1):n)[non_h]  # initialise all possible links
    
    # Extract beta values (sociability parameters) for potential contacts
    beta_poss_links <- beta[poss_links]
    
    # Compute the probability of forming a contact between individual i and others
    prob_betas <- coeff*beta[i]*beta_poss_links
    
    # Generate links based on the assigned connection probabilities
    link_idx <- rbinom(length(prob_betas), 1, prob_betas) == 1
    links <- poss_links[link_idx]
    
    # Store the links in the contacts list
    contacts[[i]] <- c(contacts[[i]], links)
    for (j in links) {
      contacts[[j]] <- c(contacts[[j]], i)
    }
  }
  
  # Return the full contact network as list
  return(contacts)
}

#system.time(alink <- get.net(beta,h,nc))


# SEIR Model

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
  # initialize denominator for probability calculation
  coeff <- (alpha_r*nc)/(mean(beta)^2 * (n-1))
  # loop over days
  for (t in 2:nt){
    # generate n random deviates for mutually exclusive daily transitions
    u <- runif(n)
    x[x==2 & u<delta] <- 3 # infected to recovered
    x[x==1 & u<gamma] <- 2 # exposed to infected
    # susceptible to exposed
    current_I <- which(x==2)
    for (i in current_I){
      # if x is susceptible, and in same household as i, and random number less than alpha_h probability
      household <- which(h==h[i] & x==0)
      # uses runif of length of household to select which household members become exposed
      x[household[runif(length(household)) < alpha_h]] <- 1
      # contacts mixing daily probability
      # if x is susceptible, and in contact list of i, and random number less than alpha_c probability
      # finds all susceptible contacts of infected person i
      contacts_i <- alink[[i]]
      sus_contacts <- contacts_i[x[contacts_i]==0]
      # uses runif of length of susceptible contacts to select which susceptible contacts become exposed
      x[sus_contacts[runif(sus_contacts) < alpha_c]] <- 1
      # random mixing daily probability
      susceptible <- which(x==0) # initialize those susceptible to infection
      r_prob <- coeff*beta[i]*beta[susceptible]
      x[susceptible[runif(length(susceptible)) < r_prob]] <- 1
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

# build function to plot population in each S.E.I.R group over time
plot_dynamics = function(pop_states, title=""){
  # pop_states - population state counts as produced by nseir
  # title - plot title
  # calculate highest population among groups for consistent axis scaling
  ymax <- max(c(pop_states$S, pop_states$E, pop_states$I, pop_states$R))
  # plot each S.E.I.R group as lines, differentiated by color
  plot(pop_states$S,ylim=c(0,ymax), type="l", lwd=3, xlab="Time (days)", ylab="Population", main=title) # S black
  lines(pop_states$E, col=4, lwd=3) # E blue
  lines(pop_states$I, col=2, lwd=3) # I red
  lines(pop_states$R, col=3, lwd=3) # R green
  # include legend to recognize meaning of each line color
  legend(x=length(pop_states$R)/2, y=ymax/1.5, # length(pop_states$t)
         legend=c("Susceptible", "Exposed", "Infected", "Recovered"),
         col=c("black", "blue", "red", "green"),
         lty=NA, lwd=2,
         pch=16,
         cex=0.7, bty="n")
}


# Exercise 5

# initialize a 'sociability' parameter (a number between 0 and 1) for each person
beta <- runif(n, 0, 1)
# create regular network of contacts between each person
alink <- get.net(beta,h)

# simulate epidemic using default parameters as specified in nseir function
def_params = nseir(beta,h,alink)

# simulate epidemic by removing household and regular network structure, keeping average initial infections per day the same for each person
random_mixing = nseir(beta, h, alink, alpha=c(0,0,.04)) # ??!! double check I did this right, and update the 4th if its wrong

# simulate epidemic using full model with default parameters, but assigning each person the same 'sociability' parameter, which is the average of our original beta parameters
avg_beta = 1:length(beta)
avg_beta[1:length(avg_beta)] = mean(beta)
constant_beta = nseir(avg_beta,h,alink)

# simulate epidemic by combining the conditions of the previous two simulations
random_mix_constant_beta = nseir(avg_beta,h,alink, alpha=c(0,0,.04))

# set plotting layout/margins for epidemic plots
par(mfcol=c(2,2), mar=c(4,4,2,1))
# specify titles to be added to plots
titles <- c("Default Parameters", "Random Mixing", "Constant Beta", "Random Mix & Constant Beta")
plots <- list(def_params, random_mixing, constant_beta, random_mix_constant_beta)
# use plot_dynamics to plot each of the four epidemic conditions
for (i in seq_along(plots)) {
  plot_dynamics(plots[[i]], titles[i])
}

# Comment on plots



