# Jackson Cramer (s2274544), Maximillian McCourt (s2145762), Natalia Montalvo Cabornero (s1969053)

# Our group worked collaboratively on the development and optimisation of all exercises.
# Although tasks were divided among us, we maintained a balanced workload through 
# regular code reviews, debugging sessions, and problem-solving meetings.
# Jackson focused on the plotting functions and household creation. 
# Max focused on implementing the nseir function.
# Natalia focused on the get.net function and commenting.



######### GENERAL DESCRIPTION #########

# This project extends the SEIR (Susceptible-Exposed-Infectious-Recovered) disease 
# model by incorporating household and social network structures. Individuals are 
# grouped into households and connected through probabilistic contact networks 
# based on sociability parameters. The model simulates disease transmission through 
# household, network, and random interactions to study how social structure and 
# individual variability affect epidemic dynamics.


######### SETUP ##########

# In this section, we define the population size and create the household and social 
# network structure. By generating these networks, this will enable us to model how 
# social structure influences the spread of infectious diseases.


## --- Distribution into households --- ##

n <- 10000 # Total population size
h_max <- 5 # Maximum number of people per household

# Randomly assign each person to a household of size 1 to h_max.
# Household IDs are repeated to give varying sizes, then shuffled
# to create an n-length vector representing a heterogeneous 
# household structure.
h <- rep(1:n, sample(1:h_max, n, replace=TRUE))[1:n] |> sample()


## --- Creation of contact network model --- ## 

# get.net() function constructs a contact network model among individuals.
# Links between individuals are formed probabilistically according to their
# sociability parameters (beta), excluding any individuals in the same household 
# (h). The output is a list where the i-th element contains the indices of person 
# i’s regular contacts.
get.net <- function(beta,h,nc=15){
  
  # Inputs:
  #   beta - vector of sociability parameters for each individual
  #   h - vector of household IDs for each individual
  #   nc - average number of contacts per person
  
  # Locally define population size
  n <- length(beta)
  
  # Initialise an empty contact list for all individuals
  contacts <- vector("list", n)
  
  # Compute the coefficient of our contact network formula
  coeff <- nc/(mean(beta)**2 * (n - 1))
  
  # Loop through all individuals to create network links
  for (i in 1:(n-1)){
    
    # For individual i, define all possible contacts of individuals i+1 to n
    # who are not in the same household (as previous links have already been generated)
    non_h <- h[(i+1):n] != h[i] 
    poss_links <- c((i+1):n)[non_h] 
    
    # Extract beta values (sociability parameters) for potential contacts
    beta_poss_links <- beta[poss_links]
    
    # Compute the probabilities of forming a contact between individual i and others using contact network formula
    prob_betas <- coeff*beta[i]*beta_poss_links
    
    # Generate links based on the assigned contact probabilities
    link_idx <- rbinom(length(prob_betas), 1, prob_betas) == 1
    links <- poss_links[link_idx]
    
    # Store the link in the i-th person's contact list
    contacts[[i]] <- c(contacts[[i]], links)
    
    # Append link with i-th person to the j-th person's contact list
    for (j in links) {
      contacts[[j]] <- c(contacts[[j]], i)
    }
  }
  
  # Return the full contact network as list
  return(contacts)
}


#########  SEIR MODEL  ##########

## --- SEIR simulation function --- ##

# nseir() function simulates a disease model over time in a population with 
# household and network-based social structure. Each individual is in one of four 
# states: Susceptible (S=0), Exposed (E=1), Infected (I=2), or Recovered (R=3). 
# The function uses constant daily probabilites for transitions from state E to I 
# and state I to R. For moving individuals from state S to E, we determine exposure 
# of susceptibles to the infectious population in state I by using the household 
# and regular contact networks previously defined. As a final method of exposure, 
# irrespective of household or regular network relations, we generate a daily probability 
# of exposure between each infected (I) and each susceptible (S) individual.
nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  
  # Inputs:
  #   beta - vector of sociability parameters for each individual
  #   h -  vector of household IDs for each individual
  #   alink - list defining the regular contacts of each person
  #   alpha - vector of household, regular network and random mixing parameters
  #   delta - daily probability of moving from state I to state R
  #   gamma - daily probability of moving from state E to state I
  #   nc - average number of contacts per person
  #   nt -  total number of days of simulation
  #   pinf - proportion of initial population to randomly start in the I state
  
  # Locally define population size
  n <- length(beta)
  
  # Defining alpha parameters
  alpha_h <- alpha[1] # Household parameter
  alpha_c <- alpha[2] # Contact network parameter
  alpha_r <- alpha[3] # Random daily mixing parameter
  
  # Recalling susceptible S=0, initialise total population as vector of zeros
  x <- rep(0, n)
  
  # Initialise infected individuals at random (using pinf input)
  x[sample(n, max(1,round(n*pinf)))] <- 2 # Use max() to ensure at least 1 person is initialised as infected
  
  # Initialise SEIR population counts for each day
  S <- E <- I <- R <- rep(0,nt)
  S[1] <- sum(x == 0) # Initialise population of susceptible individuals on day 1
  I[1] <- sum(x == 2) # Initialise population of infectious individuals on day 1
  
  # Compute the coefficient of our random daily mixing formula
  coeff <- (alpha_r*nc)/(mean(beta)^2 * (n-1))
  
  # Loop over days
  for (t in 2:nt){
    
    # Generate a random deviates for each individual for mutually exclusive transitions from states I to R and E to I
    u <- runif(n)
    x[x==2 & u<delta] <- 3 # Infectious to recovered state
    x[x==1 & u<gamma] <- 2 # Exposed to infectious state
    
    # Determine IDs of currently infectious population
    current_I <- which(x==2)
    
    # Ensure there is a current infectious population
    if (length(current_I) > 0) { 
      
      # HOUSEHOLD TRANSMISSION
      # Count infectious individuals per household using tabulate(), then map counts to each person
      num_I_h <- tabulate(h[x==2], nbins=max(h))[h]
      
      # Find indices of susceptible individuals living in households with infectious individuals
      sus_I_h <- which(x==0 & num_I_h>0)
      
      # Generate probabilities of each susceptible person becoming exposed using alpha_h parameter, adjusting for number of infectious people in household
      probs_h <- 1 - (1 - alpha_h)^num_I_h[sus_I_h]
      
      # Update population x for susceptible individuals who get exposed based on these generated probabilities
      x[sus_I_h[runif(length(sus_I_h)) < probs_h]] <- 1
      
      # CONTACT NETWORK TRANSMISSION
      # Combine all infectious individuals’ contact lists into one vector
      contacts_I <- unlist(alink[which(x==2)])
      
      # Count how many infectious contacts each person has with tabulate()
      num_I_c <- tabulate(contacts_I, nbins=length(x))
      
      # Identify indices of susceptible individuals in contact with at least one infectious person
      sus_w_I_c <- which(x==0 & num_I_c>0)
      
      # Generate probabilities of each susceptible person becoming exposed using alpha_c parameter, adjusting for number of infectious people in the susceptible individuals' contact network
      probs_c <- 1 - (1 - alpha_c)^num_I_c[sus_w_I_c]
      
      # Update population x for susceptible individuals who get exposed based on these generated probabilities
      x[sus_w_I_c[runif(length(sus_w_I_c)) < probs_c]] <- 1
      
      # RANDOM DAILY MIXING TRANSMISSION
      # Indices of susceptible individuals
      current_S <- which(x==0)
      
      # Ensure there is a current susceptible population
      if (length(current_S) > 0) {
        # Use outer() function to create all pairwise products of beta parameters for susceptible–infectious pairs
        beta_matrix <- outer(beta[current_S], beta[current_I])
        
        # Calculate probabilities of each susceptible avoiding infection from each infectious individual using random daily mixing formula
        prob_matrix_r <- 1 - coeff * beta_matrix
        
        # Combine row-wise to find total probability of each susceptible individual avoiding infection from all infectious individuals
        prob_not <- apply(prob_matrix_r, 1, prod)
        
        # Convert to infection probabilities for each susceptible
        probs_r <- 1 - prob_not
        
        # Update population x for susceptible individuals who get exposed based on these generated probabilities
        x[current_S[runif(length(current_S)) < probs_r]] <- 1
      }
    }
    
    # Record SEIR population for day t
    S[t] <- sum(x==0)
    E[t] <- sum(x==1)
    I[t] <- sum(x==2)
    R[t] <- sum(x==3)
  }
  
  # Return list of each SEIR population count from day 1 to day nt
  return(list(S=S,E=E,I=I,R=R,t=1:nt))
}


#########  VISUALISATION OF DISEASE MODEL ##########

## --- Plotting function --- ##

# plot_dynamics() function visualises the simulated evolution of the SEIR states 
# over time as a line plot, with each line in the plot coloured by their 
# corresponding SEIR state.
plot_dynamics <- function(pop_states, title=""){
  
  # Inputs:
  #   pop_states - population state counts as produced by nseir() function
  #   title - desired plot title
  
  # Calculate highest population among groups for consistent axis scaling
  ymax <- max(c(pop_states$S, pop_states$E, pop_states$I, pop_states$R))
  
  # Initialise plot and susceptible population over time as a black line
  plot(pop_states$S,ylim=c(0,ymax), type="l", lwd=3, xlab="Time (days)", ylab="Population", main=title)
  
  # Add horizontal grid lines at axis tick marks for easier visualisation
  abline(h=axTicks(2), col="gray80")
  
  # Plot remaining SEIR groups as unique colored lines
  lines(pop_states$E, col=4, lwd=3) # E - blue
  lines(pop_states$I, col=2, lwd=3) # I - red
  lines(pop_states$R, col=3, lwd=3) # R - green
  
  # Add legend to identify which line corresponds to which SEIR group
  legend(x=length(pop_states$R)/2, y=ymax/1.5,
         legend=c("Susceptible", "Exposed", "Infected", "Recovered"),
         col=c("black", "blue", "red", "green"),
         lty=NA, lwd=2,
         pch=16,
         cex=0.7, bty="n")
}


## --- Comparing Different Disease Model Scenarios --- ##

# We simulate and compare four scenarios to investigate the effects of household and contact network structure, 
# as well as the effect of individual variability in sociability (beta).

# First initialise a 'sociability' parameter (a random number between 0 and 1) for each individual
beta <- runif(n, 0, 1)

# Using get.net() function, create regular network of contacts between each individual
alink <- get.net(beta,h)

# 1. Default Parameters
# Simulate disease model with default parameters for household structure, network contacts, and random mixing
def_params <- nseir(beta,h,alink)

# 2. Random Mixing
# Simulate disease model by removing household and regular network structure,
# while keeping average initial infections per day the same for each person
# To do so, set alpha_h and alpha_c to 0, while setting alpha_r to .04
random_mixing <- nseir(beta, h, alink, alpha=c(0,0,.04))

# 3. Constant Beta
# Simulate disease model using default parameters, but assigning each person the same 
# sociability parameter, which is the average of our previously defined beta parameters
avg_beta <- rep(mean(beta), length(beta))
constant_beta <- nseir(avg_beta,h,alink)

# 4. Random Mixing & Constant Beta
# Simulate disease model by combining the conditions of the previous two simulations (random mixing and constant beta)
random_mixing_constant_beta <- nseir(avg_beta,h,alink, alpha=c(0,0,.04))


## --- Plotting the Four Simulations --- ##

# Set plotting layout/margins for disease model plots
par(mfcol=c(2,2), mar=c(4,4,2,1))

# Define list of plots and corresponding vector of titles
plots <- list(def_params, random_mixing, constant_beta, random_mixing_constant_beta)
titles <- c("Default Parameters", "Random Mixing", "Constant Beta", "Random Mixing & Constant Beta")

# Loop through each simulation and plot their disease model dynamics
for (i in seq_along(plots)) {
  plot_dynamics(plots[[i]], titles[i])
}

## --- Comments --- ##

# Including household and regular-contact structure and maintaining heterogeneity 
# in sociability noticeably reduces the speed, peak size and final attack rate 
# of the epidemic compared with a homogeneous random-mixing model. Removing 
# structure (random mixing) and removing β-heterogeneity both make the outbreak 
# more explosive; together they produce the fastest and largest epidemics. These 
# results are consistent with epidemic theory: clustering and individual heterogeneity 
# limit the efficiency of population-wide spread. To support quantitative claims, 
# the conclusions should be backed by replicate simulations and summary statistics 
# (peak size, time-to-peak and final attack rate).




Jackson's 

The default parameter plot represents a situation where individuals’ household and 
social network contacts affect their susceptibility to exposure. In contrast, the 
random mixing plot shows a case where these factors have no effect, and individuals
become exposed only based on their personal “sociability” parameters. When comparing 
the two plots, we see that the final population distributions are quite similar, 
meaning the household and contact network did not greatly change how many people 
ended up in each group. One key difference, though, is that the rate at which people
change states is higher in the random mixing case than in the default parameter case. 
This seems unusual, since we would generally expect that if there are more possible
exposure pathways, the rate of change from susceptible to recovered would be higher. 
However, we need to note that alpha_r was set to 0.01 in the default model and 0.04 
in the random mixing model. Therefore, the higher rate may be due to the increased 
random mixing probability rather than the presence of household and contact probabilities. 
This suggests that random mixing might have a stronger impact on the spread of the 
epidemic than the household or contact network.

We can also observe the effect of a constant beta on disease spread. With a constant
beta, the gap between the green and black lines at the end of the simulation is larger
than in the non-constant beta cases. This indicates that lower variance in sociability
increases the number of people who end up recovered and decreases the number still
susceptible. In other words, the disease spreads to more people when sociability is constant.


