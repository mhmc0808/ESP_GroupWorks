

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
