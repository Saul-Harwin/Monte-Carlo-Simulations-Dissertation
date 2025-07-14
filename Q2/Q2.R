library(markovchain)
library(igraph)
library(ggplot2)
library(pracma)

# Set seed for reproducibility
set.seed(963)

##### Section 1 - Define Markov Chain #####
"First section defines the transition_matrix"

n_states <- 5

# Define the transition matrix
transition_matrix <- matrix(runif(n_states * n_states), nrow = n_states)
transition_matrix <- t(apply(transition_matrix, 1, function(x) x / sum(x))) # Normalize each row to sum to 1

rownames(transition_matrix) <- colnames(transition_matrix) <- as.character(1:n_states)
# Print the row sums of the transition matrix
cat("Row sums =", rowSums(transition_matrix), "\n")



##### Section 2 - Define Simulate Markov Chain Function ######

simulate_chain <- function(P, start, n) {
  # Defines a chain as a vector with length n.
  chain <- numeric(n)
  
  # Defines the first element of the vector as the starting state.
  chain[1] <- start
  
  # Loops through the chain from 2 to n+1. k = step integer.
  for (k in 2:(n+1)) {
    states = 1:nrow(P)
    next_state_probs = P[chain[k-1], ]
    
    # Generates 1 random sanple from states each with probabilities of 
    # next_state_probs. Defines the chain at step k as this state
    chain[k] <- sample(states, 1, prob = next_state_probs)
    
  }
  return(chain)
}

##### Section 4 - Markov Chain Tests #####
"In this section we test irreducibility and aperiodicity"
"Both irreducibility and aperiodicity are satisfied"

is_irreducible <- function(P, max_steps = 100) {
  "
  Testing Irreducibility
  
  Parameters:
    P : matrix
      The transition matrix
      
    max_steps : int
      The amount of steps to simulate.
      
  Returns:
    all_reachable : bool
      True if the probability of reaching any state from any other state is > 0.
  "
  # Define the number of states
  n <- nrow(P)
  
  # Define the reach matrix as the identity matrix representing reachability with 0 steps
  reach_matrix <- diag(n)
  
  # Starts with 1 step transition matrix
  current_P <- P
  
  # Loops from 1 step to max_steps
  for (i in 1:max_steps) {
    # Add current_P to reach_matrix. Current P represents the probability of 
    # reaching every other state from current state. 
    reach_matrix <- reach_matrix + current_P
    
    # Define current_P as the transition matrix of one more step ahead
    current_P <- current_P %*% P
  }
  
  # Define all_reachable. every element in the reach_matrix > 0 then the 
  # the probability of reaching any state from any other state is > 0.
  all_reachable <- all(reach_matrix > 0)
  return(all_reachable)
}
print(is_irreducible(transition_matrix))


is_aperiodic <- function(P, max_steps = 100) {
  "
  Testing Aperiodicity
  
  Parameters:
    P : matrix
      The transition matrix
      
    max_steps : int
      The amount of steps to simulate.
      
  Returns:
    all_periods_1 : bool
      True if the probability of reaching any state from any other state is > 0.
  "
  
  # Define the number of states
  n <- nrow(P)
  
  # Defines array to store the period for each state
  periods <- rep(0, n)
  
  # Loop through each state
  for (i in 1:n) {
    
    # Define an array which will record which steps returned to state i. Defined at start every step false.
    visited <- rep(FALSE, max_steps)
    
    # Defines prob as a vector of zeros 
    prob <- rep(0, n)
    
    # Sets the state i position to 1. Representing 100% probability.
    prob[i] <- 1

    # Simulate markov chain upto max steps.
    for (k in 1:max_steps) {
      
      # Each step updates probabilities.
      prob <- prob %*% P
      
      # If there is a non-zero probability of moving to state i set the step k element of the  vistited array to True.
      if (prob[i] > 0) visited[k] <- TRUE
    }
    
    # Calculate the period of state i with the GCD function
    periods[i] <- Reduce(gcd, which(visited))
  }
  
  # If for all state the period is 1 then the Markov chain is aperiodic.
  all_periods_1 = all(periods == 1)
  return(all_periods_1)
}

print(is_aperiodic(transition_matrix))



##### Section 5 - Compute the invariant distribution ######
"Compute the invariant distribution. See chapter 5.4.3"

get_invariant_distribution <- function(P) {
  # Define the transpose of P
  P_transposed <- t(P)
  
  # Computes the eigenvalues and eigenvectors
  eig <- eigen(P_transposed)
  
  # Finds the index of the eigenvalue closest to 1. Isn't finding index of the 
  # eigenvalue which is equal to 1 to avoid precision errors.
  index_of_1 <- which.max(Re(eig$values))
  
  # Get eigenvector.
  invariant_vector <- eig$vectors[, index_of_1]
  
  # Removes any imaginary part.
  invariant_vector <- Re(invariant_vector)
  
  # Normalise the eingenvector.
  invariant_distribution <- invariant_vector / sum(invariant_vector)

  return(invariant_distribution)
}


invariant_distribution = get_invariant_distribution(transition_matrix)

print(invariant_distribution)

##### Section 6 - Compare the Distributions #####
"Compare the Empirical Distribution with Theoretical Invariant Distribution"

# Define the simulated chains parameters
start = 3
n = 1000

states <- simulate_chain(transition_matrix, start, n)
empirical_distribution <- table(states) / length(states)

# Plot empirical distribution vs theoretical invariant distribution
barplot(rbind(invariant_distribution, empirical_distribution),
        beside = TRUE, col = c("red", "blue"),
        names.arg = paste0("State ", 1:5),
        legend.text = c("Theoretical", "Empirical"),
        main = "Comparison of Theoretical and Empirical Distribution")
