library(markovchain)
library(igraph)
library(ggplot2)
library(reshape2)

# Set seed for reproducibility
set.seed(963)

##### Section 1 - Define 1D Transition Matrix Creation Function #####
"First section defines the transition_matrix"

create_transition_matrix_1d <- function(Nx) {
  "
  Creates a 1-dimmensional transition matrix
  Parameters: 
    Nx : int
      Length of the chain
  
  Returns:
    P : matrix
      Transition matrix.
  "
  
  N <- Nx
  P <- matrix(0, nrow = N, ncol = N)
  
  for (i in 1:Nx) {
    k <- i
    
    if (i == 1 || i == Nx) {
      # Boundary: absorbing state
      P[k, k] <- 1
    } else {
      
      # Interior point → connect to 4 neighbors (no self-loop)
      neighbors <- c(
        i - 1,  # left
        i + 1  # right
      )
      P[k, neighbors] <- 1 / length(neighbors)
    }
  }
  
  return(P)
}

##### Section 2 - Define 1D Simulate Chain Function #####
"This section simulates a single chain how it moves between states over N steps"

# To generate a sample chain we need to generate random walks on the graph. We will now define a function to generate sample chains:
simulate_chain_1d <- function(start_k, P, Nx, f, h) {
  "
  Parameters:
    start_k : int
      Starting state for the chain
    P : matrix
      Transition matrix
    Nx : int
      Size of lattice in x direction.
    f : function
      The Source Function. For Leplace f(k) = 0
    h : float
      Real space distance between lattice points
      
  Returns:
    [boundary_point, F_w, path] : list
      boundary_point : int
        The Boundary point that the path stopped at.
      F_w : float
        Total accumulated contribution of the source term along the walk.
      path : array
        Array of states that the chain walked.
  "
  
  
  k <- start_k
  F_w <- 0
  path <- c(k)
  
  repeat {
    i=k
    
    # Check if boundary point (absorbing state)
    if (i == 1 || i == Nx) {
      break
    }
    
    # Accumulate the f(i,j) * h^2 term
    F_w <- F_w + as.numeric(f(k)) * as.numeric(h)*as.numeric(h)
    
    # Get transition probabilities from current node
    probs <- P[k, ]
    
    # Sample next state based on P
    k <- sample(1:length(probs), 1, prob = probs)
    
    path <- c(path, k)
  }
  
  return(list(boundary_point = k, F_w = F_w, path = path))
}

##### Section 3 - Define Estimate U(A) Function #####
estimate_u_1d <- function(start_k, P, Nx, f, h, g, N = 400) {
  "
  Estimates the potential for each given state.
  
  Parameters:
    start_k : int
      Starting state for the chain
    P : matrix
      Transition matrix
    Nx : int
      Size of lattice in x direction.
    f : function
      The Source Function. For Leplace f(k) = 0
    h : float
      Real space distance between lattice points
    g : function
      Boundary condition function
    N : int
      Number of random walks
      
  Returns:
    [u_est, ga_pa, F_avg, boundary_hits] : list
      u : float
        Estimate for the potential of the state
      ga_pa : float
        Sum of g(a) * P_A(a)
      F_avg : float
        Average accumulated contribution of the source term along a walk
      path_data : array
        Array of boundary states that the simulated walks ended at.
  
  "
  boundary_hits <- numeric(N)  # Store end points
  F_vals <- numeric(N)         # Store F(w) for each walk
  
  for (i in 1:N) {
    res <- simulate_chain_1d(start_k, P, Nx, f, h)
    boundary_hits[i] <- res$boundary_point
    F_vals[i] <- res$F_w
  }
  
  # Estimate P_A(a) — empirical distribution of boundary hits
  boundary_counts <- table(boundary_hits)
  P_A <- boundary_counts / N  # P_A(a)
  
  # Compute sum of g(a) * P_A(a)
  ga_pa <- 0
  for (a_str in names(P_A)) {
    a <- as.integer(a_str)
    ga_pa <- ga_pa + g(a) * P_A[a_str]
  }
  
  # Compute average of F(w)
  F_avg <- mean(F_vals)
  
  # Final estimate of u(A)
  u_est <- ga_pa - F_avg
  
  return(list(u = u_est, ga_pa = ga_pa, F_avg = F_avg, path_data = boundary_hits))
}

##### Section 4 - Define Solve Poisson MC in 1D Function #####
solve_poisson_mc_1d <- function(Nx, P, f, h, g, N = 400) {
  "
  Solves the poisson equation with 1-dimmension Markov Chain Monte Carlo approach
  
  Parameters:
    Nx : int
      Size of lattice in x direction.
    P : matrix
      Transition matrix
    f : function
      The Source Function. For Leplace f(k) = 0
    h : float
      Real space distance between lattice points
    g : function
      Boundary condition function
    N : int
      Number of random walks
      
    Returns:
      u_line : array
        Estimated potential in a horizontal line
  "
  
  u_line <- numeric(Nx)
  for (i in 1:Nx) {
    k <- i
    
    # If boundary point, assign directly from g
    if (i == 1 || i == Nx) {
      u_line[i] <- g(k)
    } else {
      res <- estimate_u_1d(start_k = k, P, Nx, f, h, g, N)
      u_line[i] <- res$u
    }
  }
  return(u_line)
}


##### Section 5 - Problem Setup ######
"
A capacitor is constructed out of two infinite metal plates
spaced 10 cm apart from one another. One plate is placed at 5V
potential with respect to the other and the space between the
plates is free of charges. The potential is measured in one dimension only.
"

charge_density = 1e-14
epsilon_naught = 8.854e-12

# Define our f and g function for our capacitor
f_capacitor <- function(k) {
  return(0)                   # For leplace f(i,j) = 0
}  

g_capacitor <- function(k) {
  if (k == Nx) return(5)      # right edge
  return(0)                   # left edge
}

# Define true potential values
u_true <- function(x) {
  5 * (x / 0.1)
}

##### Section 6 - Plot 1: 1-D Capacitor (Leplace Equation) #####
"
The number of random walks is taken to be 400 and the number of lattice points is 
taken to be 15. Then we find the solution for the potential and plot the 
potential vs the distance from negative plate.
"
# Define the dimensions of our transition_matrix. Our problem is in 1D so (15, 1)
Nx <- 15

# Define number of walks as array so we can investigate how the error depends on N later
N <- 400

# Distance between matrix points in real space
h <- 0.1 / (Nx-1) 

# Create transition_matrix P
P <- create_transition_matrix_1d(Nx)

# Solve our poisson equation for N = 400
u_capacitor_mc <- solve_poisson_mc_1d(Nx, P, f_capacitor, h, g_capacitor, N = N)

print(u_capacitor_mc)

x_vals <- seq(0, 0.1, length.out = Nx) 
u_capacitor_true <- u_true(x_vals)

error <- abs(u_capacitor_mc - u_capacitor_true)

MAE <- mean(error)
RMSE <- sqrt(mean((u_capacitor_mc - u_capacitor_true)^2))

print(MAE)
print(RMSE)
print(max(error))

# Define data frame for plotting
df <- data.frame(
  x = (0:(Nx - 1)) * h,
  u = u_capacitor_mc # Index 21 corresponding to N=400
)

ggplot(df, aes(x = x, y = u)) +
  # Plot our Calculated voltage in blue
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue") +
  
  # Plot actual voltage
  geom_abline(slope = 5/0.1, intercept = 0, color = "red", linetype = "dashed") +
  
  # Graph Styling
  labs(
    title = "Monte Carlo Laplace — 1D Capacitor",
    x = "Distance from negative plate (m)",
    y = "Potential (Volts)"
  ) +
  theme_minimal()

##### Section 7 - Plot 2: Total Absolute Error vs Number of Random Walks #####
"
The number of lattice points is taken to be 15 and we vary the number of random 
walks from 20 to 1000 in incremental steps of 20. We then find a solution for 
the potential for each number of random walks finding the total absolute error 
and plot the total absolute error vs the number of random walks.
"
# Define the dimensions of our transition_matrix. Our problem is in 1D so (15, 1)
Nx <- 15

# Define number of walks as array so we can investigate how the error depends on N later
N <- c(10, seq(from = 50, to = 1000, by = 50))

# Distance between matrix points in real space
h <- 0.1 / (Nx-1) 

# Solve our poisson equation in 1D to find the potential for different N
u_capacitor_mc <- matrix(nrow = 0, ncol = Nx)

for (i in 1:length(N)) {
  cat("Approximating U with", N[i], "steps \n")
  u_capacitor_mc <- rbind(u_capacitor_mc, solve_poisson_mc_1d(Nx, P, f_capacitor, h, g_capacitor, N = N[i]))
}

print(u_capacitor_mc)

x_vals <- seq(0, 0.1, length.out = Nx)
u_capacitor_true <- u_true(x_vals)

u_capacitor_true <- matrix(u_capacitor_true, nrow = nrow(u_capacitor_mc), ncol = Nx, byrow = TRUE)

error <- matrix(nrow = 0, ncol = Nx)
AE <- c() 

error <- abs(u_capacitor_mc - u_capacitor_true)

for (i in 1:length(N)) {
  AE <- c(AE, sum(error[i,]))
}

# Plot the error
plot(N, AE, type = "o", col = "red", pch = 19, lwd = 2,
     xlab = "Number of Random Walks",
     ylab = "Total Absolute Error (Volts)",
     main = "Total Absolute Error vs Number of Random Walks")

grid()

##### Section 9 - Plot 3: Total Absolute Error vs Number of Lattice Points #####
"
The number of random walks N = 400 and we vary the number of lattice points from 
2 to 30. We then find a solution for the potential for each number of lattice 
points finding the total absolute error for each. We then plot the total absolute
error vs the number of lattic points.
"

# Define the dimensions of our transition_matrix. Our problem is in 1D so (15, 1)
Nx_array <- c(2, seq(from = 4, to = 30, by = 2))

# Define number of walks as array so we can investigate how the error depends on N later
N <- 400

AE <- c() 
for (i in 1:length(Nx_array)) {
  Nx <- Nx_array[i]
  cat("Approximating U with", Nx, "lattice points \n")
  
  # Distance between matrix points in real space
  h <- 0.1 / (Nx-1) 
  
  x_vals <- seq(0, 0.1, length.out = Nx)
  u_capacitor_true <- u_true(x_vals)
  
  # Create transition_matrix P
  P <- create_transition_matrix_1d(Nx)
  
  # Solve poisson equation
  u_capacitor_mc <- solve_poisson_mc_1d(Nx, P, f_capacitor, h, g_capacitor, N = N)
  
  error <- abs(u_capacitor_mc - u_capacitor_true)
  AE <- c(AE, sum(error))
}

# Plot the error
plot(Nx_array, AE, type = "o", col = "red", pch = 19, lwd = 2,
     xlab = "Number of Lattice Points",
     ylab = "Total Absolute Error (Volts)",
     main = "Total Absolute Error vs Number of Lattice Points")

grid()

