set.seed(963)

##### Section 1 - Define Functions #####
f=function(x, mu, a, k) {
  ifelse(x >= mu, (a * ((k - mu + x) / k)^(-a - 1)) / k, 0)
}

# CDF
F <- function(x, mu, a, k) {
  ifelse(x < mu, 0, 1 - ( (k - mu + x)/k )^(-a))
}

# Inverse CDF
F_inv=function(u, mu, a, k) { 
  return(k * (1 - u)^(-1 / a) - k + mu)
}


g <- function(x, lambda, mu) {
  ifelse(x >= mu, lambda * exp(-lambda * (x - mu)), 0)
}

##### Section 2 - Set Parameters #####
"Set parameters"
a <- 3
k <- 2
mu <- 1
sample_sizes <- c(100, 500, 1000, 5000, 7500, 10000)

##### Section 3 - Inverse Transform Method #####
"We are going to use the inverse transform method"
inverse_transform_method <- function(n, a, k, mu) {
  # Generate Uniform Samples
  u <- runif(n)
  samples <- F_inv(u, mu, a, k)
  
  # Plot histogram with theoretical curve overlay
  hist(samples, breaks = 50, probability = TRUE, main = bquote("Generated Samples + PDF  ( n =" ~ .(n) ~ ")"), col = "green")
 
  # Theoretical PDF for comparison
  curve(f(x, mu, a, k), from = mu, to = max(samples), add = TRUE, col = "red", lwd = 2)
   
  # Make your comparison quantitative by using a goodness-of-fit method
  ks.test(samples, function(x) F(x, mu, a, k))
}

statistics <- sapply(sample_sizes, function(n) inverse_transform_method(n, a, k, mu))
print(statistics)

##### Section 4 - Importance Sampling #####
importance_sampling <- function(n, a, k, mu) {
  # Importance Sampling Estimator of E[|X|]
  
  lambda = 0.5
  samples <- rexp(n, lambda) + mu
  
  # Calculate weights
  weights <- f(samples, a, k, mu) / g(samples, lambda, mu)
  
  # Estimate E[|X|]
  estimate <- sum(abs(samples) * weights) / sum(weights)
  
  return(estimate)
}

# Run for various sample sizes
estimates <- sapply(sample_sizes, function(n) importance_sampling(n, a, k, mu))
print(estimates)

# Plot convergence
plot(sample_sizes, estimates, type = "b", pch = 19, col = "red",
     xlab = "Sample Size (n)", ylab = "Estimate of E[|X|]",
     main = "Convergence of Importance Sampling Estimate")
grid()


