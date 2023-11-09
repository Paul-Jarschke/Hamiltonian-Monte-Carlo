# Calculate the log-likelihood of a model.
loglik <- function(m) {
  y <- m$y
  location <- fitted(m, "location")
  scale <- fitted(m, "scale")

  # Sum the logarithm of the probability density function (pdf) of the observed data.
  # A higher log-likelihood indicates a better fit.
  sum(dnorm(y, location, scale, log = TRUE))
}

# Calculate the negative log-likelihood of a model.
neg_loglik <- function(m) {
  y <- m$y
  location <- fitted(m, "location")
  scale <- fitted(m, "scale")

  # Sum the logarithm of the probability density function (pdf) of the observed data.
  # A lower negative log-likelihood indicates a better fit.
  -sum(dnorm(y, location, scale, log = TRUE))
}

# Calculate the proportional change between the last two elements of a vector x_bar.
calculate_delta <-  function(x_bar) {
  i <- length(x_bar)
  abs((x_bar[i] - x_bar[i-1])/x_bar[i-1])
}

find_reasonable_epsilon = function(m, FUN, current_q, given_other, U, K, epsilon = 1, verbose = TRUE, flag) {
  # Save the length of the position vector
  # Calculate the log ratio of the proposal and current state's Hamiltonian:
  # H(q_proposal, p_proposal) - H(q_current, p_current).
  # Recall that P(q, p) ∝ exp(-H(q, p) and H(q, p) = U(q) + K(p)
  log_ratio <- log(HMC_Dual_Averaging(m,FUN, current_q, given_other, epsilon, L = 1, flag)$alpha)

  # Determine whether to double or halve the step size based on the acceptance probability.
  zeta  <- ifelse(!is.nan(log_ratio) && !is.nan(exp(log_ratio)) && exp(log_ratio) > 0.5, 1, -1)

  # Keep adjusting epsilon, doubling or halving it, until stopping criterion is met
  for (count in 1:100) {
    if(is.nan(log_ratio) || zeta  * log_ratio > (-zeta ) * log(2)){
      epsilon <- 2 ** zeta  * epsilon

      # Make 1 leapfrog step with new epsilon (step size)
      # Evaluate potential and kinetic energies at start and end of trajectory
      # Calculate the log ratio of the proposal and current state's Hamiltonian:
      # H(q_proposal, p_proposal) - H(q_current, p_current).
      # Recall that P(q, p) ∝ exp(-H(q, p) and H(q, p) = U(q) + K(p),

      log_ratio <- log(HMC_Dual_Averaging(m, FUN, current_q, given_other, epsilon, L = 1, flag)$alpha)
    }
    else{
      if (verbose) message("Reasonable epsilon = ", epsilon, " found after ", count, " steps")
      return(epsilon)
    }
  }
  stop("Could not find a reasonable epsilon in 10000 iterations!")
}


# The joint log density (jld) is proportional to the Hamiltonian function H(q, p) = U(q) + K(p).
joint_log_density <- function(m,q, given_other, p, flag) -U(m,q, given_other, flag) - K(p)


# Function to estimate initial values for parameters using Maximum Likelihood (ML) estimation
ML_estimates <- function(y, X, Z){

  # Calculate initial_beta using the formula for linear regression
  initial_beta <- solve(t(X) %*% X) %*% t(X) %*% y

  # Calculate s_hat using vectorized operations
  s_hat <- log(abs(y - X %*% initial_beta)) + 0.6351814

  # Calculate initial_gamma
  initial_gamma <- solve(t(Z) %*% Z) %*% t(Z) %*% s_hat

  # Combine initial_beta and initial_gamma into a single parameter vector initial_q
  initial_q <- c(initial_beta, initial_gamma)

  # Return the vector of initial parameter estimates
  return(initial_q)
}


# Function to set up a hierarchical modeling framework for Hamiltonian Monte Carlo (HMC)
setup <- function(location, scale = ~ 1, data = environment(location), light = TRUE, call = NULL) {

  # Update the 'scale' formula to include the model's response variable.
  scale <- update(scale, paste(location[[2]], "~."))

  # Evaluate the model's response variable 'y'
  y <- eval(location[[2]], data, environment(location))

  # Create the model matrix 'x' using the 'location'
  x <- model.matrix(location, data)

  # Create the model matrix 'z' using the 'scale'
  z <- model.matrix(scale, data)

  # Calculate the number of observations 'nobs'
  nobs <- length(y)

  # Calculate the total degrees of freedom 'df'
  df <- ncol(x) + ncol(z)

  # Calculate the residual degrees of freedom 'df.residual'
  df.residual <- nobs - df

  # Create a model structure 'm' containing relevant information:
  m <- structure(
    list(
      y               = y,            # Model response variable
      x               = x,            # Model matrix for 'location'
      z               = z,            # Model matrix for 'scale'
      nobs            = nobs,         # Number of observations
      df              = df,           # Total degrees of freedom
      df.residual     = df.residual,  # Residual degrees of freedom
      light           = light,        # Flag for lightweight setup
      call            = call,         # Function call
      coefficients    = list(location = NULL, scale = NULL)  # Placeholder for coefficients
    ),
    class = "hmc"  # Assign the class "hmc" to this structure
  )
  # Return the model structure 'm'
  return(m)
}



