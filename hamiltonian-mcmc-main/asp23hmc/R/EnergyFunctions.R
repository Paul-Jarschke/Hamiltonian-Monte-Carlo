# Function to calculate the potential energy (U) based on position vector (q)
U <- function(m, q, given_other, flag) {
  # Determine the number of columns in X and Z matrices
  cols_beta <- ncol(m$x)   # Number of columns in the X matrix
  cols_gamma <- ncol(m$z)  # Number of columns in the Z matrix

  # Extract beta and gamma from the input parameter q based on the flag
  # The flag displays if beta should be used for the q-vector (flag = TRUE).
  if (flag == TRUE) {
    beta <- q
    gamma <- given_other
  } else {
    beta <- given_other
    gamma <- q
  }

  # Calculate and return the negative log-likelihood of the data
  return(-sum(dnorm(m$y, m$x %*% beta, exp((m$z %*% gamma)), log = TRUE)))
}

# Function to calculate the kinetic energy (K) based on momentum vector (p)
K <- function(p) {
  # Calculate the squared elements of the momentum vector and sum them up
  return(sum(p ^ 2 / 2))
}

# Function to calculate the gradient of U with respect to beta
grad_U_beta <- function(m, beta, gamma) {
  # Compute the residuals and exp_term for all observations at once
  residuals <- m$y - m$x %*% beta
  exp_term <- exp(2 * (m$z %*% gamma))

  # Calculate the coefficient matrix for all observations
  coef_matrix <- residuals / exp_term

  # Compute the gradient as a matrix product
  grad <- -t(m$x) %*% coef_matrix

  return(grad)
}

# Function to calculate the gradient of U with respect to gamma
grad_U_gamma <- function(m, gamma, beta) {
  # Compute the residuals for all observations at once
  residuals <- m$y - m$x %*% beta

  # Compute the exp_term for all observations at once
  exp_term <- exp(m$z %*% gamma)

  # Calculate the squared coefficients and their squares for all observations
  coef_sq <- (residuals / exp_term)^2

  # Compute the gradient as a matrix product
  grad <- t(m$z) %*% (1 - coef_sq)

  return(grad)
}

# Calculate the gradient of the potential energy function with respect to the parameters in q.
grad_U <- function(m) {
  # This function combines gradients for beta and gamma components.
  return(c(grad_U_beta(m), grad_U_gamma(m)))
}

# Calculate the gradient of the kinetic energy function
grad_K <- function(p) {
  return(p)
}
