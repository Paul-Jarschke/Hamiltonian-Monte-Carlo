# ML-estimates -----------------------------------------------------------------

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

# log-likelihood --------------------------------------------------------------

#' @importFrom stats dnorm fitted

loglik <- function(m) {
  y <- m$y
  location <- m$x %*% m$coefficients$location
  scale <- exp(m$z %*% m$coefficients$scale)

  # Sum the logarithm of the probability density function (pdf) of the observed data.
  # A higher log-likelihood indicates a better fit.
  sum(dnorm(y, location, scale, log = TRUE))
}

#' @importFrom stats dnorm fitted

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





# Plot helpers -----------------------------------------------------------------
plot_trace <- function(m, num_coef, var_names, data){
  for (i in 1:num_coef) {
    x_values <- 0:(length(data[,i]) - 1)  # Create x-values starting from 0
    plot(x_values, data[,i], type = 'l', xlab = "Iteration", ylab = 'value', main = var_names[i],
         xlim = c(0, length(data[,i]) - 1))  # Set x-axis limits to start at 0
    abline(v = m$hmc$total_adapt, col = "red", lty = 2)  # Add a vertical red line at num_adapt
    abline(v = m$hmc$end_warmup, col = "blue", lty = 3)  # Add a vertical blue dotted line at end_warm_up
    abline(h = data[1, i], col = "black", lty = 3)  # Add a horizontal dashed line for initialization
  }

  legend("bottomright", bty = "n", inset=c(0,0),
         legend = c("End Adaptation", "End Warmup", "Initialization"),
         col = c("red", "blue", "black"), lty = c(2, 2, 3), cex = 0.6)
}



plot_hist <- function(m, num_coef, var_names, data){
  for (i in 1:num_coef) {
    hist(data[, i], xlab = var_names[i], ylab = 'Frequency', main = var_names[i])
  }
}



plot_density <- function(m, predictor, num_coef, var_names, data){
  for (i in 1:num_coef){
    plot(density(data[, i]), main = var_names[i])
    # Add a vertical line for the mean
    abline(v = m$coefficients[[predictor]][i], col = "red")}
}



plot_acf <- function(num_coef, var_names, data){
  for (i in 1:num_coef) {
    acf(data[,i], main = var_names[i])
  }
}


plot_roll_mean <- function(m, predictor, num_coef, var_names, data){
  for(i in 1:num_coef){
    # Compute variables based on predictor
    if(predictor == 'location'){
      rolling_average <- m$hmc$rolling_average_location[-1, i, drop = FALSE]
      rolling_average_without_warm_up <- m$hmc$rolling_average_location_without_warmup[-1, i, drop = FALSE]

      # Combine both rolling averages and initial estimation values to determine y-axis limits
      combined_data <- c(rolling_average, rolling_average_without_warm_up, data[1, i])
    }else{
      rolling_average <- m$hmc$rolling_average_scale[-1, i, drop = FALSE]
      rolling_average_without_warm_up <- m$hmc$rolling_average_scale_without_warmup[-1, i, drop = FALSE]

      # Combine both rolling averages and initial estimation values to determine y-axis limits
      combined_data <- c(rolling_average, rolling_average_without_warm_up, data[1, i])
    }

    y_min <- min(combined_data, na.rm = TRUE)
    y_max <- max(combined_data, na.rm = TRUE)

    # Set y-axis limits explicitly to ensure the horizontal dashed line is visible
    ylim <- c(y_min - 0.01 * abs(y_min), y_max + 0.01 * abs(y_max))

    matplot(1:m$hmc$num_samples, cbind(rolling_average, rolling_average_without_warm_up),
            type = "l", xlab = "Iteration", ylab = "Rolling Average",
            main = var_names[i], ylim = ylim, col = c("blue", "darkgreen"), lty=c(1))

    abline(v = m$hmc$total_adapt, col = "red", lty = 2)  # Add a vertical red line at num_adapt
    abline(v = m$hmc$end_warmup, col = "blue", lty = 3)  # Add a vertical blue dotted line at end_warm_up
    abline(h = data[1, i], col = "black", lty = 3)  # Add a horizontal dashed line for initialization
  }
  legend("bottomright", bty = "n", inset = 0.01,
         legend = c("End Adaptation", "End Warmup", "Initialization", "Rolling Average (All samples)", "Rolling Average (Without Warmup)"),
         col = c("red", "blue", "black", "blue", "darkgreen"), lty = c(2, 3, 3, 1, 1), cex = 0.6) # Create a legend
}
