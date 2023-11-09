# hmc-setup

# Function to set up a hierarchical modeling framework for Hamiltonian Monte Carlo (HMC)
setup <- function(location, scale = ~ 1, data = environment(location), light = TRUE, call = NULL) {

  # Update the 'scale' formula to include the model's response variable.
  scale <- update(scale, paste(location[[2]], "~."))

  # Define differnt metrics for the initial structure
  y <- eval(location[[2]], data, environment(location))
  x <- model.matrix(location, data)
  z <- model.matrix(scale, data)
  nobs <- length(y)
  df <- ncol(x) + ncol(z)
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


# hmc function -----------------------------------------------------------------


#' Hamiltonian Markov Chain Monte Carlo (HMC) Sampling
#'
#' This function performs HMC sampling for a given location-scale regression model.
#'
#' @param location A formula specifying the model's location component (e.g., y ~ x).
#' @param scale A formula specifying the model's scale component (e.g., ~ 1).
#' @param data The dataset containing the model's variables.
#' @param light A logical indicating whether to use the "light" version of the model (default is FALSE).
#' @param num_samples The number of MCMC samples to generate (default is 1000).
#' @param verbose A logical indicating whether to display progress messages (default is TRUE).
#' @param threshold The convergence threshold for adaptive step size tuning (default is 1/10000).
#' @param include_warmup A logical indicating whether to include warm-up samples in the final results (default is FALSE).
#' @param dual_average A list of parameters for dual averaging step size adaptation.
#' @param num_adapt The maximum number of adaptation iterations (default is 500).
#' @param num_warmup The number of warm-up iterations (default is 500).
#' @param max_L The maximum number of leapfrog steps for HMC (default is 250).
#'
#' @return An HMC model object containing MCMC samples and results.
#'
#' @export
#'
#' @examples
#' library(gamlss.data)
#' data(abdom, package = "gamlss.data")
#' output <- hmc(y ~ x, ~x, data = abdom)
#'
hmc <- function(location = y ~ x, scale = ~ 1, data = environment(location), light = FALSE,
                num_samples = 1000, verbose = TRUE, threshold = 1/10000, include_warmup = FALSE,
                dual_average = list(DELTA = 0.65, LAMBDA = 5, KAPPA = 0.75, GAMMA = 0.05, t0 = 10),
                num_adapt = 500, num_warmup = 500, max_L = 250) {

  # Setup the model structure
  m <- setup(location, scale, data, light, match.call())

  # Check if num_adapt is longer than num_warmup and adjust if needed
  if (num_adapt > num_warmup) {
    message('The maximum number of tuning iterations is longer than the warmup phase.\nSetting warmup length same length as maximum number of warmup phase.')
    num_warmup <- max(num_adapt, num_warmup)
  }

  # Store the iteration number at which the warm-up phase concludes for plotting
  warmup_over <- FALSE

  # Access dual_average parameters and provide default values if incomplete
  DELTA <- ifelse(is.null(dual_average$DELTA), 0.65, dual_average$DELTA)
  LAMBDA <- ifelse(is.null(dual_average$LAMBDA), 5, dual_average$LAMBDA)
  KAPPA <- ifelse(is.null(dual_average$KAPPA), 0.75, dual_average$KAPPA)
  GAMMA <- ifelse(is.null(dual_average$GAMMA), 0.05, dual_average$GAMMA)
  t0 <- ifelse(is.null(dual_average$t0), 10, dual_average$t0)

  # Input validation checks
  if (DELTA <= 0 || DELTA > 1) {
    stop("DELTA must be in (0, 1]!")
  }
  if (LAMBDA <= 0) {
    stop("LAMBDA must be greater than 0!")
  }
  if (KAPPA <= 0 || KAPPA > 1) {
    stop("KAPPA must be in (0, 1]!")
  }
  if (GAMMA <= 0) {
    stop("GAMMA must be greater than 0!")
  }
  if (attributes(m)$class[1] != "hmc") {
    stop("Cannot run HMC, provided R object is not of class 'lmls'")
  }

  # Calculate initial parameter estimates
  initial_params <- ML_estimates(m$y, m$x, m$z)

  if (verbose) paste("Initial values: ", initial_params)

  # Define length of parameter vectors
  n_beta <- ncol(m$x)
  n_gamma <- ncol(m$z)

  # Create a structure to store MCMC results
  m$hmc <- list(
    location = matrix(
      nrow = num_samples + 1,
      ncol = ncol(m$x),
      dimnames = list(NULL, colnames(m$x))
    ),
    scale = matrix(
      nrow = num_samples + 1,
      ncol = ncol(m$z),
      dimnames = list(NULL, colnames(m$z))
    ),
    accepted_beta = vector(mode = 'integer', length = num_samples + 1),
    accepted_gamma = vector(mode = 'integer', length = num_samples + 1),
    rolling_average_location = matrix(nrow = num_samples + 1, ncol = n_beta),
    num_samples = num_samples,
    end_warmup = num_warmup,
    total_adapt = num_adapt,
    rolling_average_scale = matrix(nrow = num_samples + 1, ncol = n_gamma)
  )

  # Set starting values as the first iteration value in the sample matrix
  m$hmc$location[1,] <- head(initial_params, n_beta)
  m$hmc$scale[1,] <- tail(initial_params, n_gamma)

  m$hmc$rolling_average_location_without_warmup <- matrix(c(NA), nrow = num_samples + 1, ncol = n_beta, byrow = TRUE)
  m$hmc$rolling_average_scale_without_warmup <- matrix(c(NA), nrow = num_samples + 1, ncol = n_gamma, byrow = TRUE)

  # Set initial step sizes
  eps_beta <- c(find_reasonable_epsilon(m, grad_U_beta, current_q = m$hmc$location[1,], given_other = m$hmc$scale[1,], U, K, epsilon = 1, verbose = TRUE, flag = TRUE))
  eps_gamma <- c(find_reasonable_epsilon(m, grad_U_gamma, current_q = m$hmc$scale[1,], given_other = m$hmc$location[1,], U, K, epsilon = 1, verbose = TRUE, flag = FALSE))

  # Constants for Dual Averaging
  MU_beta <- log(t0 * eps_beta[1])
  eps_beta_bar <- c(1)
  H_beta_bar <- c(0)

  MU_gamma <- log(t0 * eps_gamma[1])
  eps_gamma_bar <- c(1)
  H_gamma_bar <- c(0)

  # Calculate the MCMC while adapting the step sizes for the adaptive period
  eps_beta_delta = 1
  eps_gamma_delta = 1

  # Initialize a progress bar
  cli_progress_bar("Parameter tuning: ")

  for (i in 2:(num_adapt + 1)) {
    # Generate beta samples and save acceptance probability
    result <- HMC_Dual_Averaging(m, grad_U_beta, current_q = m$hmc$location[i - 1,], given_other = m$hmc$scale[i - 1,], epsilon = eps_beta[i - 1], L = min(max_L, max(1, round(LAMBDA / eps_beta[i - 1]))), flag = TRUE)
    m$hmc$location[i,] <- result$q
    m$hmc$accepted_beta[i] <- result$accepted
    alpha_beta <- result$alpha

    # Generate gamma samples
    result <- HMC_Dual_Averaging(m, grad_U_gamma, current_q = m$hmc$scale[i - 1,], given_other = m$hmc$location[i,], epsilon = eps_gamma[i - 1], L = min(max_L, max(1, round(LAMBDA / eps_gamma[i - 1]))), flag = FALSE)
    m$hmc$scale[i,] <- result$q
    m$hmc$accepted_gamma[i] <- result$accepted
    alpha_gamma <- result$alpha

    # Adapt step size for location
    H_beta_bar[i] <- (1 - 1 / ((i - 1) + t0)) * H_beta_bar[i - 1] + (1 / ((i - 1) + t0)) * (DELTA - alpha_beta)
    log_eps <- MU_beta - (sqrt(i - 1) / GAMMA) * H_beta_bar[i]
    eps_beta[i] <- exp(log_eps)
    eps_beta_bar[i] <- exp(((i - 1) ^ (-KAPPA)) * log_eps + (1 - (i - 1) ^ (-KAPPA)) * log(eps_beta_bar[i - 1]))

    # Adapt step size for scale
    H_gamma_bar[i] <- (1 - 1 / ((i - 1) + t0)) * H_gamma_bar[i - 1] + (1 / ((i - 1) + t0)) * (DELTA - alpha_gamma)
    log_eps <- MU_gamma - (sqrt(i - 1) / GAMMA) * H_gamma_bar[i]
    eps_gamma[i] <- exp(log_eps)
    eps_gamma_bar[i] <- exp(((i - 1) ^ (-KAPPA)) * log_eps + (1 - (i - 1) ^ (-KAPPA)) * log(eps_gamma_bar[i - 1]))

    # Calculate delta for parameter tuning
    eps_beta_delta <- calculate_delta(eps_beta_bar)
    eps_gamma_delta <- calculate_delta(eps_gamma_bar)

    # Compute new rolling average for the chain
    m$hmc$rolling_average_location[i,] <- colMeans(m$hmc$location[1:i,], na.rm = TRUE)
    m$hmc$rolling_average_scale[i,] <- colMeans(m$hmc$scale[1:i,], na.rm = TRUE)

    # Update the progress bar
    cli_progress_update()

    # Check for convergence
    if ((eps_beta_delta < threshold && eps_gamma_delta < threshold)) {
      cli_progress_done()
      break
    }
  }

  # Save the number of adaptive iterations completed
  m$hmc$total_adapt = i - 1

  # Step size for the remaining non-adapting HMC steps
  epsilon_beta <- eps_beta_bar[i]
  epsilon_gamma <- eps_gamma_bar[i]

  if (verbose) message("Finished the step size adaptations after ", i - 1, " iterations.\nFor beta: epsilon = ", epsilon_beta, "\nFor gamma: epsilon = ", epsilon_gamma)
  if (verbose) message("For beta: suggested L = ", round(LAMBDA / epsilon_beta), "\nFor gamma: suggested L = ", round(LAMBDA / epsilon_gamma))

  # Number of steps for the remaining non-adapting HMC steps suggested per Chain
  L_beta <- min(max(1, round(LAMBDA / epsilon_beta)), max_L)
  L_gamma <- min(max(1, round(LAMBDA / epsilon_gamma)), max_L)

  if (verbose) message("max_L = ", max_L, "\nFor beta: L = ", L_beta, "\nFor gamma: L = ", L_gamma)

  # Initialize a progress bar for sampling
  cli_progress_bar("Sampling: ")

  # Perform Gibbs updates by sequentially updating beta then gamma
  for (i in (i + 1):(num_samples + 1)) {
    # Check if warm-up phase is over
    if (i == num_warmup + 2) warmup_over <- TRUE
    cli_progress_update()

    # Generate beta samples
    result <- HMC_Dual_Averaging(m, grad_U_beta, current_q = m$hmc$location[i - 1,], given_other = m$hmc$scale[i - 1,], epsilon = epsilon_beta, L_beta, flag = TRUE)
    m$hmc$location[i,] <- result$q
    m$hmc$accepted_beta[i] <- result$accepted

    # Generate gamma samples
    result <- HMC_Dual_Averaging(m, grad_U_gamma, current_q = m$hmc$scale[i - 1,], given_other = m$hmc$location[i,], epsilon = epsilon_gamma, L_gamma, flag = FALSE)
    m$hmc$scale[i,] <- result$q
    m$hmc$accepted_gamma[i] <- result$accepted

    # Compute new rolling average for the chain
    m$hmc$rolling_average_location[i,] <- colMeans(m$hmc$location[1:i,], na.rm = TRUE)
    m$hmc$rolling_average_scale[i,] <- colMeans(m$hmc$scale[1:i,], na.rm = TRUE)

    if (warmup_over) {
      # Calculate rolling averages, catch the case where there's only 1 value
      m$hmc$rolling_average_location_without_warmup[i,] <- colMeans(m$hmc$location[(num_warmup + 2):i,, drop = FALSE])
      m$hmc$rolling_average_scale_without_warmup[i,] <- colMeans(m$hmc$scale[(num_warmup + 2):i,, drop = FALSE ])
    }
  }

  # Use warm-up or non-warm-up rolling averages based on the include_warmup flag
  if (include_warmup) {
    m$coefficients$location[1:n_beta] <- tail(m$hmc$rolling_average_location, 1)
    m$coefficients$scale[1:n_gamma] <- tail(m$hmc$rolling_average_scale, 1)
  } else {
    m$coefficients$location[1:n_beta] <- tail(m$hmc$rolling_average_location_without_warmup, 1)
    m$coefficients$scale[1:n_gamma] <- tail(m$hmc$rolling_average_scale_without_warmup, 1)
  }

  # Name the coefficient vectors based on column names in the data
  names(m$coefficients$location) <- colnames(m$x)
  names(m$coefficients$scale) <- colnames(m$z)

  # Store acceptance ratios and hyperparameters in the model object
  m$acceptance_ratios <- acceptance_ratios(m, num_adapt, num_samples)
  m$hyperparameter <- hyperparameter(epsilon_beta, epsilon_gamma, L_beta, L_gamma, LAMBDA)

  # Return the model object
  return(m)
}





