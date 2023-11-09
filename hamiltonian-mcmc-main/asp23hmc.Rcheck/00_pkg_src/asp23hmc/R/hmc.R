hmc <-  function(location,scale = ~1, data = environment(location), light = FALSE, maxit = 100,
                num_samples = 1000, verbose = TRUE, threshold = 1/10000, include_warmup = FALSE,
                dual_average = list(DELTA = 0.65, LAMBDA = 5, KAPPA = 0.75, GAMMA = 0.05, t0 = 10),
                num_adapt = 500, num_warmup = 500, max_L = 250, reltol = sqrt(.Machine$double.eps)){

  #setup m structure

  m <- setup(location, scale, data, light, match.call())



  if(num_adapt > num_warmup){
    message('The maximum number of tuning iterations is longer than the warmup phase.\nSetting warmup length same length as maximum number of warmup phase.')
    num_warmup <- max(num_adapt,num_warmup)
  }

  # Store the iteration number at which the warmup phase concludes for plotting purposes.
  warmup_over <- FALSE

  # Access the parameters from the dual_average list, use default value if incomplete list was given as argument
  DELTA <- ifelse(is.null(dual_average$DELTA),0.65, dual_average$DELTA)
  LAMBDA <- ifelse(is.null(dual_average$LAMBDA),5, dual_average$LAMBDA)
  KAPPA <- ifelse(is.null(dual_average$KAPPA),0.75, dual_average$KAPPA)
  GAMMA <- ifelse(is.null(dual_average$GAMMA),0.05, dual_average$GAMMA)
  t0 <- ifelse(is.null(dual_average$t0),10, dual_average$t0)

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
  if(attributes(m)$class[1] != "hmc"){
    stop("Cannot run HMC, provided R object is not of class 'lmls'")
  }

  if(m$light){
    stop("Cannot run MCMC, lmls() called with argument 'light = TRUE'")
  }

  initial_params <- ML_estimates(m$y, m$x, m$z)

  if(verbose) paste("Initial values: ", initial_params)

  # Define length of parameter vectors
  n_beta <- ncol(m$x) #length(coef(m, "location"))
  n_gamma <- ncol(m$z) #length(coef(m, "scale"))    # or

  m$hmc <-list(
    location = matrix(nrow = num_samples+1,
                      ncol = ncol(m$x),
                      dimnames = list(NULL, colnames(m$x))),
    scale = matrix(nrow = num_samples+1,
                   ncol = ncol(m$z),
                   dimnames = list(NULL, colnames(m$z))),
    accepted_beta = vector(mode = 'integer', length = num_samples+1),
    accepted_gamma = vector(mode = 'integer', length = num_samples+1),
    rolling_average_location <- matrix(nrow = num_samples+1, ncol = n_beta),
    rolling_average_scale <- matrix(nrow = num_samples+1, ncol = n_gamma),
    num_samples = num_samples
    )

  # Store warmup iterations for plotting purposes inside m
  m$hmc$end_warmup <- num_warmup

  # Set starting values as first iteration value in sample matrix
  m$hmc$location[1,] <- head(initial_params, n_beta)
  m$hmc$scale[1,] <- tail(initial_params, n_gamma)

  m$hmc$rolling_average_location_without_warmup <- matrix(c(NA), nrow = num_samples+1, ncol = n_beta, byrow = TRUE)
  m$hmc$rolling_average_scale_without_warmup <- matrix(c(NA), nrow = num_samples+1, ncol = n_gamma, byrow = TRUE)

  # Set initial step sizes
  eps_beta <- c(find_reasonable_epsilon(m,grad_U_beta, current_q = m$hmc$location[1,], given_other = m$hmc$scale[1,], U, K, epsilon = 1, verbose = TRUE, flag = TRUE))
  eps_gamma <- c(find_reasonable_epsilon(m, grad_U_gamma, current_q = m$hmc$scale[1,], given_other = m$hmc$location[1,], U, K, epsilon = 1, verbose = TRUE, flag = FALSE))

  # Constants for Dual Averaging
  MU_beta <- log(t0*eps_beta[1])
  eps_beta_bar <- c(1)
  H_beta_bar <- c(0)

  MU_gamma <- log(t0*eps_gamma[1])
  eps_gamma_bar <- c(1)
  H_gamma_bar <- c(0)

  # Calculate the MCMC whilst adapting the step sizes for the adaptive period

  eps_beta_delta = 1
  eps_gamma_delta = 1

  # In R, vector indexing starts from 1 instead of 0.
  # To align with vector indexing set 'i' equal to 'iteration_index + 1'.

  cli_progress_bar("Parameter tuning: ")

  for(i in 2:(num_adapt+1))  {
    # Generate beta samples and save acceptance probability
    result <- HMC_Dual_Averaging(m, grad_U_beta, current_q = m$hmc$location[i-1,], given_other = m$hmc$scale[i-1,], epsilon = eps_beta[i-1], L = min(max_L,max(1, round(LAMBDA/eps_beta[i-1]))), flag = TRUE)
    m$hmc$location[i,] <- result$q
    m$hmc$accepted_beta[i] <- result$accepted
    alpha_beta <- result$alpha

    # Generate gamma samples
    result <- HMC_Dual_Averaging(m, grad_U_gamma, current_q = m$hmc$scale[i-1,], given_other = m$hmc$location[i,],  epsilon = eps_gamma[i-1], L = min(max_L,max(1, round(LAMBDA/eps_gamma[i-1]))), flag = FALSE)
    m$hmc$scale[i,] <- result$q
    m$hmc$accepted_gamma[i] <- result$accepted
    alpha_gamma <- result$alpha

    # Adapt step size for location
    H_beta_bar[i] <- (1 - 1/((i-1) + t0))*H_beta_bar[i-1] + (1/((i-1) + t0))*(DELTA - alpha_beta)
    log_eps <- MU_beta - (sqrt(i-1)/GAMMA)*H_beta_bar[i]
    eps_beta[i] <- exp(log_eps)
    eps_beta_bar[i] <- exp(((i-1)^(-KAPPA))*log_eps + (1 - (i-1)^(-KAPPA))*log(eps_beta_bar[i-1]))

    # Adapt step size for scale !!!!!!! ANPASSEN!
    H_gamma_bar[i] <- (1 - 1/((i-1) + t0))*H_gamma_bar[i-1] + (1/((i-1) + t0))*(DELTA - alpha_gamma)
    log_eps <- MU_gamma - (sqrt(i-1)/GAMMA)*H_gamma_bar[i]
    eps_gamma[i] <- exp(log_eps)
    eps_gamma_bar[i] <- exp(((i-1)^(-KAPPA))*log_eps + (1 - (i-1)^(-KAPPA))*log(eps_gamma_bar[i-1]))

    # Calculate delta for parameter tuning
    eps_beta_delta <- calculate_delta(eps_beta_bar)
    eps_gamma_delta <- calculate_delta(eps_gamma_bar)

    # Compute new rolling average for the chain
    m$hmc$rolling_average_location[i,] <- colMeans(m$hmc$location[1:i, ], na.rm = TRUE)
    m$hmc$rolling_average_scale[i,] <- colMeans(m$hmc$scale[1:i, ], na.rm = TRUE)

    cli_progress_update()

    if ((eps_beta_delta < threshold && eps_gamma_delta < threshold)){
      cli_progress_done()
      break
    }
  }
  # Save the number of adaptive iterations completed.
  m$hmc$total_adapt = i - 1

  # Step size for the remaining non-adapting HMC steps
  epsilon_beta <- eps_beta_bar[i]
  epsilon_gamma <- eps_gamma_bar[i]


  if(verbose) message("Finished the step size adaptations after ", i-1," iterations.\nFor beta: epsilon = ", epsilon_beta, "\nFor gamma: epsilon = ", epsilon_gamma)

  if(verbose) message("For beta: suggested L = ", round(LAMBDA/epsilon_beta), "\nFor gamma: suggested L = ", round(LAMBDA/epsilon_gamma))

  # Number of steps for the remaining non-adapting HMC steps suggested per Chain
  L_beta <- min(max(1, round(LAMBDA/epsilon_beta)), max_L)
  L_gamma <- min(max(1, round(LAMBDA/epsilon_gamma)), max_L)

  if(verbose) message("max_L = ", max_L, "\nFor beta: L = ", L_beta, "\nFor gamma: L = ", L_gamma)

  cli_progress_bar("Sampling: ")

  # Perform gibbs updates by sequentially updating beta then gamma
  for (i in (i+1):(num_samples+1)) {
    # Check if warmup phase is over
    if(i == num_warmup + 2) warmup_over <- TRUE
    cli_progress_update()

    # Generate beta samples
    result <- HMC_Dual_Averaging(m,grad_U_beta, current_q = m$hmc$location[i-1,], given_other = m$hmc$scale[i-1,], epsilon = epsilon_beta, L_beta, flag = TRUE)
    m$hmc$location[i,] <- result$q
    m$hmc$accepted_beta[i] <- result$accepted

    # Generate gamma samples
    result <- HMC_Dual_Averaging(m,grad_U_gamma, current_q = m$hmc$scale[i-1,], given_other = m$hmc$location[i,],  epsilon = epsilon_gamma, L_gamma, flag = FALSE)
    m$hmc$scale[i,] <- result$q
    m$hmc$accepted_gamma[i] <- result$accepted

    # Compute new rolling average for the chain
    m$hmc$rolling_average_location[i,] <- colMeans(m$hmc$location[1:i, ], na.rm = TRUE)
    m$hmc$rolling_average_scale[i,] <- colMeans(m$hmc$scale[1:i, ], na.rm = TRUE)

    if(warmup_over){
      # Calculate rolling averages catch the case were there only is 1 value
      m$hmc$rolling_average_location_without_warmup[i,] <- colMeans(m$hmc$location[(num_warmup+2):i, ,drop = FALSE])
      m$hmc$rolling_average_scale_without_warmup[i,] <- colMeans(m$hmc$scale[(num_warmup+2):i, ,drop = FALSE ])
    }
  }
  if(include_warmup){
    m$coefficients$location[1:n_beta]<- tail(m$hmc$rolling_average_location,1)
    m$coefficients$scale[1:n_gamma] <- tail(m$hmc$rolling_average_scale,1)

  } else {
    m$coefficients$location[1:n_beta]<- tail(m$hmc$rolling_average_location_without_warmup,1)
    m$coefficients$scale[1:n_gamma] <- tail(m$hmc$rolling_average_scale_without_warmup,1)
  }

  # Calculate acceptance rates
  m$hmc$acceptance_rate_beta_adapt_phase <- sum(m$hmc$accepted_beta[1:num_adapt+1])/(num_adapt)
  m$hmc$acceptance_rate_gamma_adapt_phase <- sum(m$hmc$accepted_gamma[1:num_adapt+1])/(num_adapt)

  m$hmc$acceptance_rate_beta_tuned <- sum(m$hmc$accepted_beta[(num_adapt+1):(num_samples)+1])/(num_samples-num_adapt)
  m$hmc$acceptance_rate_gamma_tuned <- sum(m$hmc$accepted_gamma[(num_adapt+1):(num_samples)+1])/(num_samples-num_adapt)

  m$hmc$acceptance_rate_beta <- sum(m$hmc$accepted_beta[1:num_samples+1])/(num_samples)
  m$hmc$acceptance_rate_gamma <- sum(m$hmc$accepted_gamma[1:num_samples+1])/(num_samples)

  m$acceptance_ratios <- acceptance_ratios.hmc(m)

  return(m)
}
