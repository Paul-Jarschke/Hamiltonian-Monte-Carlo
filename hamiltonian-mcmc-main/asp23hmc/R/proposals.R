HMC_Dual_Averaging <- function(m, FUN, current_q, given_other, epsilon, L, flag) {
  # Remember position from last state of Markov Chain
  q <- current_q

  # Save the length of the position vector
  n <- length(current_q)

  # Sample random momenta
  p <- rnorm(length(q), 0, 1)   # default: independent standard normal variables
  current_p <- p

  # Compute deterministic trajectory depending on method
  phase <- leapfrog(m, FUN, p, q, given_other, epsilon, L)

  q <- phase[1:n]
  p <- phase[(n + 1):(2 * n)]

  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(m, current_q, given_other, flag)
  current_K = K(current_p)
  proposed_U = U(m, q, given_other, flag)
  proposed_K = K(p)

  alpha <-
    min(1, exp(current_U - proposed_U + current_K - proposed_K))
  if (is.nan(alpha))
    return(list(
      q = current_q,
      alpha = 0,
      accepted = FALSE
    ))

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < alpha) {
    return(list(q = q, alpha = alpha, accepted = TRUE))
    # accept
  } else {
    return(list(q = current_q, alpha = alpha, accepted = FALSE))
    # reject
  }
}
