update_gamma <- function(U, grad_U, epsilon, L, current_beta, current_gamma){
  q = current_gamma

  # Sample random momenta
  p = rnorm(length(current_gamma),0,1)  # default: independent standard normal variables
  current_p = p

  # Compute deterministic trajectory depending on method
  phase <- leapfrog(p, q, epsilon, L, grad_U_gamma, current_beta)

  q <- phase[1:length(p)]
  p <- phase[(length(p) + 1):(2 * length(p))]

  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_beta, current_gamma)
  current_K = K(current_p)
  proposed_U = U(current_beta, q)
  proposed_K = K(p)

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position

  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
    current_gamma <- q
    return (q)
    # accept
  } else {
    return (current_gamma)
    # reject
  }
}
