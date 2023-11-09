update_beta <- function(current_beta){
  new_beta <- HMC(q = current_beta)
}



update_beta <- function(U, grad_U_beta, epsilon, L, current_beta, current_gamma, method = "leapfrog", M){

  q = current_beta

  # Sample random momenta
  p = rnorm(length(current_beta),0,1)  # default: independent standard normal variables
  current_p = p

  # Compute deterministic trajectory depending on method
  phase <- leapfrog(p, q, epsilon, L, grad_U_beta, current_gamma)

  q <- phase[1:length(p)]
  p <- phase[(length(p) + 1):(2 * length(p))]

  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_beta, current_gamma)
  current_K = K(current_p)
  proposed_U = U(q, current_gamma)
  proposed_K = K(p)

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position

  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)) {
    current_beta <- q
    return (q)
    # accept
  } else {
    return (current_beta)
    # reject
  }
}















