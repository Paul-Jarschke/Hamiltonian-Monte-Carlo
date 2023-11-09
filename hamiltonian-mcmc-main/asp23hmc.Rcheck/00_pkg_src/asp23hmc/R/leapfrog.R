leapfrog <- function(m, FUN, p, q, given_other, epsilon, L) {
  # Catch the special case that there is only 1 step
  if (L == 1) return(leapfrog_single_step(m, FUN, p, q, given_other, epsilon))

  # Make a half step for momentum at the beginning
  p = p - epsilon * FUN(m, q, given_other) / 2

  # Alternate full steps for position and momentum
  for (i in 1:(L-1)) {
    # Make a full step for the position
    q = q + epsilon * grad_K(p)

    # Make a full step for the momentum, except at end of trajectory
    p = p - epsilon * FUN(m, q, given_other)
  }

  # Make a full step for the position at the end.
  q = q + epsilon * grad_K(p)

  # Make a half step for momentum at the end.
  p = p - epsilon * FUN(m, q, given_other) / 2

  return(c(q, p))
}

leapfrog_single_step <- function(m, FUN, p, q, given_other, epsilon) {
  # Make a half step for momentum at the beginning
  p = p - epsilon * FUN(m, q, given_other) / 2

  # Make a full step for the position
  q = q + epsilon * grad_K(p)

  # Make a half step for momentum at the end.
  p = p - epsilon * FUN(m, q, given_other) / 2

  return(c(q, p))
}
