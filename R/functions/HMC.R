#' @title The function to fit Hamiltonian Monte Carlo sampling
#' @description Takes the the likelihood and the gradient of parameter of interest
#' and updates the prior to generate a new proposal 
#' @references 
#'
#' @param U which returns the potential energy given a value for q
#' @param grad_U which returns the vector of partial derivatives of U given q
#' @param epsilon is the stepsize
#' @param L is the number of leapfrog steps in the trajectory
#' @param current_q is the current position that the trajectory starts from
#' @param arc is a control parameter incremented by one if the new proposal is accepted

#' @return HMC returns a list containing either the position at the end of the trajectory 
#' or the initial position (up) and the control parameter to record acceptance of proposal (arc)

HMC = function (U, grad_U, epsilon, L = 30, current_q, arc)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  
  for (j in 1:L)
  {
    # Make a full step for the position
    
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    
    if (j!=L) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  
  R <- current_U-proposed_U+current_K-proposed_K
  if(is.nan(R)){R = -Inf}
  if (log(runif(1)) < R)
  {
    up = q  # accept
    arc <- arc + 1
  }
  else
  {
    up = current_q
  }
  return(list(up = up, arc = arc))  # reject
}
