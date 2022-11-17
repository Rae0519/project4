
# The function is to implement Newton’s method for minimization of functions
# Input containing: initial values for the optimization parameters(theta), objective function to minimize(func), gradient function(grad),
# Hessian matrix function(hess), convergence tolerance(tol), a rough estimate of the magnitude of func near the optimum(fscale), 
# maximum number of Newton iterations(maxit), maximum number of times a step should be halved(max.half), finite difference intervals(eps)
# Output containing: the minimum value of the objective function(f), the minimum value of the parameters(theta),
# the number of iterations taken to reach the minimum(iter), gradient vector at the minimum(g), the inverse of the Hessian matrix at the minimum(Hi)

newt <- function(theta,func, grad, hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  # judge if the objective or derivatives are finite at the initial theta
  if ( FALSE %in% is.finite(func(theta,...)) | FALSE %in% is.finite(grad(theta,...)) | (if(is.null(hess) == FALSE){FALSE %in% is.finite(hess(theta,...))}else{FALSE})){
    warning('The objective or derivatives are not finite at the initial theta!')
  }
  iter = 0
  for (w in 1:maxit){
    
    nll0 <- func(theta,...)
    grad0 <- grad(theta,...)
    
    if (is.null(hess) == TRUE){
      hess0 <- matrix(0,2,2)
      for (i in 1:length(theta)){
        theta_c <- theta
        theta_c[i] <- theta_c[i] + eps
        grad_c <- grad(theta_c,...)
        hess0[i,] <- (grad_c - grad0) / eps
      }
    }
    else{
      hess0 <- hess(theta,...)
    }
   
    if (all(grad0 == 0) & all(eigen(hess0)$values >= 0)){
     break
    }
   
    if (all(class(try(chol(hess0))) == 'try-error')){
      hess0 <- (t(hess0) + hess0)/2
    }
    
    delta = -chol2inv(chol(hess0))%*% grad0
    count = 0
    for (i in 1:(max.half+1)) {
      theta_n = theta + delta
      if (func(theta_n,...) < nll0){
        theta = theta_n
        break
      }
      else {
        delta = delta / 2
        count = count + 1
      }
    }
    if (count == max.half + 1){
      warning('The step fails to reduce the objective despite trying max.half step halvings!')
    } 
    
    iter = iter + 1
  }
  # judge convergence
  if (all(abs(grad0) >= tol*(abs(nll0) + fscale))){
    warning('Maxit is reached without convergence!')
  }
  else{
    if(all(class(try(chol(hess0))) == 'try-error')){
    warning('Hessian is not positive definite!')
    }
  }
  
  list(f = nll0 ,theta = theta, iter = iter, g = grad0, Hi = chol2inv(chol(hess0)))
} 




# 
# rb <- function(th, k=2) {    #   objective function
#   k*(th[2]-th[1]^2)^2 + (1-th[1])^2
# }
# 
# gb <-  function(th,k=2) {   # gradient function
#   c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2), k*2*(th[2]-th[1]^2))
# }
# 
# hb <- function(th,k=2) {    # hessian function
#   h <- matrix(0,2,2)
#   h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
#   h[2,2] <- 2*k
#   h[1,2] <- h[2,1] <- -4*k*th[1]
#   h
# }
# 
# theta = c( 1.3, 0.4)
# newt(theta, rb, gb)




