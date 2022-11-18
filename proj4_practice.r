# Group member1: Rui Zeng      ID:s2298181
# Group member2: Lifu Zheng    ID:s2314868
# Group member3: Sihan Liu     ID:s2337553



# The function is to implement Newton’s method for minimization of functions
# Input containing: initial values for the optimization parameters(theta), objective function to minimize(func), gradient function(grad),
# Hessian matrix function(hess), convergence tolerance(tol), a rough estimate of the magnitude of func near the optimum(fscale), 
# maximum number of Newton iterations(maxit), maximum number of times a step should be halved(max.half), finite difference intervals(eps)
# Output containing: the minimum value of the objective function(f), the minimum value of the parameters(theta),
# the number of iterations taken to reach the minimum(iter), gradient vector at the minimum(g), the inverse of the Hessian matrix at the minimum(Hi)

newt <- function(theta,func, grad, hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  # judge if the objective or derivatives are finite at the initial theta
  if ( FALSE %in% is.finite(func(theta,...)) | FALSE %in% is.finite(grad(theta,...)) | (if(is.null(hess) == FALSE){FALSE %in% is.finite(hess(theta,...))}else{FALSE})){
    stop('The objective or derivatives are not finite at the initial theta!')
  }
  iter = 0                 # define iteration
  m <- length(theta)       # the number of elements in theta
  # go through 'maxit' Newton iterations to find the minimum
  for (w in 1:maxit){
    
    nll0 <- func(theta,...)                    # the initial value of objective function - nll0
    grad0 <- grad(theta,...)                   # the initial gradient vector - grad0
    # if hessian is not provided
    if (is.null(hess) == TRUE){
      hess0 <- matrix(0,m,m)                   # create a m*m zero matrix - hess0
      for (i in 1:length(theta)){              # for each element in 'theta' vector
        theta_c <- theta                       # define a new theta
        theta_c[i] <- theta_c[i] + eps         # theta_c changes theta by adding eps
        grad_c <- grad(theta_c,...)            # grad_c is the gradient vector of changed theta
        hess0[i,] <- (grad_c - grad0) / eps    # 
      }
    }
    # if initial Hessian is provided
    else{
      hess0 <- hess(theta,...)                 # the initial Hessian matrix at theta - hess0
    }
    # if initial gradient vector is 0 and initial hessian matrix is positive definite
    # then it satisfies the minimum condition and terminate
    if (all(class(try(chol(hess0))) == 'try-error')){
      multiple <- 1e-6
      while ((all(class(try(chol(hess0), silent = TRUE)) == 'try-error')) == TRUE){
        hess0 <- hess0 + multiple*norm(hess0)*diag(dim(hess0)[1])
        multiple <-  10 * multiple
      }
    }
    else{
      if (all(grad0 == 0)){
        break
      }
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
      stop('The step fails to reduce the objective despite trying max.half step halvings!')
    } 
    
    iter = iter + 1
  }
  # judge convergence by seeing whether all absolute values of the gradient vector have less than 'tol' mutiple the absolute value of the objective function + 'fscale'
  if (all(abs(grad0) >= tol*(abs(nll0) + fscale))){
    warning('Maxit is reached without convergence!')
  }
  # if satisfy the convergence condition
  else{
    if(all(class(try(chol(hess0))) == 'try-error')){
    warning('Hessian is not positive definite!')
    }
  }
  # return a list at minimum consist of: value of objective function, value of parameters, the iteration, gradient vector and inverse Hessian matrix 
  list(f = nll0 ,theta = theta, iter = iter, g = grad0, Hi = chol2inv(chol(hess0)))
} 

















