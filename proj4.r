## Group member1: Rui Zeng      ID:s2298181
## Group member2: Lifu Zheng    ID:s2314868
## Group member3: Sihan Liu     ID:s2337553


## Github Repo: https://github.com/Rae0519/project4.git

## Contribution:
## We both add comments and integrate the code.
## Every group member has made contribution to the final code. 
## We all completed the practical work in our own way.
## The proportion of contribution is roughly equal in our group.

## Overview:
## We are writing a function about implementing Newton’s method for minimization of functions.
## Newton's method is about using Taylor's theorem to minimize successive quadratic approximations 
## to objective function (func) by minimizing parameters (theta). The small perturbation(delta) from Taylor's theorem should be:
## delta = - hess^(-1)%*%grad where 'grad' is the gradient vector and 'hess' is the Hessian Matrix of 'func'.
## (1) evaluate 'func', 'grad' and 'hess' by theta
## (2) we need to check the minimum of theta by checking whether grad = 0 and 'hess' is positive definite.
## If it meets, terminate. 
## (3) If not and 'hess' is not positive definite, we perturb 'hess' with adding a multiple of the identity matrix to it, 
## large enough to force positive definiteness.
## (4) we need to halve 'delta' until func(theta + delta) < func(theta)
## (5) Let theta = theta + delta, return to step 1.


newt <- function(theta,func, grad, hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
## 'newt' is the function of using Newton's method for minimization of functions.
## Input containing: 'theta'; 'func'; 'grad'; 'hess'- Hessian Matrix, if not supplied obtaining by finite differencing of 'grad'; 
## 'tol' - convergence tolerance; 'fscale' - a rough estimate of the magnitude of func near the optimum for convergence testing;
## 'maxit'- maximum number of Newton iterations; 'max.half' - maximum number of times a step should be halved
## 'eps' - finite difference intervals for hessian matrix if not supplied;
## Output containing: the minimum value of the objective function(f); the minimum value of the parameters(theta);
## the number of iterations taken to reach the minimum(iter); gradient vector at the minimum(g);
## the inverse of the Hessian matrix at the minimum(Hi)
  
## judge if the objective or derivatives are finite at the initial theta            ##   check if 'hess' is null, if null, FALSE; if not,check whether it is finite.
  if ( FALSE %in% is.finite(func(theta,...)) | FALSE %in% is.finite(grad(theta,...)) | (if(is.null(hess) == FALSE){FALSE %in% is.finite(hess(theta,...))}else{FALSE})){
    stop('The objective or derivatives are not finite at the initial theta!')    
  }
  iter = 0                 ## define iteration time
  m <- length(theta)       ## the number of elements in theta
  
## loop through 'maxit' Newton iterations to find the minimum: we need to restric the time for return to step 1
  for (w in 1:maxit){
    nll0 <- func(theta,...)    ## the initial value of objective function - nll0
    grad0 <- grad(theta,...)   ## the initial gradient vector - grad0
    
## if hessian is not provided => hess = NULL
    if (is.null(hess) == TRUE){ 
      hess0 <- matrix(0,m,m)       ## create a mxm zero matrix - hess0
      for (i in 1:length(theta)){  ## construct 'theta_c' which has the same length as 'theta'
        theta_c <- theta      
        theta_c[i] <- theta_c[i] + eps   ## 'theta_c' adds eps to 'theta'
        grad_c <- grad(theta_c,...)      ## 'grad_c' is the gradient vector of 'theta_c'
        hess0[i,] <- (grad_c - grad0) / eps   ## differentiating 'grad_c' to derive the hessian matrix
      }
    }
## if initial Hessian is provided
    else{
      hess0 <- hess(theta,...)   ## the initial Hessian matrix at theta - hess0
    }
## check if 'hess0' is positive definite by Cholesky decomposition
    if (all(class(try(chol(hess0), silent = TRUE)) == 'try-error')){     ## if not positive definite
      multiple <- 1e-6   ## set the initial multiple                                ## try() all code for fail and silent for not return error
      while ((all(class(try(chol(hess0), silent = TRUE)) == 'try-error')) == TRUE){  ## run until positive definite
        hess0 <- hess0 + multiple*norm(hess0)*diag(dim(hess0)[1])  ##  hess = hess + sigma*itentity matrix  
        multiple <-  10 * multiple          ## if still not positive definite, enlarge multiple
      }
    }
    else{  ## when 'hess' is positive definite
      if (all(grad0 == 0)){  ## check ifinitial gradient vector is 0
        break                ## then it satisfies the minimum condition and terminate
      }
    }
    ## step 4
    delta = -chol2inv(chol(hess0))%*% grad0  ## 'delta' formula which will construct a mx1 matrix
    count = 0   ##  record the time for halving delta
    for (i in 1:(max.half+1)) {   ##  loop through max.half + 1 for check the delta for max.half times
      theta_n = theta + delta
      if (func(theta_n,...) < nll0){  ## when func(theta + delta) < func(theta)
        theta = theta_n   ## step 5: change for the new theta
        break  ## stop halving delta
      }
      else {  ## when func(theta + delta) >= func(theta)
        delta = delta / 2  ## half
        count = count + 1  ##  record the time of halving
      }
    }
## If the step fails to reduce the objective despite trying max.half step halvings
    if (count == max.half + 1){  ## since we loop through max.half + 1 times above
      stop('The step fails to reduce the objective despite trying max.half step halvings!') 
    } 
    
    iter = iter + 1  ## record the iteration time
  }
## judge convergence by seeing whether all absolute values of the gradient vector have less than 'tol' mutiple the absolute value of the objective function + 'fscale'
  if (all(abs(grad0) >= tol*(abs(nll0) + fscale))){
    warning('Maxit is reached without convergence!')
  }
  else{ ## when satisfy the convergence condition
## check if the Hessian is not positive definite at convergence.
    if(all(class(try(chol(hess0), silent = TRUE)) == 'try-error')){
    warning('Hessian is not positive definite!')
    }
  }
## return a list at minimum consist of: value of objective function, value of parameters, the iteration, gradient vector and inverse Hessian matrix   
  list(f = nll0 ,theta = theta, iter = iter, g = grad0, Hi = chol2inv(chol(hess0)))
} 







