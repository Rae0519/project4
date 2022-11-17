
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <-  function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}



newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
     # judge if the objective or derivatives are finite at the initial theta
     if ( Inf | -Inf %in% grad(th = theta,...) == TRUE){
       warning('the objective or derivatives are not finite at the initial theta')
     }

    iter = 0
    for (k in 1:maxit){
    
        nll0 <- func(th = theta,...)
        #if 
        #else
        grad0 <- grad(th = theta,...)
        hess0 <- hess(th = theta,...)
      
        if (class(try(chol(hess0))) == 'try-error'){
          warning('Hessian is not positive definite')
          hess0 <- (t(hess0) + hess0)/2
        }
 
        r = chol(hess0)
        I = diag(nrow(hess0))
        delta = (-1) * backsolve(r,forwardsolve(t(r),I))%*% grad0
      
        count = 0
        for (i in 1:(max.half+1)) {
            theta_n = theta + delta
            if (func(th = theta_n,...) < func(th = theta,...)){
              theta = theta_n
              break
            }
            else {
              delta = delta / 2
              count = count + 1
            }
        }
        if (count == 20){
        warning('step fails to reduce the objective despite')
        } 
        iter = iter + 1
    }
    # judge convergence
    if ((abs(grad(theta,...)) < tol*(func(theta,...) + fscale)) == FALSE){
      warning('maxit is reached without convergence')
    }
   list(f = func(th = theta,...) ,theta = theta, iter = iter, g = grad(th = theta,... ,Hi = hess(th = theta,...)))
}                          


