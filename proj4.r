

# 
# newt <- function(theta, func, grad,hess= NULL, ..., tol = 1e-8, fscale = 1,maxit = 100, max.half = 20, eps = 1e-6){
#   
# delta <- -backsolve(hess(theta), forwardsolve(t(hess(theta)), grad(theta)))
# 
# 
# 
# while((func(theta) + t(delta)%*%grad(theta) + 1/2*t(delta)%*%hess(theta)%*%delta) >= func(theta)){
#   delta = 1/2*delta
# }
#   
# if (all(grad(theta) == 0) & all(eigen(hess(theta))$values >= 0)){
#   print(list(f, theta, iter, g, Hi))
# }
# 
# 
# } 


rb <- function(th, k=2) {    #   objective function
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <-  function(th,k=2) {   # gradient function
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2), k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {    # hessian function
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

#installed.packages(priority="recommended")

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  # judge if the objective or derivatives are finite at the initial theta
  if ( FALSE %in% is.finite(func(theta,...)) | FALSE %in% is.finite(grad(theta,...)) | FALSE %in% is.finite(hess(theta,...))){
    warning('the objective or derivatives are not finite at the initial theta')
  }
  iter = 0
  for (k in 1:maxit){
    
    nll0 <- func(th = theta,...)
    grad0 <- grad(th = theta,...)
    if (hess = NULL){
      hess0 <- matrix(0,2,2)
      for (i in 1:length(theta)){
        theta_c <- theta
        theta_c[i] <- theta_c[i] + eps
        grad_c <- grad(th = theta_c,...)
        hess0[i,] <- (grad_c - grad0) / eps
      }
    }
    else{
      hess0 <- hess(th = theta,...)
    }
   
    if (all(grad(theta) == 0) & all(eigen(hess(theta))$values >= 0)){
     break
    }
   
    if (all(class(try(chol(hess0))) == 'try-error')){
      hess0 <- (t(hess0) + hess0)/2
    }
    
    delta = -chol2inv(chol(hess0))%*% grad0
    count = 0
    
    for (i in 1:(max.half+1)) {
      theta_n = theta + delta
      #func_theta_delta <- nll0 + dela
      if (func(th = theta_n,...) < nll0){
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
  if (all((abs(grad(theta,...)) < tol*(func(theta,...) + fscale)) == FALSE)){
    warning('maxit is reached without convergence')
  }
  if (all(class(try(chol(hess0))) == 'try-error')){
    warning('Hessian is not positive definite')
  }
  list(f = func(th = theta,...) ,theta = theta, iter = iter, g = grad(th = theta,... ), Hi = chol2inv(chol(hess0)))
} 

# grad <- gb(theta) 
# eps =  1e-6
# hess = NULL
# for (i in 1:length(grad)){
#   hess[i] <- grad[i] + eps
# }
# hess1 = hb(theta)
# hess1
# 
# theta <- c(2, 3)
# newt(theta, rb, gb, hb)
# 
# hess0 <- hb(theta)

# r = chol(hess0)
# I = diag(nrow(hess0))
# delta = (-1) * backsolve(r,forwardsolve(t(r),I))%*% gb(theta)
# delta
# delta2 = -chol2inv(chol(hess0))%*% gb(theta)
# delta2
# start.time <- Sys.time()
# delta = slove(hess0)%*% gb(theta)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# x1 <- cbind(1, 1:10)
# norm(x1)
# norm(x1, "I")



# grad0 <- gb(theta) 
# hess0 <- matrix(0,2,2)
# for (i in 1:length(theta)){
#   theta_c <- theta
#   theta_c[i] <- theta_c[i] + eps
#   grad_c <- gb(theta_c)
#   hess0[i,] <- (grad_c - grad0) / eps
# }
# # hess <- hb(theta)
# hess0
# hess

# multiple <- 1e-6
# multiple <- multiple * 10
# hess <- hess + multiple * norm(hess) * diag(dim(hess))
# hess <- matrix(1:4, 2,2)
# hess1 <- (t(hess) + hess)/2
# hess1
# hess
# U <- eigen(hess)$vector
# gamma <- diag(abs(eigen(hess)$values))
# hess2 <- U%*%gamma%*%t(U)
# hess2
# hess1
