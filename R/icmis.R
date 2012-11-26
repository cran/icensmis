icmis <-
function(subject, testtime, result, data, sensitivity, specificity, formula = NULL,
		 negpred = 1, betai = NULL){
  
  	## Basic input check
  	if (sensitivity > 1 | sensitivity < 0) stop("invalid sensitivity")
    if (specificity > 1 | specificity < 0) stop("invalid specificity")
    if (negpred > 1 | negpred < 0) stop("invalid negpred")    

  	## Obtain result, time, id
  	result <- eval(substitute(result), data, parent.frame())
  	time   <- eval(substitute(testtime), data, parent.frame())
  	id     <- eval(substitute(subject), data, parent.frame())

	## Obtain C matrix
  	boundry <- c(sort(unique(time)), Inf)
  	o       <- t(sapply(time, function(x) x>=boundry))
  	r       <- (result==1)
  	likmat  <- (r & o)*sensitivity + (!r & o)*(1-sensitivity) +
               (r & !o)*(1-specificity) + (!r & !o)*specificity
  	Cm      <- apply(likmat, 2, function(x) tapply(x, id, prod)) 

  	## Create transformation matrix, then obtain D matrix
    ntest   <- length(boundry) - 1
  	mat     <- c(1,rep(0, ntest), rep(c(-1, 1, rep(0,ntest)), ntest-1), -1, 1)
  	S2theta <- matrix(mat, ncol=ntest+1)

  	C2C     <- negpred*diag(ntest+1) + (1-negpred)*matrix(c(1,rep(0,ntest)),
  	           nrow=ntest+1, ncol=ntest+1)   ## Adjust for baseline misclassification
  	S2theta1 <- C2C%*%S2theta           ## 10-14 change S2theta -> S2theta1
  	Dm      <- Cm%*%S2theta1            ## 10-14 change S2theta -> S2theta1

  	if (is.null(formula)) {    ##### If there is no covariate
     	
     	## Log-likelihood function
     	loglik <- function(parm) {
         	return(-sum(log(Dm%*%c(1,parm))))
       	}
     	
     	## Gradient function
     	gradlik <- function(parm) {
         	Sur   <- c(1, parm)
         	deriv <- colSums(Dm/as.vector(Dm%*%Sur))
         	return(-deriv[-1])
       	}
     	
     	## Initial values for parameters and contraints
     	parmi <- (ntest:1)/(ntest+1)
     	ui    <- S2theta[, -1]
     	ci    <- c(-1,rep(0, ntest))
		
		## Optimization
     	q <- constrOptim(parmi, loglik, gradlik, ui, ci, hessian=T)

     	## Create output
     	loglik   <- -q$value
     	time     <- boundry[1:ntest]
     	surv     <- q$par
     	survival <- data.frame(time,surv)
     	cov      <- solve(q$hessian)
     	name     <- paste("time=",time,sep="")
     	rownames(cov) <- colnames(cov) <- name
     	return(list(loglik=loglik,survival=survival,covariance=cov))
     	
     } else {                      ##### If there are covariates
     
     	## Obtain design matrix
     	designm <- model.matrix(formula, data=data)  
     	beta.nm <- colnames(designm)[-1]
     	design  <- apply(as.matrix(designm[ ,-1]), 2, function(x) tapply(x, id, function(y) y[1]))
     	nbeta   <- dim(design)[2]
		if (!is.null(betai)) {
			if (length(betai)!=nbeta) stop ("betai should be NULL or same length as the number of covariates")}
			
     	## Log-likelihood function
     	loglik <- function(parm) {
        	beta   <- parm[1:nbeta]
        	Sur    <- c(1, parm[-(1:nbeta)])
        	prog.score <- design%*%beta
        	Surv   <- sapply(Sur, function(x) x^exp(prog.score))
        	likmat <- rowSums(Dm*Surv)
        	return(-sum(log(likmat)))
      	}

     	## Gradient function
     	gradlik <- function(parm) {
        	beta   <- parm[1:nbeta]
        	Sur    <- c(1, parm[-(1:nbeta)])
        	prog.score <- design%*%beta
        	eprog  <- as.vector(exp(prog.score))
        	Surv   <- sapply(Sur, function(x) x^eprog)
        	Surv.d <- sapply(Sur, function(x) x^(eprog-1))
        	Deriv.Surv <- colSums(Dm*Surv.d*(eprog/rowSums(Dm*Surv))) 
        	Deriv.beta <- colSums(rowSums(Dm*Surv*log(Surv))/rowSums(Dm*Surv)*design)
        	return(-c(Deriv.beta, Deriv.Surv[-1]))
       	}
   
     	## Initial values
     	if (is.null(betai)) {
     		betai <- rep(0, nbeta)
     	}
     	parmi  <- c(betai, (ntest:1)/(ntest+1))
     	Sui    <- S2theta[, -1]
     	betaui <- matrix(0, nrow=ntest+1, ncol=nbeta)
     	ui     <- cbind(betaui, Sui)
     	ci     <- c(-1,rep(0, ntest))  

     	## Optimization
     	q <- constrOptim(parmi, loglik, gradlik, ui, ci, hessian=T)

     	## Create output
     	loglik   <- -q$value
     	time     <- boundry[1:ntest]
     	surv     <- q$par[-(1:nbeta)]
     	survival <- data.frame(time, surv)
     	cov      <- solve(q$hessian)
     	name     <- c(beta.nm, paste("time=", time, sep=""))
     	rownames(cov) <- colnames(cov) <- name
     	beta.fit <- q$par[1:nbeta]
     	beta.sd  <- sqrt(diag(cov)[1:nbeta])
     	beta.z   <- beta.fit/beta.sd
     	p.value  <- 2*(1-pnorm(abs(beta.z)))
     	coef     <- data.frame(coefficient=beta.fit, SE=beta.sd, z=beta.z, p.value=p.value)
     	return(list(loglik=loglik, coefficient=coef, survival=survival, covariance=cov)) 
   	}
}
