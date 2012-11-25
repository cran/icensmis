icpower <- function(HR, sensitivity, specificity, survivals, N = NULL, power = NULL,
                   rho = 0.5, alpha = 0.05, pmiss = 0, design = "MCAR", negpred = 1) {

   	## Basic input check
   	if (!(HR>0)) stop("Check input for HR")
   	if (!(sensitivity<=1 & sensitivity>0)) stop("Check input for sensitivity")
   	if (!(specificity<=1 & specificity>0)) stop("Check input for specificity")
   	if (sum(diff(survivals)>0)>0 | sum(survivals<0)>0 | sum(survivals>1)>0) stop("Check input for survivals")
   	if (rho<0 | rho>1) stop("Check input for rho")
   	if (alpha<0 | alpha>1) stop("Check input for alpha")
   	if (sum(pmiss<0)>0 | sum(pmiss>1)>0) stop("Check input for pmiss")
   	if ((is.null(N)+is.null(power))!=1) stop("Need input for one and only one of power and N")
   	if (!is.null(power)) {if (any(power > 1) | any(power < 0)) stop("Check input for power")}
   	if (length(pmiss)!=1 & length(pmiss)!=length(survivals))
       stop("length of pmiss must be 1 or same length as survivals")
  	if (!(design %in% c("MCAR", "NTFP"))) stop("invalid design")
  	if (negpred<0 | negpred>1) stop("Check input for negpred")
  	
  	## Use dedicated function for perfect test, so when perfect test use alternative instead
   	if (sensitivity==1 & specificity==1) {
   		cat("For perfect test, use icpowerpf for improved computational efficiency... \n \n")
   		design <- "NTFP"
   		}
   	
	miss 	<- any(pmiss!=0)  
   	beta  	<- log(HR)
   	J 		<- length(survivals)
   	Surv  	<- c(1, survivals)  	                   	
	if (length(pmiss)==1) pmiss <- rep(pmiss, J)

	if (design=="MCAR" & !miss) {	
		## MCAR design with no missing tests
		data 	<- as.matrix(do.call("expand.grid", rep(list(0:1), J))) 
		o		<- 2*outer(1:J, 1:(J+1), ">=") + 1
		smat 	<- c(specificity, 1-specificity, 1-sensitivity,  sensitivity)
		Cm 		<- sapply(1:J, function(x) smat[outer(data[, x], o[x, ], "+")])
		Cm 		<- apply(Cm, 1, prod)
		Cm		<- matrix(Cm, ncol=J+1)
		prob 	<- 1
	
	} else if (design=="MCAR" & miss) {
		## MCAR design with missing tests
		data 	<- as.matrix(do.call("expand.grid", rep(list(0:2), J)))[-3^J, ]
		o 		<- 3*outer(1:J, 1:(J+1), ">=") + 1 
		smat 	<- c(specificity, 1-specificity, 1, 1-sensitivity, sensitivity, 1)
		Cm	 	<- sapply(1:J, function(x) smat[outer(data[, x], o[x, ], "+")])
		Cm	 	<- apply(Cm, 1, prod)
		Cm	 	<- matrix(Cm, ncol=J+1)
		prob 	<- apply(pmiss^(t(data==2))*(1-pmiss)^(t(data!=2)), 2, prod)

	} else if (design=="NTFP" & !miss) {
		## NTFP without missing tests
		data <- sapply(0:(J-1), function(x) c(rep(0, x), 1, rep(2, J-x-1)))
		data <- t(cbind(data, rep(0, J)))
		o 		<- 3*outer(1:J, 1:(J+1), ">=") + 1 
		smat <- c(specificity, 1-specificity, 1, 1-sensitivity,  sensitivity, 1)
		Cm	 	<- sapply(1:J, function(x) smat[outer(data[, x], o[x, ], "+")])
		Cm	 	<- apply(Cm, 1, prod)
		Cm	 	<- matrix(Cm, ncol=J+1)
		prob	<- 1
			
	} else if (design=="NTFP" & miss) {
		## NTFP with missing tests
		smat <- c(specificity, 1-specificity, 1, 1-sensitivity,  sensitivity, 1)
		crmat <- function(n) {
			if (n==0) {
				data <- matrix(1)		  
			} else if (n==J) {
				data <- as.matrix(do.call("expand.grid", rep(list(c(0, 2)), n)))[-2^J, ]
			} else {
				data <- cbind(as.matrix(do.call("expand.grid", rep(list(c(0, 2)), n))), 1)
			}
			n1 	<- min(n+1, J)
			o 	<- 3*outer(1:n1, 1:(J+1), ">=") + 1
			mats <- sapply(1:n1, function(x) smat[outer(data[, x], o[x, ], "+")])
			mats <- apply(mats, 1, prod)
			mats <- matrix(mats, ncol=J+1)
			p 	<- pmiss[1:n1]
			probs <-  apply(p^(t(data==2))*(1-p)^(t(data!=2)), 2, prod) 
			mats <- cbind(mats, probs)
			return(mats)
		}
		Cm <- do.call("rbind", sapply(0:J, crmat))
		prob	<- Cm[, J+2]
		Cm		<- Cm[, -(J+2)]
	}

	## Create transformation matrix and obtain D matrix
   	mat     <- c(1, rep(0, J), rep(c(-1, 1, rep(0, J)), J-1), -1, 1)
   	S2theta <- matrix(mat, ncol=J+1)
   	C2C     <- negpred*diag(J+1) + (1-negpred)*matrix(c(1, rep(0,J)),
   		       nrow=J+1, ncol=J+1)   ## Adjust for baseline misclassification  
   	Dm      <- Cm%*%C2C%*%S2theta

	I1 <- I2 <- matrix(NA, nrow=J+1, ncol=J+1)
	
	## Group 1 information matrix
    DS         	<- c(prob/(Dm%*%Surv))
    I1[1, ] 	<- I1[, 1] 	<- 0
    I1[-1, -1] 	<- (t(DS*Dm)%*%Dm)[-1, -1]
    
    ## Group 2 information matrix	
	a 		<- HR
	Surv2 	<- Surv^a
	Surv21	<- Surv^(a-1)
	lSurv 	<- log(Surv)
	DS2		<- c(1/(Dm%*%Surv2))
	
	I2[1, 1] <- sum(prob*((Dm%*%(Surv2*lSurv))^2*a^2*DS2 - Dm%*%(Surv2*lSurv*(1+a*lSurv))*a))
	I2[-1, -1] <- (t(Dm*prob*DS2)%*%Dm*outer(Surv21, Surv21, "*")*a^2)[-1, -1] 	
	I2[cbind(2:(J+1), 2:(J+1))] <- colSums(prob*t(Surv21^2*t(Dm^2*DS2)*a^2 - t(Dm)*Surv^(a-2)*a*(a-1)))[-1]	
	I2[1, 2:(J+1)] <- I2[2:(J+1), 1] <- colSums(prob*t(t(c(Dm%*%(Surv2*lSurv))*DS2*Dm)*Surv21*a^2-t(Dm)*Surv21*(1+a*lSurv)*a))[-1]
	
	If       <- rho*I1+(1-rho)*I2
   	inv.If   <- solve(If)
   	beta.var <- inv.If[1]

   	## Calculate Sample size or Power
   	if (!is.null(power)) {
    	N  <- ceiling((qnorm(1-alpha/2)+qnorm(power))^2*beta.var/beta^2)
        N1 <- round(rho*N)
        N2 <- N-N1
        return(list(result=data.frame(N,N1,N2,power), I1=I1, I2=I2))
   	} else {
        beta.varN <- beta.var/N
        za    <- qnorm(1-alpha/2)
        zb    <- beta/sqrt(beta.varN)
        power <- pnorm(-za,zb,1)+1-pnorm(za,zb,1)        
        N1 <- round(rho*N)
        N2 <- N-N1        
        return(list(result=data.frame(N,N1,N2,power), I1=I1, I2=I2))
  	}
}	
