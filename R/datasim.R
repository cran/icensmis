datasim <-
function(N, blambda, testtimes, sensitivity, specificity,  
					  betas = NULL, twogroup = NULL, pmiss = 0, design = "MCAR", negpred = 1){
    
    ## Basic input check
    if (N < 0) stop("invalid N")
    if (sensitivity > 1 | sensitivity < 0) stop("invalid sensitivity")
    if (specificity > 1 | specificity < 0) stop("invalid specificity")
    if (!is.null(twogroup)) {if (twogroup > 1 | twogroup < 0) stop("invalid twogroup")}
    if ((length(pmiss)!=1 & length(pmiss)!=length(testtimes)) | (sum(pmiss>1 | pmiss<0)>0)) {stop("invalid pmiss")}
    if (!(design %in% c("MCAR", "NTFP"))) stop("invalid design")
    if (negpred > 1 | negpred < 0) stop("invalid negpred")    
    
	nbeta <- length(betas)
	ntest <- length(testtimes)
 
	## Create subject profiles
	ID            <- 1:N
	if (is.null(betas)) {
		  prog <- rep(0, N)
	  } else if (is.null(twogroup)) {
		  covm           <-  matrix(rnorm(nbeta*N), ncol=nbeta)
		  prog           <-  covm%*%betas 
	      colnames(covm) <- paste("cov", 1:nbeta, sep="")
	  } else if (twogroup > 0 & twogroup < 1) {
		  N1 <- round(N*twogroup)
		  N2 <- N - N1
		  covm           <- matrix(c(rep(0, N1), rep(1, N2)), ncol = 1)
		  prog           <- covm*betas[1]
		  colnames(covm) <- "group"
	  } else {stop("invalid input twogroup")}

	lambdas       <- blambda*exp(prog)
	EventTime     <- rexp(N, lambdas)
	basepos       <- rbinom(N, 1, 1-negpred)==1  ## Whether baseline is positive or not
	EventTime[basepos] <- 0                      ## Set baseline positive's event time to be 0
	if (is.null(betas)) {
 		  profile       <- data.frame(ID, prog, lambdas, EventTime)
 	  } else {
 		  profile       <- data.frame(ID, prog, lambdas, EventTime, covm)
 	}

	## Exapand profiles and test times
	IDs           <- data.frame(ID=rep(ID, each=ntest))
	exprofile     <- merge(IDs,profile, by="ID", all=T)
	testtime      <- rep(testtimes, times=N)

	## Obtain test results
	occur   <- testtime > exprofile$EventTime
	probs   <- ifelse(occur, sensitivity, 1-specificity)
	result  <- rbinom(length(occur), 1, probs)

  	## combine data and apply random missing
	sim.data    <- data.frame(exprofile, testtime, result)[, -(2:4)]
	notmiss     <-  rbinom(dim(sim.data)[1], 1, 1-pmiss)
	sim.data    <- sim.data[notmiss==1, ]

	## Apply missing after first positive (if design is NTFP)
	if (design=="NTFP") {
		stopositive <- function(x){
     	   idx <- which(x$result==1)[1]
     		if (!is.na(idx)) {x <- x[1:idx, ]} 
     		return(x)
    	}
		sim.data  <- do.call("rbind", by(sim.data, sim.data$ID, stopositive))
	}
	
	## return simulated data
	rownames(sim.data) <- NULL
	return(sim.data)
}
