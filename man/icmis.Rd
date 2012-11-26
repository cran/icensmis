\name{icmis}
\alias{icmis}
\title{
Maximum likelihood estimation for settings of error-prone diagnostic tests and self-reported outcomes
}
\description{
This function estimates the baseline survival function evaluated at each test time in the presence of error-prone diagnostic tests and self-reported outcomes. If there are covariates included in the dataset, it also estimates their coefficients assuming proportional hazards. The function can also be used to incorporate misclassification of disease status at baseline (due to an error-prone diagnostic procedure).
}
\usage{
icmis(subject, testtime, result, data, sensitivity, specificity, formula = NULL,
		 negpred = 1, betai = NULL)
}
\arguments{
  \item{subject}{
variable in data for subject id.
}
  \item{testtime}{
variable in data for test time.
}
  \item{result}{
variable in data for test result.
}
  \item{data}{
the data to analyze.
}
  \item{sensitivity}{
the sensitivity of test.
}
  \item{specificity}{
the specificity of test.
}
  \item{formula}{
a formula to specify what covariates to be included in the model. If there is no covariate or one sample setting, set it to NULL. Otherwise, input like ~x1 + x2 + factor(x3).
}
  \item{negpred}{
baseline negative predictive value, i.e. the probability of being truely disease free for those who were tested (reported) as disease free at baseline. If baseline screening test is perfect, then negpred = 1. 
}
  \item{betai}{
a vector of initial values for the regression coefficients corresponding to the vector of covariates. If betai=NULL, then 0s are used for the initial values. Otherwise, the length of betai must equal the number of covariates. The initial values of survival function at each test time are automatically generated assuming equal interval probabilities.	
}
}
\details{
The input data should be in longitudinal form with one row per test time. Use \code{\link{datasim}} to simulate a dataset to see the sample data structure. If an error is generated due to the optimization procedure, then we recommend trying different initial values. This likelihood-based approach is a function of the survival function evaluated at each unique test time in the dataset and the vector of regression coefficients as model parameters. Therefore, it works best for situations where there is a limited number of unique test times in the dataset. If there are a large number of unique test times, one solution is to group several test times together. 
}
\value{
A list of fitting results is returned. 
\item{loglik}{log-likelihood function at MLE}
\item{coefficient}{a data frame with estimated coefficients for the covariates if covariates are inluded in the model}
\item{survival}{a data frame with estimated baseline survival function at each test time}
\item{covariance}{estimated covariance matrix of all parameters}
}
\author{
 Xiangdong Gu and Raji Balasubramanian
}

\examples{
## One sample setting
simdata1 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7, specificity = 0.98,  
		betas = NULL, twogroup = NULL, pmiss = 0.3, design = "MCAR", negpred = 1)
fit1 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata1, sensitivity = 0.7, 
		specificity= 0.98, formula = NULL, negpred = 1)				  
					  
## Two group setting, and the two groups have same sample sizes
simdata2 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7, specificity = 0.98,  
		betas = 0.7, twogroup = 0.5, pmiss = 0.3, design = "MCAR", negpred = 1)					  
fit2 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata2, sensitivity = 0.7,
		specificity= 0.98, formula = ~group, negpred = 1)						  
					  
## Three covariates with coefficients 0.5, 0.8, and 1.0
simdata3 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7, specificity = 0.98,  
		betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3, design = "MCAR", negpred = 1)					  
fit3 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata3, sensitivity = 0.7, 
		specificity= 0.98, formula = ~cov1+cov2+cov3, negpred = 1)
			 
## Fit data with NTFP missing mechanism (the fitting is same as MCAR data)
simdata4 <- datasim(N = 1000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7, specificity = 0.98,  
		betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3, design = "NTFP", negpred = 1)						  
fit4 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata4, sensitivity = 0.7,
		specificity= 0.98, formula = ~cov1+cov2+cov3, negpred = 1)					  
					  
## Fit data with baseline misclassification
simdata5 <- datasim(N = 2000, blambda = 0.05, testtimes = 1:8, sensitivity = 0.7, specificity = 0.98,  
		betas = c(0.5, 0.8, 1.0), twogroup = NULL, pmiss = 0.3, design = "MCAR", negpred = 0.97)						  
fit5 <- icmis(subject = ID, testtime = testtime, result = result, data = simdata5, sensitivity = 0.7, 
		specificity= 0.98, formula = ~cov1+cov2+cov3, negpred = 0.97)	
}
