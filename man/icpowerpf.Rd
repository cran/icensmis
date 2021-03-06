% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/icpowerpf.R
\name{icpowerpf}
\alias{icpowerpf}
\title{Study design in the presence of interval censored outcomes (assuming perfect
diagnostic tests)}
\usage{
icpowerpf(
  HR,
  survivals,
  N = NULL,
  power = NULL,
  rho = 0.5,
  alpha = 0.05,
  pmiss = 0
)
}
\arguments{
\item{HR}{hazard ratio under the alternative hypothesis.}

\item{survivals}{a vector of survival function at each test time for 
baseline(reference) group. Its length determines the number of tests.}

\item{N}{a vector of sample sizes to calculate corresponding powers. If one 
needs to calculate sample size, then set to NULL.}

\item{power}{a vector of powers to calculate corresponding sample sizes. If 
one needs to calculate power, then set to NULL.}

\item{rho}{proportion of subjects in baseline(reference) group.}

\item{alpha}{type I error.}

\item{pmiss}{a value or a vector (must have same length as survivals) of the 
probabilities of each test being randomly missing at each test time. If 
pmiss is a single value, then each test is assumed to have an identical 
probability of missingness.}
}
\value{
same form as returned value of \code{\link{icpower}}
}
\description{
This function implements power and sample size calculations for interval
censored time-to-event outcomes, when the diagnostic tests are assumed to be
perfect (i.e. sensitivity=1 and specificity=1). This is a special case of the
more general study design function \code{\link{icpower}}. However, for the
special case of perfect diagnostic tests, this function can be used with
significantly improved computational efficiency.
}
\note{
See \code{\link{icpower}} for more details in a general situation.
}
\examples{
powpf1 <- icpowerpf(HR =2 , survivals = seq(0.9, 0.1, by=-0.1), N = NULL,
   power = 0.9, pmiss = 0)
   
powpf2 <- icpowerpf(HR =2 , survivals = seq(0.9, 0.1, by=-0.1), N = NULL,
   power = 0.9, pmiss = 0.7)
   
## Different missing probabilities at each test time
powpf3 <- icpowerpf(HR =2 , survivals = seq(0.9, 0.1, -0.1), N = NULL, 
   power = 0.9, pmiss = seq(0.1, .9, 0.1))    
   
}
