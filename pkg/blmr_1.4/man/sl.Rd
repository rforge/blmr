\name{sl,Cblmr-method}
\alias{sl,Cblmr-method}
\alias{sl}

\title{Significance Level for Changepoint}

\description{
  Significance level of a postulate value for the changepoint's x-coordinate or (x,y)-coordinates.
}

\usage{
\S4method{sl}{Cblmr}( theta0,  method ="clr", accuracy =0.001, output =TRUE )
sl( theta0, alpha0,  method ="clr", accuracy =0.001, output =TRUE )
}

\arguments{
  \item{theta0}{postulate value for theta, the changepoint's x-coordinate }
  \item{alpha0}{postulate value for alpha, the changepoint's y-coordinate }
  \item{method}{"clr", "mc" or "af" which stand for conditional likelihood-ratio, conditional likelihood-ratio by Monte Carlo, or approximate-F, details below}
  \item{accuracy}{maximum absolute error in numerical integration by the "clr" method, or in Monte Carlo evaluation by the "mc" method, has no effect for the "af" method}
  \item{output}{TRUE or FALSE}
}

\details{

The default evaluation method "clr" evaluates the conditional likelihood-ratio significance level using a geometric-expectation formula derived by Siegmund-Zhang, which slightly over-estimates.

Method "mc" evaluates the conditonal likelihood-ratio significance level by direct Monte Carlo simulation, which avoids the over-estimate of the "clr" method.

Method "af" evaluates significance level by estimating the distribution of the likelihood-ratio statistic from the related F-distribution (or chi-squared distribution if variance is known) that would be exact if the model were linear.  This method is not exact.  It is common in nonlinear regression.

}

\note{
The 'accuracy' argument does not include the slight over-estimate that is inherent in the "clr" method, nor the approximation inherent in the Approximate-F method.
}

\value{
Returns the numerical significance level if 'output' is set explicitly to TRUE or FALSE, otherwise only prints the result.
}

\examples{
##  Data for Patient B from Smith and Cook (1980):
y = c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x = 1:10
sc <- blmr(y,x)

sc$sl(6.1)
sc$sl(6.1,"af")
\donttest{sc$sl(6.1,"mc",0.00001)}
tmp <- sc$sl(6.1,"clr",0.001,FALSE)
tmp

sc$sl(6.1,88.2)
sc$sl(6.1,88.2,"clr",0.00001)
sc$sl(6.1,88.2,"af")
\donttest{sc$sl(6.1,88.2,"mc",0.00001)}
tmp2 <- sc$sl(6.1,88.2,"af",.001,TRUE)
tmp2

}
