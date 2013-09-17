\name{cr}
\alias{cr,Cpp_Clmbr-method}
\alias{cr}

\title{Confidence Region for the Changepoint}

\description{
  Joint confidence region for ( theta, alpha ),  the changepoint's (x,y)-coordinates.
}

\usage{
\S4method{cr}{Cpp_Clmbr}(  CL =0.95 ,  method ="clr",  incr  )
}

\arguments{
  \item{CL}{confidence level, between 0 and 1}
  \item{method}{"clr" or "af" which stand for conditional likelihood-ratio  or  approximate-F}
  \item{incr}{increment of theta values in the confidence region printout}
}

\details{
This subroutine scans to determine the postulate values of (theta,alpha) that have significance levels greater than 1-CL.  It scans first along the (theta, alpha-MLE) ridge to determine 'theta' boundaries.
}

\examples{
##  Data for Patient B from Smith and Cook (1980)
library(lm.br)
y = c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x = 1:10
sc <- lm.br(y~x)
\donttest{sc$cr()}
\donttest{sc$cr( .9, incr = 0.1 )}
sc$cr(.99,"af")

}
