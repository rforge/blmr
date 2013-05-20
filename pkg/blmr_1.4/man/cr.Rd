\name{cr,Cblmr-method}
\alias{cr,Cblmr-method}
\alias{cr}

\title{Confidence Region for the Changepoint}

\description{
  Joint confidence region for ( theta, alpha ),  the changepoint's (x,y)-coordinates.
}

\usage{
\S4method{cr}{Cblmr}(  CL =0.95 ,  method ="clr",  incr =0.2 )
}

\arguments{
  \item{CL}{confidence level, between 0 and 1}
  \item{method}{"clr" or "af" which stand for conditional likelihood-ratio (exact)  or  approximate-F (fast)}
  \item{incr}{increment of theta values in the confidence region printout}
}

\details{
This subroutine scans to determine the postulate values of (theta,alpha) that have significance levels greater than 1-CL.  It scans first along the (theta, alpha-MLE) ridge to determine theta boundaries.
}

\examples{
##  Data for Patient B from Smith and Cook (1980):
y = c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x = 1:10
sc <- blmr(y,x)
sc$cr()
\donttest{sc$cr(.9,"clr",0.1)}
sc$cr(.99,"af")

}
