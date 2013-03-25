\name{ci,Cblmr-method}
\alias{ci,Cblmr-method}
\alias{ci}

\title{Confidence Interval for the Changepoint}

\description{
  Confidence interval for theta, the changepoint's x-coordinate.
}

\usage{
\S4method{ci}{Cblmr}(  CL =0.95, method ="clr" )
}

\arguments{
  \item{CL}{confidence level, between 0 and 1}
  \item{method}{"clr" or "af" which stand for conditional likelihood-ratio, which is exact, or approximate-F}
}

\details{
This subroutine scans to determine the postulate values of theta that have significance levels greater than 1-CL.
}

\examples{
## Data for Patient B from Smith and Cook (1980):

x=c(1., 2., 3., 4., 5., 6., 7., 8., 9., 10.)
y=c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
sc <- blmr(x,y)
sc$ci()
sc$ci(0.90)
sc$ci(.99,"af")


}
