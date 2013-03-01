\name{mle,Cblmr-method}
\alias{mle}

\title{Maximum Likelihood Estimates}

\description{
  Printout the maximum likelihood estimates of parameters.  Estimates are without bias correction except for the variance.
}

\usage{
\S4method{mle}{Cblmr}()
}

\examples{
## Data for Patient B from Smith and Cook(1980):

x=c(1., 2., 3., 4., 5., 6., 7., 8., 9., 10.)
y=c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
sc <- blmr(x,y)
sc$mle()


}
