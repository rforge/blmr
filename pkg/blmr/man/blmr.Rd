\name{blmr}
\alias{blmr-package}
\alias{blmr}

\docType{package}

\title{ Broken Line Model Regression }

\description{
Confidence intervals and confidence regions with exact coverage probabilities for the changepoint in broken line regression.

Version 1.2,  2013-03-04,  license:  GPL (>=2)
}


\usage{
blmr(  x,  y,  S,   model =1,  var_known =FALSE )
}

\arguments{
  \item{x}{vector of x values, non-decreasing}
  \item{y}{vector of y values}
  \item{S}{(optional) weights or covariance matrix, positive-definite, default is the identity matrix}
  \item{model}{ 1 or 2, defined below }
  \item{var_known}{FALSE if variance is unknown, TRUE if known}
}


\details{

A broken line model consists of two straight lines joined at a changepoint.  The versions supported here are

1:   y = alpha + B' * min( x-theta , 0 ) + B * max( x-theta, 0 ) + e

2:   y = alpha +  B * max( x-theta , 0 ) + e

e ~ N( 0 , var*S ).

Inferences about the changepoint  theta  or  (theta,alpha)  are based on the distribution of its likelihood-ratio statistic, conditional on sufficient statistics for the other parameters.

Function 'blmr' creates a data+model Cblmr object.  Then method 'ci' gets confidence intervals, 'cr' confidence regions, and 'sl' significance levels for the changepoint's x-coordinate or (x,y)-coordinates.  Other methods are 'mle' to get maximum likelihood estimates and 'sety' to set new y-values.

Test for a changepoint by the significance level of a postulate value outside the x-values.  Thus,  sl(x[1]-1)  gives the exact significance level of the null hypothesis H0:"single line" versus the alternate hypothesis H1:"broken line."

}

\note{
S is a weights matrix if variance is unknown, a covariance matrix if variance is known.  var=1  in the algebraic expressions of the model if variance is known.
}



\references{
Siegmund, D. and Zhang, H.P.  (1994),  "Confidence regions in broken line regression,"  in  _Change-point Problems_,  IMS Lecture Notes - Monograph Series, vol. 23,  eds. E. Carlstein, H. Muller and D. Siegmund,  Hayward, CA: Institute of Mathematical Statistics,  pp. 292-316.

Draper, N.R. and Smith, H.  (1998), _Applied Regression Analysis_, 3rd ed.,  New York:  Wiley.

}

\keyword{ package }

\seealso{
ci, cr, mle, sety, sl
\cr  'sety' help page gives an example simulation test.
\cr  'sl' help page gives details about the evaluation method choices.
}

\examples{
## Data for Patient B from Smith and Cook (1980):

x = c( 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. )
y = c( 37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4 )
sc <- blmr(x,y)

sc$mle()
sc$ci()
sc$cr(.99,"af")
\donttest{sc$sl(6.1,"mc",0.0001)}

## Test for the presence of a changepoint:

sc$sl( x[1] - 1 )



## An example with variance known, in the Normal approximations
## of binomial random variables.
##      Drinking-and-driving surveys 1998-2007 by TIRF, Ottawa,
## with seasonal adjustment based on monthly surveys by CAMH, Toronto,
## log-odds vs year:

x = c( 0.92, 3.25, 4.29, 5.37, 6.37, 7.71, 8.71, 9.71 )
y = c( -1.194, -2.023, -2.285, -1.815, -1.673, -1.444, -1.237, -1.228 )

Sigma = matrix( c(   0.0361, 0, 0, 0, 0, 0, 0, 0,
          0, 0.0218, 0.0129, 0, 0, 0, 0, 0,
          0, 0.0129, 0.0319, 0, 0, 0, 0, 0,
          0, 0, 0, 0.0451, 0.0389, 0, 0, 0,
          0, 0, 0, 0.0389, 0.0445, 0, 0, 0,
          0, 0, 0, 0, 0, 0.0672, 0.0607, 0.0607,
          0, 0, 0, 0, 0, 0.0607, 0.0664, 0.0607,
          0, 0, 0, 0, 0, 0.0607, 0.0607, 0.0662 ) ,  nrow = 8, ncol = 8 )

tirf_camh <- blmr(x,y,Sigma,1,TRUE)

tirf_camh$ci()
tirf_camh$ci(.95,"af")
\donttest{tirf_camh$cr()}



## A model 2 example, from figure 1 of Chiu et al.(2006),
## log(salmon abundance) vs year:
x = c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
18, 19, 20, 21 )
y = c( 2.499, 2.93, 2.938, 2.827, 2.431, 2.841, 3.064, 2.965, 2.937, 2.645,
2.915, 2.711, 2.928, 2.595, 2.122, 2.076, 1.813, 2.446, 1.714, 0.547, 1.298 )
chiu <- blmr(x,y,2)
chiu$ci(0.90)



## An example that shows different confidence regions by CLR and by AF:
x = c( 1, 2, 3, 4, 5, 6, 7, 8 )
y = c( 1.55, 3.2, 6.3, 4.8, 4.3, 4.0, 3.5, 1.8 )
eg <- blmr(x,y)
eg$cr()
eg$cr(.95,"af")

}

