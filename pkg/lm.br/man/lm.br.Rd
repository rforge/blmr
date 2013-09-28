
\name{'lm.br'}
\alias{lm.br-package}
\alias{lm.br}
\alias{print.lm.br}


\docType{package}

\title{ Fit a Linear Model with a Breakpoint }

\description{
Exact significance tests for a coefficient changepoint in
linear or multivariate linear regression, assuming continuity.
Confidence intervals and confidence regions with exact coverage
probabilities for the changepoint.
}


\usage{
lm.br( formula, type = "LL", data, subset, weights, inverse = FALSE,
var.known = FALSE, offset, contrasts, na.action, ... )
}


\arguments{
  \item{formula}{a formula expression as for regression models, of
the form \code{response ~ predictors}, see \code{\link{formula}} 
for more}
  \item{type}{ "LL", "LT" or "TL" which stand for line-line,
line-threshold or threshold-line, defined below }
  \item{data}{an optional data frame that defines variables in
\code{formula} }
  \item{subset}{expression saying which subset of the data to use}
  \item{weights}{vector or positive-definite matrix}
  \item{inverse}{if TRUE, 'weights' specifies the inverse of the
weights vector or matrix, as for a covariance matrix}
  \item{var.known}{TRUE or FALSE}
  \item{offset}{a constant vector to be subtracted from the
responses vector}
  \item{contrasts}{an optional list, see 'contrasts.arg' in
\code{\link{model.matrix.default}}  }
  \item{na.action}{a function to filter missing data}
  \item{\dots}{ other arguments to \code{\link{lm.fit}} or
\code{\link{lm.wfit}} }
}


\details{
A broken-line model consists of two straight lines joined at a
changepoint.  Three versions are

LL:   y = alpha +  B * min( x-theta , 0 ) + B'* max( x-theta, 0 )
+ e

LT:   y = alpha +  B * min( x-theta , 0 ) + e

TL:   y = alpha +  B'* max( x-theta , 0 ) + e

where  e ~ N( 0 , var * inv(weights) ).  All parameters are
unknown except for 'weights'.  The LT and TL models
omit 'alpha' if the formula is without intercept, such as 'y~x+0'.

The same models apply for a multivariate formula such as  'y ~
x1 + x2 + ... + xn'  where 'alpha' becomes the coefficient of the
"1"-vector and 'theta' becomes the changepoint for the cofficient 
of the first predictor term, 'x1'.

Exact inferences about the changepoint  'theta'  or
'(theta,alpha)'  are based on the distribution of its
likelihood-ratio statistic, conditional on sufficient statistics
for the other parameters.

Test for the presence of a changepoint by the significance level
of a postulate value outside of the x-values.  Thus,
'sl( min(x1) - 1 )'  would give the exact significance level of the null
hypothesis "single line" versus the alternate hypothesis
"broken line," in the LL model.

}


\value{
'lm.br' returns a list that includes a C++ object with accessor
functions.

Function 'ci' gets confidence intervals, 'cr' confidence regions,
and 'sl' significance levels for the changepoint's x-coordinate or
(x,y)-coordinates.

Other functions are 'mle' to get maximum likelihood estimates
and 'sety' to set new y-values.

The returned list also contains the components of an 'lm' list
including 'coefficients', 'fitted.values' and 'residuals'.

}


\note{
If variance is known,  'weights'  is the inverse of the variances
vector or variance-covariance matrix, and  'var'=1  in the
algebraic expressions of the model above.
}


\references{
Siegmund, D. and Zhang, H.P. (1994),  "Confidence regions in
broken line regression,"  in  _Change-point Problems_,  IMS
Lecture Notes -- Monograph Series, vol. 23,  eds. E. Carlstein, H.
Muller and D. Siegmund,  Hayward, CA: Institute of Mathematical
Statistics,  pp. 292-316.
}

\keyword{ package }

\seealso{
\code{\link{ci}}, \code{\link{cr}}, \code{\link{sl}}, 
\code{\link{mle}}, \code{\link{sety}}

The 'sety' help page gives an example simulations test.

vignette( "\href{../doc/lm.br.pdf}{lm.br}" )

\code{\link{lm}} describes some of the arguments and output values.
}

\examples{
library( lm.br )

#  Data for Patient B from Smith and Cook (1980),
#  reciprocal of blood creatinine L/micromol  versus
#  day after kidney transplant:
crea <- c( 37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4 )
day <- c( 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. )
sc <- lm.br( crea ~ day )
sc$mle()
sc$ci()
sc$cr(.90,'af')
sc$sl( crea[1] - 1.5 )      # test for the presence of a changepoint



#  A 'TL' example, from figure 1 of  Chiu et al. (2006),
#  log(salmon abundance) vs year:
salmon <- c( 2.50, 2.93, 2.94, 2.83, 2.43, 2.84, 3.06, 2.97, 2.94, 2.65,
 2.92, 2.71, 2.93, 2.60, 2.12, 2.08, 1.81, 2.45, 1.71, 0.55, 1.30 )
year <- 1980:2000
chiu <- lm.br( salmon ~ year, 'tl' )
chiu$ci(0.90)



#  A multivariate example, using one of R's included datasets,
#  automobile miles-per-gallon versus weight and horsepower:

lm.br( mpg ~ wt + hp, data = mtcars )



#  An example with variance known, for the Normal approximations
#  of binomial random variables.
#       Drinking-and-driving yearly surveys by TIRF, Ottawa,
#  with seasonal adjustment based on monthly surveys by CAMH, Toronto,
#  log-odds vs year:

log_odds <- c( -1.194, -2.023, -2.285, -1.815, -1.673, -1.444, -1.237, -1.228 )
year <- c(1998.92, 2001.25, 2002.29, 2003.37, 2004.37, 2005.71, 2006.71, 2007.71)

VarCov <- matrix(  c(   0.0361, 0, 0, 0, 0, 0, 0, 0,
          0, 0.0218, 0.0129, 0, 0, 0, 0, 0,
          0, 0.0129, 0.0319, 0, 0, 0, 0, 0,
          0, 0, 0, 0.0451, 0.0389, 0, 0, 0,
          0, 0, 0, 0.0389, 0.0445, 0, 0, 0,
          0, 0, 0, 0, 0, 0.0672, 0.0607, 0.0607,
          0, 0, 0, 0, 0, 0.0607, 0.0664, 0.0607,
          0, 0, 0, 0, 0, 0.0607, 0.0607, 0.0662 ) ,  nrow = 8, ncol = 8 )

dd <- lm.br( log_odds ~ year, w = VarCov, inv = TRUE, var.known = TRUE )

dd$ci()
dd$ci(.95,'af')
\donttest{dd$cr()}



#  An example that shows different confidence regions
#  by conditional likelihood-ratio inference and by approximate-F inference:
y <- c( 1.55, 3.2, 6.3, 4.8, 4.3, 4.0, 3.5, 1.8 )
x <- 1:8
eg <- lm.br(y~x)
\donttest{eg$cr()}
eg$cr( m = 'af' )

}

