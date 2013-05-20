\name{blmr}
\alias{blmr-package}
\alias{blmr}

\docType{package}

\title{ Broken Line Model Regression }

\description{
Confidence intervals and confidence regions with exact coverage probabilities for the changepoint in broken line regression.
}


\usage{
blmr(  y,  x,   model =1,  var_known =FALSE,  weights  )
}


\arguments{
  \item{y}{vector of y values}
  \item{x}{vector of x values, non-decreasing}
  \item{model}{ 1, 2 or -2, defined below }
  \item{var_known}{FALSE if variance is unknown, TRUE if known}
  \item{weights}{weights vector or positive-definite matrix, default is the identity matrix}
}


\details{
A broken-line model consists of two straight lines joined at a changepoint.  The versions here are

1:   y = alpha + B' * min( x-theta , 0 ) + B * max( x-theta, 0 ) + e

2:   y = alpha +  B * max( x-theta , 0 ) + e

-2:   y = alpha +  B' * min( x-theta , 0 ) + e

where  e ~ N( 0 , var*weights ).  All parameters are unknown except for 'weights'.  Exact inferences about the changepoint  theta  or  (theta,alpha)  are based on the distribution of its likelihood-ratio statistic, conditional on sufficient statistics for the other parameters.

Function 'blmr' creates a data+model object.  Then methods 'ci' gets confidence intervals, 'cr' confidence regions, and 'sl' significance levels for the changepoint's x-coordinate or (x,y)-coordinates.  Other methods are 'mle' to get maximum likelihood estimates and 'sety' to set new y-values.

Test for a changepoint by the significance level of a postulate value outside the x-values.  Thus,  sl(x[1]-1)  gives an exact significance level of the null hypothesis H0:"single line" versus the alternate hypothesis H1:"broken line".

}


\note{
If variance is known, 'weights' is a variances vector or a variance-covariance matrix, and 'var'=1  in the algebraic expressions of the model above.
}



\references{
Siegmund, D. and Zhang, H.P. (1994),  "Confidence regions in broken line regression,"  in  _Change-point Problems_,  IMS Lecture Notes -- Monograph Series, vol. 23,  eds. E. Carlstein, H. Muller and D. Siegmund,  Hayward, CA: Institute of Mathematical Statistics,  pp. 292-316.
}

\keyword{ package }

\seealso{
ci, cr, mle, sety, sl
\cr  The 'sety' help page gives an example simulations test.
\cr  The 'sl' help page gives details about the evaluation method choices.
}

\examples{
##  Data for Patient B from Smith and Cook (1980),
##  reciprocal of blood creatinine level  vs  day after kidney transplant:
y = c( 37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4 )
x = c( 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. )
sc <- blmr(y,x)

sc$ci()
sc$cr(.90,'af')
sc$sl( x[1] - 1.5 )      ## test for the presence of a changepoint
sc$mle()



##  An example with variance known, in the Normal approximations
##  of binomial random variables.
##       Drinking-and-driving yearly surveys by TIRF, Ottawa,
##  with seasonal adjustment based on monthly surveys by CAMH, Toronto,
##  log-odds vs year:

y = c( -1.194, -2.023, -2.285, -1.815, -1.673, -1.444, -1.237, -1.228 )
x = c( 1998.92, 2001.25, 2002.29, 2003.37, 2004.37, 2005.71, 2006.71, 2007.71 )

Sigma = matrix(  c(   0.0361, 0, 0, 0, 0, 0, 0, 0,
          0, 0.0218, 0.0129, 0, 0, 0, 0, 0,
          0, 0.0129, 0.0319, 0, 0, 0, 0, 0,
          0, 0, 0, 0.0451, 0.0389, 0, 0, 0,
          0, 0, 0, 0.0389, 0.0445, 0, 0, 0,
          0, 0, 0, 0, 0, 0.0672, 0.0607, 0.0607,
          0, 0, 0, 0, 0, 0.0607, 0.0664, 0.0607,
          0, 0, 0, 0, 0, 0.0607, 0.0607, 0.0662 ) ,  nrow = 8, ncol = 8 )

tirf_camh <- blmr(y,x,1,TRUE,Sigma)

tirf_camh$ci()
tirf_camh$ci(.95,'af')
\donttest{tirf_camh$cr()}



##  A model 2 example, from figure 1 of Chiu et al.(2006),
##  log(salmon abundance) vs year:
y = c( 2.50, 2.93, 2.94, 2.83, 2.43, 2.84, 3.06, 2.97, 2.94, 2.65, 2.92,
 2.71, 2.93, 2.60, 2.12, 2.08, 1.81, 2.45, 1.71, 0.55, 1.30 )
x = 1980:2000
chiu <- blmr(y,x,2)
chiu$ci(0.90)



#  A model -2 example with repeat observations, from Robbins, Saxton and Southern (2006). 
#  pig weight-gain in grams-per-day  versus   % dietary isoleucine :
y = c( 533, 452, 433, 418,   686, 543, 648, 530,
       819, 657, 561, 712,   757, 743, 560, 721,
       691, 757, 657, 809,   667, 721, 651, 862  )
x = c( 0.38, 0.38, 0.38, 0.38,  0.42, 0.42, 0.42, 0.42, 
       0.46, 0.46, 0.46, 0.46,  0.50, 0.50, 0.50, 0.50, 
       0.54, 0.54, 0.54, 0.54,  0.58, 0.58, 0.58, 0.58 )

rss <- blmr(y,x,-2)
rss$ci()



##  An example that shows different confidence regions by CLR and by AF:
y = c( 1.55, 3.2, 6.3, 4.8, 4.3, 4.0, 3.5, 1.8 )
x = 1:8
eg <- blmr(y,x)
eg$cr()
eg$cr(.95,"af")

}

