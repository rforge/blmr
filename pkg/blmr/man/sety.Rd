\name{sety,Cblmr-method}
\alias{sety,Cblmr-method}
\alias{sety}

\title{Set y-Values}

\description{
  Set the y-values in a Cblmr object.
}

\usage{
\S4method{sety}{Cblmr}( irSy )
}

\arguments{
\item{irSy}{vector of y-values, pre-multiplied by the inverse-square-root of S}
}

\note{
  The pre-multiplied vector is more convenient as input during simulation tests.  This vector is simply the y-vector if S is the identity matrix.
}


\examples{
## Data for Patient B from Smith and Cook (1980):

x= c(1., 2., 3., 4., 5., 6., 7., 8., 9., 10.)
y= c(36.0, 45.5, 50.0, 60.0, 73.3, 71.0, 66.7, 60.0, 30.5, 18.3)
sc <- blmr(x,y)
sc$ci()
newy= c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8,
19.4)
sc$sety(newy)
sc$ci()



## An example simulation test:

simtest <- function( x, theta, alpha, Bp, B, var, N =10000) {

   y <- x
   mod <- blmr(x,y)     ## construct a Cblmr object with arbitrary y values.

   cat("\n   Monte Carlo simulation test:\n")
   cat("for model   theta=",theta," alpha=",alpha," B'=",Bp," B=",B," var=",var,"\n")
   cat("with observations at  x=",x,"\n\n")
   cat("     no. of     coverage frequency of the 0.95-confidence interval by\n")
   cat("   iterations               CLR                    AF\n")
   flush.console()

   sigma = sqrt(var)
   n = length(x)
   countCLR = countAF = 0

   for( i in seq(1,N) ) {

      err <- rnorm(n,0,sigma)
      for( j in seq(1,n) )
         y[j] = alpha + Bp*min(x[j]-theta,0) + B*max(x[j]-theta,0) + err[j]
      mod$sety(y)
      stest = mod$sl(theta,"clr",.0001,FALSE)
      if(stest>0.05) countCLR = countCLR + 1
      stest = mod$sl(theta,"af",.0001,FALSE)
      if(stest>0.05) countAF = countAF + 1

      if(i/1000 - floor(i/1000) == 0) {
         cat( format(i,width=10),
              format(countCLR/i,digits=4,nsmall=4,width=22),
              format(countAF/i,digits=4,nsmall=4,width=22), "\n" )
         flush.console()
      }
   }

   cat("\n")
}


x = c( 1.0, 1.1, 1.3, 1.7, 2.4, 3.9, 5.7, 7.6, 8.4, 8.6 )

simtest(x,3,0,-1,0.5,1)

}

