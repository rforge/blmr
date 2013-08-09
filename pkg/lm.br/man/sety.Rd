\name{sety,Cpp_Clmbr-method}
\alias{sety,Cpp_Clmbr-method}
\alias{sety}

\title{Set y-Values}

\description{
  Reset the y-values in the S4 C++ object.  
}

\usage{
\S4method{sety}{Cpp_Clmbr}( rWy )
}

\arguments{
\item{rWy}{vector of y-values, pre-multiplied by the square-root of 'weights'}
}

\details{
This method changes the y-values for the accessor functions 'sl', 'ci', 'cr' and 'mle'.  The output list of '$coefficients' and '$fitted' and so on, however, retains its values from the original call.
}

\note{
The pre-multiplied vector is more convenient as input during simulation tests.  This vector is simply the y-vector if 'weights' is omitted.  

The square-root of a vector 'W' is the vector with the square root of each element.  The square-root of a matrix 'W' here is the matrix 'rW' such that rW*rW = W (a stricter definition than  rW*transpose(rW) = W).
}


\examples{
##  Data for Patient B from Smith and Cook (1980):
library(lm.br)
y= c(36.0, 45.5, 50.0, 60.0, 73.3, 71.0, 66.7, 60.0, 30.5, 18.3)
x= 1:10
sc <- lm.br(y~x)
sc$ci()
newy= c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
sc$sety(newy)
sc$ci()


## An example simulation test:

simtest <- function( x, theta, alpha, B, Bp, var, N =10000) {

   y <- x
   mod <- lm.br(y~x)     ## construct an 'lm.br' object with arbitrary y values.

   cat("\n   Monte Carlo simulation test:\n")
   cat("for model   theta=",theta," alpha=",alpha," B=",B," B'=",Bp," var=",var,"\n")
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
         y[j] = alpha + B*min(x[j]-theta,0) + Bp*max(x[j]-theta,0) + err[j]
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

