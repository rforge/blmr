
% File  man/blmr-internal.Rd
% Part of the blmr package

\name{blmr-internal}

\alias{Cblmr}
\alias{Cblmr-class}
\alias{Rcpp_Cblmr-class}

\docType{class}


\title{Internal blmr class \code{"Cblmr"} }

\description{
Internal S4 class.

Class of data + model objects used for broken line regression.  An object instance stores pre-calculated quantities for a particular data set, so that methods can quickly compute exact-confidence intervals and regions for the changepoint.

}

\details{
  This class is not meant for the user.
}


\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}


\section{Methods}{
  \describe{
    \item{\code{ci}:}{ printout confidence interval }
    \item{\code{cr}:}{ printout confidence region }
    \item{\code{mle}:}{ maximum likelihood estimates of parameters }
    \item{\code{sety}:}{ set new y-values in a Cblmr object }
    \item{\code{sl}:}{ significance level for a postulate changepoint value }
  }
}


\seealso{ blmr }

\examples{
showClass( Cblmr )
}


\keyword{ classes }
\keyword{ internal }

