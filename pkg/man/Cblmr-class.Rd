\name{Cblmr-class}
\Rdversion{1.1}
\docType{class}
\alias{Cblmr}
\alias{Rcpp_Cblmr-class}


\title{Class \code{"Cblmr"}}

\description{
  Class of data + model objects used for broken line regression.  An object instance stores pre-calculated quantities for a particular data set, so that methods can quickly compute exact-confidence intervals and regions for the changepoint.
}

\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}

\seealso{ blmr }

\examples{
showClass( Cblmr )
}

\keyword{classes}

\section{Methods}{
  \describe{
    \item{\code{ci}:}{ printout confidence interval }
    \item{\code{cr}:}{ printout confidence region }
    \item{\code{mle}:}{ maximum likelihood estimates of parameters }
    \item{\code{sety}:}{ set new y-values in a Cblmr object }
    \item{\code{sl}:}{ significance level for a postulate changepoint value }
  }
}
