



.onAttach <- function(...) {
   library(help=lm.br)$info[[1]] -> version
   version <- version[pmatch("Version",version)]
   um <- strsplit(version," ")[[1]]
   version <- um[nchar(um)>0][2]
   hello <- paste( " lm.br  version ", version, 
      ",  '?lm.br' starts help", sep="" )
   packageStartupMessage( hello )
}






lm.br  <- function( formula, type ="LL", data, subset, weights,
inverse =FALSE, var.known =FALSE, offset, contrasts =NULL, 
na.action, ... )
{
   call <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match( c( "formula", "data", "subset", "weights",
"na.action", "offset" ), names(mf), 0L )
   mf <- mf[c(1L, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   offset <- as.vector(model.offset(mf))

   type <- toupper(type)

   x <- model.matrix( mt, mf, contrasts )

   y <- model.response( mf, "numeric")
   if( NCOL(y) > 1 )  
      stop( "mutiple response vectors not supported" )

   w <- model.weights(mf)
   w_ <- w
   if( is.vector(w) )  { if(inverse)
      { if(any(w==0)) 
           stop( "zero variances not allowed" )  else
         w <- 1/w }
   }  else
   if( is.matrix(w) )  {
      n <- NROW(y)
      if( any( dim(w) != c(n, n) ) )  
         stop( "dim('weights') invalid" )
      eW <- eigen(w, TRUE)
      d <- eW$values
      if( any(d <= 0) ) stop("'weights' not positive-definite")
      A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
      Ainv <- eW$vector %*% diag(d^ifelse(inverse, 0.5, -0.5))
      if(!is.null(offset)) 
         offset <- as.vector( A %*% as.matrix(offset) )
   }

   z  <-  if( is.null(w) )
         lm.fit( x, y, offset=offset, ... )
      else  {
         if( is.vector(w) )
            lm.wfit( x, y, w, offset=offset,  ... )
         else
            lm.fit( A %*% x, A %*% y, offset=offset, ... )
      }


   if( length(z$coef) > 0  &&  !is.na(z$coef[1]) )  {

# 'x1' is the term with a coefficient changepoint
      dn <- colnames(x)
      xint <- if( dn[1] == "(Intercept)" )  TRUE  else  FALSE
      if(xint) x1c=2 else x1c=1
      x1 <- as.vector( x[,x1c] )
      x1nm <- dn[x1c]


# 'xb' is the 'x' input matrix, but with two vectors for 'x1'
#  before and after changepoint,
#  and zeroes for columns that were linearly dependent
      nx <- ncol(x)
      nxb <- ncol(x)+1
      xb <- matrix( 0, nrow(x), nxb )
      xb[,1] <- x[,1]
      for(i in x1c:nx) 
         xb[,i+1] <- if( is.na(z$coef[i]) )  0  else  x[,i]

      bnm <- paste( " ", x1nm, "< theta" )
      bpnm <- paste( " ", x1nm, "> theta" )
      rownames(xb) <- rownames(x)
      colnames(xb) <-  if( xint )
         { if( x1c < nx )  
              c("1-vector", bnm, bpnm, dn[(x1c+1):nx])  
           else  
              c("1-vector", bnm, bpnm) 
         }
         else
         { if( x1c < nx )  
              c(bnm, bpnm, dn[(x1c+1):nx])  
           else
              c(bnm, bpnm) 
         }


#  construct C++ object
#  re-order for 'x1' non-decreasing,  drop rows with  w = 0,  
#  drop columns not independent
      x_ <- x[ order(x1), , drop = FALSE]
      y_ <- y
      if(!is.null(offset)) y_ <- y_ - as.vector(model.offset(mf))
      y_ <- y_[ order(x1) ]
      if( !is.null(w_) )  
         w_ <- if(is.matrix(w_))  w_[ order(x1), order(x1) ]  
                  else  w_[ order(x1) ]

      if( is.vector(w_) )  if( any(w_==0) )  {
         ok <- w_!=0
         w_ <- w_[ok]
         y_ <- y_[ok]
         x_ <- x_[ok, , drop = FALSE]
      }
      if( is.vector(w_) )  w_ <- as.matrix(w_)

      x_fullcol <- x_
      x_ <- x_[ , !is.na(z$coef), drop = FALSE ]


#  loop to drop columns that are dependent at some changepoint value
#  if x-matrix is dependent at  changepoint = 'th'  then  
#  Qf(th)=0  on lower rows,  where  f(th) = max( x1-th, 0 )

      x_dep <- TRUE

      while( x_dep ) {

         obj <- if( is.null(w) )
               new( Cpp_Clmbr, y_, x_, type, var.known )
            else
               new( Cpp_Clmbr, y_, x_, type, w_, inverse, var.known )

         thQfmin <- obj$param()[5]
         xb[ ,x1c] <- if( type=='TL' )  0  else
               { if( is.infinite(thQfmin) )  -1  else  pmin(x1 - thQfmin, 0 ) }
         xb[ ,x1c+1] <- if( type=='LT' )  0  else
               { if ( is.infinite(thQfmin) )  1  else  pmax(x1 - thQfmin, 0 ) }

         z  <-  if( is.null(w) )
               lm.fit( xb, y, offset=offset, ... )
            else  {
               if( is.vector(w) )
                  lm.wfit( xb, y, w, offset=offset,  ... )
               else
                  lm.fit( A %*% xb, A %*% y, offset=offset, ... )
            }

         xbr <- z$rank
         nb <- if(type=='LL')  ncol(x_) + 1  else  ncol(x_)
         if( xbr < nb )  {
            xb[ , is.na(z$coef) ]  <- 0
            ok <- vector( 'logical', nx )
            ok[] <- TRUE
            for(i in (x1c+1):nx)  if( is.na(z$coef[i+1]) )  ok[i]<- FALSE
            x_ <- x_fullcol[ , ok, drop = FALSE ]
         } else
            x_dep <- FALSE

      }


# output
      par <- obj$param()
      if( is.nan( par[1] ) )  {   # case of perfect line
         xb <- x
         if(xint) colnames(xb)[1] <- "  1-vector"
         for(i in 1:nx) if( is.na(z$coef[i+1]) ) xb[,i] <- 0
      }  else  {
         xb[ , x1c ]  <-  if( type=='TL' )  0  else
               { if( is.infinite(par[1]) )  1  else  pmin( x1-par[1], 0 ) }
         xb[ , x1c+1 ]  <-  if( type=='LT' )  0  else
               { if ( is.infinite(par[1]) )  1  else  pmax( x1-par[1], 0 ) }
      }

      if( is.infinite(par[1]) )  {   # can occur for models with alpha=0
         if( type=='LT' )
            { colnames(xb)[1] <- "  (Intercept)";  colnames(xb)[2] <- bnm }
         if( type=='TL' )
            { colnames(xb)[1] <- bpnm;  colnames(xb)[2] <- "  (Intercept)" }
      }


      z  <-  if( is.null(w) )
            lm.fit( xb, y, offset=offset,  ... )
         else  {
            if( is.vector(w) )
               lm.wfit( xb, y, w, offset=offset,  ... )
            else
               lm.fit( A %*% xb, A %*% y, offset=offset, ... )
         }


      if( type=='TL' )  z$coefficients[x1c] <- 0
      if( type=='LT' )  z$coefficients[x1c+1] <- 0

#  add 'theta' to the coefficients
      for(i in ncol(xb):1 )  {
         z$coef[i+1] <- z$coef[i]
         names(z$coef)[i+1] <- names(z$coef)[i]
      }
      z$coef[1] <- par[1]
      names(z$coef)[1] <- "theta"


#  accessor functions
      z$Cpp_obj  <-  obj

      z$ci <- function( CL =0.95, method ="clr" )
         (z$Cpp_obj)$ci( CL, method )

      z$cr <- function( CL =0.95 , method ="clr", incr =NULL )  {
         if( is.null(incr) )
            (z$Cpp_obj)$cr( CL, method )
         else
            (z$Cpp_obj)$cr( CL, method, incr )
      }

      z$sl <- function( theta0, alpha0 =NULL,  method ="clr", 
                  accuracy =0.001, verbose =NULL )
      # overload by if-else statements
      {
         if( !is.null(alpha0) && !is.numeric(alpha0) )  {
            if( is.numeric(method) ) {
               if( is.null(verbose) )
                  if( accuracy==0 | accuracy==1 ) verbose <- accuracy
               accuracy <- method
            }
            method <- alpha0
         }
         if( is.null(alpha0) | !is.numeric(alpha0) )  {
            if( is.null(verbose) )
               (z$Cpp_obj)$sl( theta0, method, accuracy )
            else
               (z$Cpp_obj)$sl( theta0, method, accuracy, verbose )
         } else {
            if( is.null(verbose) )
               (z$Cpp_obj)$sl( theta0, alpha0, method, accuracy )
            else
               (z$Cpp_obj)$sl( theta0, alpha0, method, accuracy, verbose )
         }
      }

      z$mle <- function( )  (z$Cpp_obj)$mle( )

      z$sety <- function( rWy )  {
         if( NCOL(rWy) > 1 | !is.numeric(rWy) )  
            stop( "'rWy' must be a numeric vector" )
         rWy <- drop(rWy)
         if( !is.null(z$offset) )  rWy <- rWy - z$rWoffset
         rWysort <- rWy[ order(z$x1) ]
         if( is.vector(z$weights) )  if( any(z$weights==0) )  {
            wsort <- z$weights[ order(z$x1) ]
            rWysort <- rWysort[ wsort!=0 ]
         }
         (z$Cpp_obj)$sety( rWysort )
      }


#  output list
      p <-  if( is.nan( par[1] ) )  z$rank  else  z$rank + 1
      z$no_of_parameters <- p
      z$df.residual <- nrow(x_) - p
      z$xint <- xint
      z$x1 <- x1
   }

   if( is.vector(w) )  {
      if(!is.null(offset)) z$rWoffset <- sqrt(w) * offset
   }  else
   if( is.matrix(w) )  {
       z$fitted.values <- drop(Ainv %*% z$fitted.values)
       z$residuals <- drop(Ainv %*% z$residuals)
      if(!is.null(offset)) z$rWoffset <- offset
   }

   class(z) <- "lm.br"
   z$call <- call
   z$na.action <- attr(mf, "na.action")
   z$contrasts <- attr(x, "contrasts")
   z$terms <- mt
   z$xlevels <- .getXlevels(mt, mf)
   z$offset <- as.vector(model.offset(mf))
   z$weights <- model.weights(mf)
   z$type <- type
   z
}




print.lm.br  <-  function ( x, digits = max(3L, getOption("digits") - 3L), ... )
{
   cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse ="\n"),
        "\n\n", sep = "")
   type <- x$type
   cat( "Broken-line type:  ", type, "\n\n" )
   if (length(coef(x))>1  && !is.na(x$coef[2])) {
      cat( "Significance Level of hypothesis \"no changepoint\" versus",
              " h. \"one changepoint\"\n" )
      cat("  ")
      mx1 <- max( x$x1 )
      mn1 <- min( x$x1 )
      ai <- (mx1-mn1)/( length( x$x1 ) - 1 )
      if( type=='LT' && x$xint )
         x$sl( round( min( mx1 + ai*1.5 ), 2) )
      else
         if( type=='TL' && !x$xint )
            x$sl( -Inf )
         else
            if( type=='LT' && !x$xint )
               x$sl( Inf )
            else
               x$sl( round( max( mn1 - ai*1.5 ), 2) )
      cat("\n")
      x$ci()
        cat( "Best-fit changepoint and coefficients:\n" )
      print.default( round( x$coef, 5 ) )
      cat("\n")
    }
    else  cat( "No coefficients\n" )
    cat("\n")
    invisible(x)
}





