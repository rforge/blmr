



.onAttach <- function(...) { 
  library(help=lm.br)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste( " lm.br  version ", version, ",  \"?lm.br\" starts help", sep="" )
  packageStartupMessage( hello )
}






lm.br  <-  function  ( formula, type = "LL", data, subset, weights, inverse = FALSE, var.known = FALSE,
    na.action, method = "qr", model = TRUE, x = FALSE, y = FALSE, contrasts = NULL, offset, ... ) 
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c( "formula", "data", "subset", "weights", "na.action", 
        "offset" ), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame") 
        return(mf)
    else if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)

    mt <- attr(mf, "terms")
    y <- model.response( mf, "numeric")
	ynm <- names(mf)[1]
    ny <- NCOL(y)
    if (is.matrix(y) && ny == 1) 
        y <- drop(y)

    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
	if ( is.null(model.weights(mf))  ||  is.vector(model.weights(mf)) )  {
       w <- as.vector( model.weights(mf) )
	}  else  {
		if( is.matrix(model.weights(mf)) )  w <- as.matrix(model.weights(mf))
	}
    if (!is.null(w) && !is.numeric(w)) 
          stop("'weights' must be a numeric vector or matrix")

    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
            3) else numeric(), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
            0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }

    }  else  {

        x <- model.matrix( mt, mf, contrasts )
		if (is.null(n <- nrow(x))) 
		    stop("'x' invalid")
		if (NROW(y) != n) 
		    stop("incompatible dimensions")
		dn <- colnames(x)
		if( dn[1] == "(Intercept)" ) xint= TRUE else xint= FALSE

		if( xint )  x1 <- as.vector( x[,2] )  else  x1 <- as.vector( x[,1] )
		dn1 <- c( " ", dn )
		if( xint )  x1nm <- dn[2]  else  x1nm <- dn[1]
		if( xint ) dn1[1] <- "1-vector"

		x <- as.matrix( x[ order(x1), ] )
		y <- if(is.matrix(y)) y[ order(x1), ]  else  y[ order(x1) ]
		if(!is.null(w))  w <- if(is.matrix(w)) w[ order(x1), order(x1) ]  else  w[ order(x1) ]

		b <- type
		if( b=='LL' ) {
			if( xint ) {  
				model_no=1
				p = 4
			}  else
				stop("'alpha'=0 is not supported for broken line type 'LL'")
		} else {
			if( b=='TL' ) {
				if( xint )  
					{ model_no=2; p=3 }  
				else  
					{ model_no=3; p=2 }
			}  else  {
				if( b=='LT' ) {
					if( xint )  { model_no=-2; p=3 }  else  { model_no=-3; p=2 } 
				}  else  
					stop(gettextf("type = '%s' is not supported",b), domain = NA)
			}
		}
        if (is.null(w)) {
			W <- matrix(0,n,n)
			diag(W) <- 1	
        }  else  {
			if(is.vector(w))  {
				if ( length(w) != n) 
					stop("incompatible dimensions")
				if (any(w < 0 | is.na(w))) 
					stop("missing or negative weights not allowed")
				ok <- w!=0
				w <- w[ok]
				if(is.matrix(y))  y <- y[ok,]  else  y <- y[ok]
				x <- as.matrix( x[ok,] )
				n <- length(w)
				W <- matrix(0,n,n)
				diag(W) <- w
			}  else  {
				if (any(dim(w) != c(n, n))) 
					stop("dim(weights) is not correct")
				W <- w
			}
		}

## drop 'x'-columns that are not independent
		nc <- ncol(x)
		xdrop <- FALSE
		if( qr( x, ... )$rank < nc ) {
			xdrop <- TRUE
			ndrop <- nc - qr(x)$rank
			cdrop <- qr(x)$pivot[ (nc-ndrop+1):nc ]
			x <- as.matrix( x[ , (1:nc)!=cdrop ] )
		}
		if (ncol(x) == 0L) {
		    z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
		        3) else numeric(), residuals = y, fitted.values = 0 * 
		        y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
		        0) else if (is.matrix(y)) nrow(y) else length(y))
		    if (!is.null(offset)) {
		        z$fitted.values <- offset
		        z$residuals <- y - offset
		    }
			class(z) <- "lm.br"
		    z$call <- cl
			z$type <- type
			return(z)
		}
		if(xdrop) if( xint && any(cdrop==2) || !xint && any(cdrop==1) )  stop("'x' invalid")

		if( is.matrix(y) )  { 
			objs <- list()
			for(i in 1:ny) {
				yi <- as.vector( y[ , i ] )
				objs[[i]] <- new( Cpp_Clmbr, yi, x, model_no, W, var.known, inverse )
			}
			names(objs) <- colnames(y)
			objs
		}  else
			obj <- new( Cpp_Clmbr, y, x, model_no, W, var.known, inverse )


# restore input model, but with two new vectors for the two halves of the changepoint term 
# in order to use 'lm.fit' and 'lm.wfit' to calculate the output list
		y <- model.response(mf,'numeric')
		if (is.matrix(y) && ny == 1)  y <- drop(y)
        x <- model.matrix( mt, mf, contrasts )
		n <- NROW(x)
		w <- model.weights( mf )
		if(xdrop) x[ ,cdrop] = 0.	
		if ( is.null(model.weights(mf))  ||  is.vector(model.weights(mf)) )  {
			w <- as.vector( model.weights(mf) )
		}  else  {
			if( is.matrix(model.weights(mf)) )  w <- as.matrix(model.weights(mf))
		}
		xb <- matrix( 0, NROW(x), nc+1 )
		xb[,1] <- x[,1]
		if(xint) xb1=2 else xb1=1
		for(i in xb1:nc) xb[,i+1] <- x[,i]
		if( is.matrix(w) ) {
			eW <- eigen(w, TRUE)
			d <- eW$values
			A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
			Ainv <- eW$vector %*% diag(d^ifelse(inverse, 0.5, -0.5))
		}

		if( is.matrix(y) ) {
			z <- list()
			z$coefficients <- list()
			z$fitted.values <- list()
			z$residuals <- list()
			z$effects <- list()
			for(i in 1:ny)  {
				par <- (objs[[i]])$param()
				if( is.nan( par[1] ) ) {	#perfect line
					xb <- x
				}  else  {
					cp <- as( round(par[1],2), "character" )
					xb[ ,xb1] <- if( b=='TL' )  0  else  pmin( x1 - par[1], 0 )
					xb[ ,xb1+1] <- if( b=='LT' )  0  else  pmax( x1 - par[1], 0 )
					bnm <- paste( "   ", x1nm, "<", cp )	
					bpnm <- paste( "   ", cp, "<", x1nm )
					if( xint ) { 
						dn1[2]<-bnm; dn1[3]<-bpnm;
					}  else  {
						if( is.infinite(par[1]) ) {
							xb[ ,xb1] <- 1
							xb[ ,xb1+1] <- 0
							if( par[1] > 0 )  dn1[2] <- bnm  else  dn1[2] <- bpnm
						}  else
							{ dn1[1]<-bnm; dn1[2]<-bpnm }
					}
				}
				yi <- as.matrix( y[ , i ] )

				if (is.null(w)) 
				    zi <- lm.fit(xb, yi, offset = offset, method = method, ... )
				else 
					if(is.vector(w)) 
						zi <- lm.wfit(xb, yi, w, offset = offset, method = method, ... )
					else  {
						if ( !is.null(offset) )  yi <- yi - offset
						zi <- lm.fit(A %*% xb, A %*% yi, method = method, ... )
						zi$residuals <- drop(Ainv %*% zi$residuals)
						zi$fitted.values <- drop(Ainv %*% zi$fitted.values)
						if ( !is.null(offset) )  zi$fitted.values <- zi$fitted.values + offset
					}

				if( is.infinite(par[1]) ) zi$coefficients[2] <- 0
				if( b=='TL' )  zi$coefficients[xb1] <- 0
				if( b=='LT' )  zi$coefficients[xb1+1] <- 0
				z$coefficients[[i]] <- zi$coefficients
				if( !is.nan(par[1]) )  names( z$coefficients[[i]] ) <- dn1
				z$fitted.values[[i]] <- zi$fitted.values 
				z$residuals[[i]] <- zi$residuals  	
				z$effects[[i]] <- zi$effects
			}
			names(z$coefficients) <- colnames(y)
			names(z$fitted.values) <- colnames(y)
			names(z$residuals) <- colnames(y)
			names(z$effects) <- colnames(y)

		} else  {
			par <- obj$param()
			if( is.nan( par[1] ) ) {
				xb <- x
			}  else  {
				cp <- as( round(par[1],2), "character" )
				xb[ ,xb1] <- if( b=='TL' )  0  else  pmin( x1 - par[1], 0 )
				xb[ ,xb1+1] <- if( b=='LT' )  0  else  pmax( x1 - par[1], 0 )
				bnm <- paste( "   ", x1nm, "<", cp )	
				bpnm <- paste( "   ", cp, "<", x1nm )
				if( xint ) { 
					dn1[2]<-bnm; dn1[3]<-bpnm;
				}  else  {
					if( is.infinite(par[1]) ) {
						xb[ ,xb1] <- 1
						xb[ ,xb1+1] <- 0
						if( par[1] > 0 )  dn1[2] <- bnm  else  dn1[2] <- bpnm
					}  else
						{ dn1[1]<-bnm; dn1[2]<-bpnm }
				}
			}
		    if (is.null(w)) 
		        z <- lm.fit(xb, y, offset = offset, method = method, ... )
		    else 
				if(is.vector(w)) 
					z <- lm.wfit(xb, y, w, offset = offset, method = method, ... )
				else  {
					if ( !is.null(offset) )  y <- y - offset
					z <- lm.fit(A %*% xb, A %*% y, method = method, ... )
					z$residuals <- drop(Ainv %*% z$residuals)
					z$fitted.values <- drop(Ainv %*% z$fitted.values)
					if ( !is.null(offset) )  z$fitted.values <- z$fitted.values + offset
				}
		
			if( !is.nan(par[1]) )  names( z$coefficients ) <- dn1
			if( is.infinite(par[1]) ) z$coefficients[2] <- 0
			if( b=='TL' )  z$coefficients[xb1] <- 0
			if( b=='LT' )  z$coefficients[xb1+1] <- 0
		}

		z$Cpp_obj  <-  if(is.matrix(y))  objs  else  obj

		z$ci <- function(...)  if(is.matrix(y))  { 
							for(i in 1:ny)  {
								cat( ynm, i, "  " )
								(z$Cpp_obj[[i]])$ci(...)
							}
						}  else  
							(z$Cpp_obj)$ci(...)

		z$cr <- function(...)  if(is.matrix(y))  { 
							for(i in 1:ny)  {
								cat( ynm, i, "  " )
								(z$Cpp_obj[[i]])$cr(...)
							}
						}  else  (z$Cpp_obj)$cr(...)

		z$sl <- function(...)  if(is.matrix(y))  { 
							for(i in 1:ny)  {
								cat( ynm, i, "  " )
								(z$Cpp_obj[[i]])$sl(...)
							}
						}  else  (z$Cpp_obj)$sl(...)

		z$mle <- function(...)  if(is.matrix(y))  { 
							for(i in 1:ny)  {
								cat( ynm, i, "  " )
								(z$Cpp_obj[[i]])$mle(...)
							}
						}  else  (z$Cpp_obj)$mle(...)

		z$sety <- function( rWy ) {
							if( NCOL(y) > 1 )  stop("matrix 'y' not supported for 'sety'")
							rWysort <- rWy[ order(z$x1) ]
							(z$Cpp_obj)$sety( rWysort )
						}


##      return 'rank' and 'df' for 'theta' within 'x1' domain, and not perfect line,
## where domain is  [ min(x1), max(x1) ] in 'LL',  [ min(x1), Inf ] in 'LT',  and  [ -Inf, max(x1) ] in 'TL'
## 'rank' and 'df' would be  p-1  and  n-p+1  for 'theta' outside of domain, and 
## the values of a linear model if perfect line
		z$rank <- p
		z$df.residual <- n - p		
		if( is.vector(w) ) z$df.residual = z$df.residual - sum( w==0 )
		z$xint <- xint
		z$x1 <- x1
    }
	
    class(z) <- "lm.br"
	z$type <- type
    z$na.action <- attr(mf, "na.action")
	z$weights <- model.weights(mf)
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- model.matrix( mt, mf, contrasts )
    if (ret.y) 
        z$y <- model.response( mf, "numeric")
    z
}




print.lm.br  <-  function ( x, digits = max(3L, getOption("digits") - 3L), ... ) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
	b <- x$type
	cat( "Broken line type:  ", b, "\n\n" )
    if (length(coef(x))) {
		cat("Significance Level of hypothesis \"no changepoint\" versus h. \"single changepoint\":\n")
		if( (b=='LT' && x$xint) || (b=='TL' && !x$xint) )
			x$sl( max( x$x1 ) + 1.5 )
		else
			x$sl( min( x$x1 ) - 1.5 )
        cat("Coefficients:\n")
		print.default( x$coef )
		cat("\n")
		if( length( x$x1 ) < 30 )  x$ci()  else  x$ci(0.95,'af')
    }
    else  cat("No coefficients\n")
    cat("\n")
    invisible(x)
}





