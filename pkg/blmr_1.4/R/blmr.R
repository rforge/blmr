


"blmr" <- function( y, x, model =1, var_known =FALSE, weights ) {


# initial checks

	if (!is.vector(x, mode="numeric")) stop("'x' must be a numeric vector")
	if (!is.vector(y, mode="numeric")) stop("'y' must be a numeric vector")
	if (length(x)!=length(y)) stop("'x' and 'y' must have same length")

	if ( length(model)>1 || is.logical(model) || !is.numeric(model) || (model!=1 && model!=2 && model!=-2) ) 
		  stop("'model' must be 1, 2 or -2")
	if ( length(var_known)>1
	   || (!is.logical(var_known) && !is.numeric(var_known))
	   || (var_known!=TRUE && var_known!=FALSE && var_known!=1 && var_known!=0) )
		  stop("'var_known' must be TRUE, FALSE, 1 or 0")

 
	if ( missing( weights )  )  { 

		n <- length(x)
		W <- matrix(0,n,n)
		for(i in 1:n) W[i,i]= 1

		out <- new(Cblmr,y,x,model,var_known,W)

	}  else  { 

		W <- weights

		if (!is.numeric(W)) stop ("'weights' must be a numeric vector or matrix")

		if ( is.matrix(W) )  {
			if (nrow(W)!=ncol(W) || nrow(W)!=length(x))  
				stop("matrix 'weights' must be  n x n  where  n = length(x)") 
		} else {
			n <- length(x)
			if (length(W)!=n) stop("vector 'weights' must have same length as 'x'")
			Wtmp <- matrix(0,nrow=n,ncol=n)
			for(i in 1:n) Wtmp[i,i] = W[i]
			W <- Wtmp
		}

		out <- new(Cblmr,y,x,model,var_known,W) 

	}


  return(out)

}


