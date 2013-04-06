


"blmr" <- function(x, y, S, model =1, var_known =FALSE) {


# initial checks

  if (!is.vector(x, mode="numeric")) stop("'x' must be a numeric vector")
  if (!is.vector(y, mode="numeric")) stop("'y' must be a numeric vector")
  if (length(x)!=length(y)) stop("'x' and 'y' must have same length")

 
# overload by if-else statements

  if (missing(S)) { 

     out <- new(Cblmr,x,y,model,var_known)

  }  else  { 

     if (!is.matrix(S)  && length(S)==1 ) {

		if ( !missing(var_known) ) stop("'S' must be a numeric vector or matrix")

        if ( !missing(model) ) var_known <- model

        model <- S

        if ( length(model)>1 || is.logical(model) || !is.numeric(model) || (model!=1 && 		model!=2 && model!=-2) ) 
              stop("'model' must be 1, 2 or -2")
        if ( length(var_known)>1
           || (!is.logical(var_known) && !is.numeric(var_known))
           || (var_known!=TRUE && var_known!=FALSE && var_known!=1 && var_known!=0) )
              stop("'var_known' must be TRUE, FALSE, 1 or 0")

        out <- new(Cblmr,x,y,model,var_known)

     } else { 

		if (!is.numeric(S)) stop ("'S' must be a numeric vector or matrix")

		if ( is.matrix(S) )  {
		   if (nrow(S)!=ncol(S) || nrow(S)!=length(x))  
			  stop("matrix 'S' must be n x n where n=length(x)") 
		} else {
		   n <- length(x)
		   if (length(S)!=n) stop("vector 'S' must have same length as 'x'")
		   Stmp <- matrix(0,nrow=n,ncol=n)
		   for(i in 1:n) Stmp[i,i] = S[i]
		   S <- Stmp
		}

        if ( length(model)>1 || is.logical(model) || !is.numeric(model) || (model!=1 && 		model!=2 && model!=-2) ) 
              stop("'model' must be 1, 2 or -2")
        if ( length(var_known)>1
           || (!is.logical(var_known) && !is.numeric(var_known))
           || (var_known!=TRUE && var_known!=FALSE && var_known!=1 && var_known!=0) )
              stop("'var_known' must be TRUE, FALSE, 1 or 0")

		out <- new(Cblmr,x,y,S,model,var_known) 

     }

  }


  return(out)

}


