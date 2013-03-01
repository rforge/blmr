


"blmr" <- function(x, y, S, model =1, var_known =FALSE) {

  require(blmr)

# INITIALIZE
  if (!is.vector(x, mode="numeric")) stop("x must be a numeric vector")
  if (!is.vector(y, mode="numeric")) stop("y must be a numeric vector")
  if (length(x)!=length(y)) stop("x and y must have same length")
  
  if (missing(S)) { 

     if ( length(model)>1 || (model!=1 && model!=2) ) 
              stop("model must be 1 or 2")
     if ( length(var_known)>1
           || (!is.logical(var_known) && !is.numeric(var_known))
           || (var_known!=TRUE && var_known!=FALSE && var_known!=1 && var_known!=0) )              stop("var_known must be TRUE or FALSE or 1 or 0")

     out <- new(Cblmr,x,y,model,var_known)

  }  else  { 

     if (!is.matrix(S)) {

        if (!missing(model)) var_known <- model

        model <- S

        if ( length(model)>1 || (model!=1 && model!=2) ) 
              stop("model must be 1 or 2")
        if ( length(var_known)>1
           || (!is.logical(var_known) && !is.numeric(var_known))
           || (var_known!=TRUE && var_known!=FALSE && var_known!=1 && var_known!=0) )              stop("var_known must be TRUE or FALSE or 1 or 0")

        out <- new(Cblmr,x,y,model,var_known)

     } else { 

	if (!is.numeric(S)) stop ("matrix S must be numeric")
	if (nrow(S)!=ncol(S)) stop("matrix S must be square")
	if (nrow(S)!=length(x)) stop("matrix S must be n x n where n=length(x)")

        if ( length(model)>1 || (model!=1 && model!=2) ) 
              stop("model must be 1 or 2")
        if ( length(var_known)>1
           || (!is.logical(var_known) && !is.numeric(var_known))
           || (var_known!=TRUE && var_known!=FALSE && var_known!=1 && var_known!=0) )              stop("var_known must be TRUE or FALSE or 1 or 0")

	out <- new(Cblmr,x,y,S,model,var_known) 

     }

  }


  return(out)

}


