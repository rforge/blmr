//
// constructor and destructor for class lmbr
//


#include "lmbr.h"




Clmbr::Clmbr( const NumericVector  yR,  const NumericMatrix  xR,  const string  type,  const NumericMatrix  wR,
						const bool  inv,  const bool  var_k )
// constructor
{
	w_in = x_in = y_in = xs = NULL;
	is = NULL;
	rS = irS = Q = tau = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = NULL;
	pq1 = pqx = NULL;
	pmq1 = pm1h = NULL; 
	py = psy = NULL;
	pqy = NULL;


	int i, j;

	n = yR.size();

	xrank = xR.ncol();

// infer model from first column of input 'x' matrix
	bool  xint = true;
	for(i=0;i<n;i++)   if( xR(i,0) != 1 )  xint = false;
	if( xint  &&  xrank==1 )  stop( _("'x' rank invalid") );

	if( type=="LL" ) {
		if( xint )  model_in = 1;  else
			stop( _("'alpha'=0 not supported for type \"LL\"") );
	} else {
		if( type=="TL" ) {
			if( xint )  model_in= 2;  else  model_in= 3;
		}  else  {
			if( type=="LT" ) {
				if( xint )  model_in= -2;  else  model_in= -3; 
			}  else  
				stop( _("'type' must be \"LL\", \"LT\" or \"TL\"") );
		}
	}


	variance_unknown= !var_k;

	inverse= inv;


	bool  cov_matrix_I = true;
	cov_matrix_diagonal = true;
	if( wR.ncol()==1 )  {
		for (i=0;i<n;i++)  if( fabs( wR(i,0) - 1. ) > zero_eq )  cov_matrix_I = false;
	}  else  {
		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			if (i==j &&  fabs( wR(i,j) - 1. )>zero_eq) cov_matrix_I = false;
			if (i!=j && fabs( wR(i,j) )>zero_eq) {cov_matrix_I = false; cov_matrix_diagonal = false;}
		}
	}

	vectorS = false;
	matrixS = false;
	if( !cov_matrix_I )  {
		if( cov_matrix_diagonal )  vectorS = true;  else  matrixS = true;
	}

 
	try{
		y_in = new double[n];  
		x_in = new double[n*xrank]; 
		if( vectorS )  w_in = new double[n];  
		if( matrixS )  w_in = new double[n*n];  

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


// store input values

	for (i=0;i<n;i++) {
		y_in[i] =  yR[i];
		for(j=0; j< xrank; j++)  *(x_in+j*n+i) = xR(i,j);
		if( vectorS )  { if( wR.ncol()==1 )  *(w_in+i) = wR(i,0);  else  *(w_in+i) = wR(i,i); }
		if( matrixS )  for(j=0; j< n; j++)  *(w_in+j*n+i) = wR(i,j);
	}



	initialize();

}



Clmbr::Clmbr( const NumericVector  yR,  const NumericMatrix  xR,  const string  type,  const bool  var_k )
// constructor
{
	w_in = x_in = y_in = xs = NULL;
	is = NULL;
	rS = irS = Q = tau = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = NULL;
	pq1 = pqx = NULL;
	pmq1 = pm1h = NULL; 
	py = psy = NULL;
	pqy = NULL;

	
	int  i;

	n = yR.size();

	xrank = xR.ncol();

// infer model from first column of x matrix
	bool  xint = true;
	for(i=0;i<n;i++)   if( xR(i,0) != 1 )  xint = false;
	if( xint  &&  xrank==1 )  stop( _("invalid x") );

	if( type=="LL" ) {
		if( xint )  model_in = 1;  else
			stop( _("'alpha'=0 not supported for type 'LL'") );
	} else {
		if( type=="TL" ) {
			if( xint )  model_in= 2;  else  model_in= 3;
		}  else  {
			if( type=="LT" ) {
				if( xint )  model_in= -2;  else  model_in= -3; 
			}  else  
				stop( _("'type' must be one of 'LL', 'LT' or 'TL'") );
		}
	}


	variance_unknown= !var_k;


	cov_matrix_diagonal = true;
	vectorS = false;
	matrixS = false;

	try{
		y_in = new double[n];  
		x_in = new double[n*xrank];

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


// store input values

	for (int i=0;i<n;i++) {
		y_in[i] =  yR[i];
		for(int j=0; j< xrank; j++)  *(x_in+j*n+i) = xR(i,j);
	}


	initialize();

}




Clmbr::Clmbr( const Clmbr  &initM )
//copy constructor
{
	w_in = x_in = y_in = xs = NULL;
	is = NULL;
	rS = irS = Q = tau = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = NULL;
	pq1 = pqx = NULL;
	pmq1 = pm1h = NULL; 
	py = psy = NULL;
	pqy = NULL;

	
	model_in = initM.model_in;

	variance_unknown= initM.variance_unknown;

	inverse= initM.inverse;

	n = initM.n;

	xrank = initM.xrank;

	cov_matrix_diagonal = initM.cov_matrix_diagonal;
	vectorS = initM.vectorS;
	matrixS = initM.matrixS;

	try{
		y_in = new double[n];  
		x_in = new double[n*xrank];  
		if( vectorS )  w_in = new double[n];  
		if( matrixS )  w_in = new double[n*n];  

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


// store input values

	for (int i=0;i<n;i++) {
		y_in[i] = (*initM.py)[i];
		for(int j=0; j< xrank; j++)  *(x_in+j*n+i) = *(initM.x_in+j*n+i);
		if( vectorS )  w_in[i] = initM.w_in[i];
		if( matrixS )  for(int j=0; j< n; j++)  *(w_in+j*n+i) = *(initM.w_in+j*n+i);
	}


	initialize();
}




Clmbr::~Clmbr()
// destructor
{

	delete[] w_in;  delete[] x_in;  delete[] y_in;  delete[] xs;
	delete[] is;
	delete[] rS;  delete[] irS;  delete[] Q; delete[] tau;
	delete[] q11;  delete[] qx1;  delete[] qxx;  delete[] ck;  delete[] qff;
	delete[] q10;  delete[] qx0;  delete[] a0;  delete[] b0;
	delete[] f01;  delete[] f0x;
	delete[] B;  delete[] C;
	delete[] px;  delete[] psig1;  delete[] psigx;  delete[] pv1h;  delete[] pxh;
	delete[] nan_m1;  delete[] pnse1;  delete[] pnuse1;  delete[] pusen;
	delete[] nan_m;  delete[] puqe1;  delete[] puqen;  delete[] puqx;
	delete[] ps1;  delete[] psx;
	delete[] pq1;  delete[] pqx;
	delete[] pmq1; delete[] pm1h;
	delete[] py;  delete[] psy;
	delete[] pqy;

}


