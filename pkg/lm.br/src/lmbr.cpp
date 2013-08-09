//
// constructor and destructor for class lmbr
//


#include "lmbr.h"




Clmbr::Clmbr(const NumericVector yR, const NumericMatrix xR, const int model, const NumericMatrix SigmaR,
						const bool var_k, const bool inv )
// constructor
{
	S_in = X_in = Y_in = xs = NULL;
	is = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = NULL;
	C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = ps1c = psxc = NULL;
	pq1 = pqx = NULL;
	pmq1 = NULL; 
	pm1h = NULL;
	py = psy = NULL;
	pqy = NULL;
	prS = pirS = pQ = NULL;

	

	model_in = model;

	if( var_k )  variance_unknown= false;  else  variance_unknown= true;

	if( inv )  inverse= true;  else  inverse= false;

	n = yR.size();

	xrank = xR.ncol();


	try{
		Y_in = new double[n];  
		X_in = new double[n*xrank];  
		S_in = new double[n*n];  

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 1 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


// store input values

	for (int i=0;i<n;i++) {
		Y_in[i] =  yR[i];
		for(int j=0; j< xrank; j++)  *(X_in+i*xrank+j) = xR(i,j);
		for(int j=0; j< n; j++)  *(S_in+i*n+j) = SigmaR(i,j);
	}



	initialize();

}




Clmbr::Clmbr(const Clmbr &initM)
//copy constructor
{
	S_in = X_in = Y_in = xs = NULL;
	is = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
    f01 = f0x = NULL;
    B = NULL;
	C = NULL;
	px = psig1 = psigx = pv1h = pxh = NULL;
	nan_m1 = pnse1 = pnuse1 = pusen = NULL;
	nan_m = puqe1 = puqen = puqx = NULL;
	ps1 = psx = ps1c = psxc = NULL;
	pq1 = pqx = NULL;
	pmq1 = NULL; 
	pm1h = NULL;
	py = psy = NULL;
	pqy = NULL;
	prS = pirS = pQ = NULL;

	
	model_in = initM.model_in;

	if (!initM.variance_unknown) variance_unknown= false;  else  variance_unknown= true;

	if (!initM.inverse) inverse= false;  else  inverse= true;

	n = initM.n;

	xrank = initM.xrank;


	try{
		Y_in = new double[n];  
		X_in = new double[n*xrank];  
		S_in = new double[n*n];  

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 2 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


// store input values

	for (int i=0;i<n;i++) {
		Y_in[i] = (*initM.py)[i];
		for(int j=0; j< xrank; j++)  *(X_in+i*xrank+j) = *(initM.X_in+i*xrank+j);
		for(int j=0; j <n; j++)  *(S_in+i*n+j) = *(initM.S_in+i*n+j);
	}


	initialize();
}




Clmbr::~Clmbr()
// destructor
{

	delete[] S_in;  delete[] X_in;  delete[] Y_in;  delete[] xs;
	delete[] is;
	delete[] q11;  delete[] qx1;  delete[] qxx;  delete[] ck;  delete[] qff;
	delete[] q10;  delete[] qx0;  delete[] a0;  delete[] b0;
	delete[] f01;  delete[] f0x;
	delete[] B;
	delete[] C;
	delete[] px;  delete[] psig1;  delete[] psigx;  delete[] pv1h;  delete[] pxh;
	delete[] nan_m1;  delete[] pnse1;  delete[] pnuse1;  delete[] pusen;
	delete[] nan_m;  delete[] puqe1;  delete[] puqen;  delete[] puqx;
	delete[] ps1;  delete[] psx;  delete[] ps1c;  delete[] psxc;
	delete[] pq1;  delete[] pqx;
	delete[] pmq1; delete[] pm1h;
	delete[] py;  delete[] psy;
	delete[] pqy;
	delete[] prS;  delete[] pirS;  delete[] pQ;

}


