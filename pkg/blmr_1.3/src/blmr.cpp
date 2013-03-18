//
// constructors and destructor for class blmr
//


#include "blmr.h"





Cblmr::Cblmr(NumericVector Xda, NumericVector Yda, int n_model , bool var_k )
// constructor
{
	n =Xda.size();

	if (n_model==1) model =M1; else model =M2;

	if (var_k) variance_unknown =false; else variance_unknown =true;

	double *Xtmp, *Ytmp;
	Xtmp = new double[n]; for (int i=0;i<n;i++) Xtmp[i] = Xda[i];
	Ytmp = new double[n]; for (int i=0;i<n;i++) Ytmp[i] = Yda[i];

	initialize( Xtmp,  Ytmp );

	delete[] Xtmp;
	delete[] Ytmp;

}




Cblmr::Cblmr(NumericVector Xda, NumericVector Yda, NumericMatrix S, int n_model , bool var_k )
// constructor
{
	n =Xda.size();

	if (n_model==1) model =M1; else model =M2;

	if (var_k) variance_unknown =false; else variance_unknown =true;

	double *Xtmp, *Ytmp, *Stmp;
	Xtmp = new double[n]; for (int i=0;i<n;i++) Xtmp[i] = Xda[i];
	Ytmp = new double[n]; for (int i=0;i<n;i++) Ytmp[i] = Yda[i];
	Stmp = new double[n*n]; for(int i=0;i<n;i++) for(int j=0;j<n;j++) *(Stmp+i*n+j) = S(i,j);

	initialize( Xtmp,  Ytmp,  Stmp );

	delete[] Xtmp;
	delete[] Ytmp;
	delete[] Stmp;

}




Cblmr::~Cblmr()
// destructor
{

	delete[] C;  delete[] is;
	delete[] Xs;  delete[] Xd;  delete[] Yd;  delete[] Sig;
	delete[] q11;  delete[] qx1;  delete[] qxx;  delete[] ck;  delete[] qff;
	delete[] q10;  delete[] qx0;  delete[] a0;  delete[] b0;
	delete[] px;  delete[] py;  delete[] psy;  delete[] pqy;
	delete[] pmq1;  delete[] pq1;  delete[] pqx;  delete[] ps1;  delete[] psx;  
	delete[] ps1c;  delete[] psxc;
	delete[] psig1;  delete[] psigx;  delete[] pv1h;  delete[] pxh;
	delete[] vund;  delete[] vund2;  delete[] pnse1;  delete[] pnuse1;  
	delete[] puqe1;  delete[] puqen;  delete[] puqx;  delete[] pusen;
	delete[] prS;  delete[] pirS;  delete[] pQ; 

}

