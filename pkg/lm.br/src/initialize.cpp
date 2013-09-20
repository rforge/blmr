//
//

#include "lmbr.h"



void Clmbr::initialize( void )
// initialize variables, printout 
// precalculate working variables in the object
{
	int  i, j;


// echo input
	const bool  echo_input = false;

	if( echo_input )  {
		Rcout << endl << _( "initializing:" ) << endl;
		Rcout << _( "model " ) << model_in;
		if (variance_unknown)  Rcout << _("  with variance unknown");  else  Rcout << _("  with variance known");
		Rcout << endl << "y = "; for (i=0;i<n;i++) Rcout << " " << y_in[i];
		Rcout << endl << _("x matrix:  "); 
		Matrix<double> X(n,xrank,0.);
		for (i=0;i<n;i++) for (j=0;j<xrank;j++) X[i][j] = *(x_in+j*n+i);
		Rcout << X;
		Rcout << endl;
		if ( vectorS )  {
			if(variance_unknown || !inverse)  Rcout << _("weights vector:  ");  else  Rcout << _("variances:  ");
			for (i=0;i<n;i++) Rcout << " " << *(w_in+i);
			Rcout << endl << endl;
		} 
		if ( matrixS )  {
			if(variance_unknown || !inverse)  Rcout << _("weights matrix:  ");  else  Rcout << _("variance-covariance matrix:  ");
			Matrix<double> W(n,n,0.);
			for (i=0;i<n;i++) for (j=0;j<n;j++) W[i][j] = *(w_in+j*n+i);
			Rcout << W;
			Rcout << endl;
		}
	}


	if(model_in==1)  Model= M1;  else 
		if(model_in==2 || model_in== -2)  Model= M2;			// treat Model -2 as M2 internally
			else  if(model_in== 3 || model_in== -3)  Model= M3;	// treat Model -3 as M3 internally


	if(Model==M1)  { m= n-2-(xrank-2);  m1= m+2;  k1=  1; } 
	if(Model==M2)  { m= n-1-(xrank-2);  m1= m+1;  k1=  0; }
	if(Model==M3)  { m= n-(xrank-1);    m1= m;    k1= -1; }


// for multivariate models, set cov_matrix flag to 'non-diagonal'  
// to invoke more general routines 
	if( m1 < n )  cov_matrix_diagonal = false;



//  'S' = 'Sigma'  is the covariate or covariance matrix,  such that
//   errors  ~ N( 0, var * Sigma ) ,   Sigma = inverse( 'weights' )

	if( vectorS || matrixS )  {

		try {
			if(vectorS)  {
				 rS = new double[n];
				irS = new double[n];
			}  else  {
				 rS = new double[n*n];
				irS = new double[n*n];
			}

		} catch( bad_alloc &ex ) {
			Rcout << _("message: ") << ex.what() << endl;
			stop( _("memory allocation failed") );
		}

		set_Sigma();
	}



	try{
		px = new (Vector<double>[n]);
		pv1h = new (Vector<double>[m1]);
		pxh = new (Vector<double>[m1]);
		psig1 = new (Vector<double>[m1]);
		psigx = new (Vector<double>[m1]);
		Q = new double[n*xrank];
		tau = new double[xrank];
		is = new int[n];

		nan_m1 = new (Vector<double>[m1]);
		nan_m = new (Vector<double>[m]);
		pnse1 = new (Vector<double>[m1]);
		pnuse1 = new (Vector<double>[m1]);
		pusen = new (Vector<double>[m1]);
		puqe1 = new (Vector<double>[m]);
		puqen = new (Vector<double>[m]);
		puqx = new (Vector<double>[m]);

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


	set_x();		// allocates memory for arrays with size 'ns' or 'ns*m'


	try{
		py = new (Vector<double>[n]);
		psy = new (Vector<double>[m1]);
		pqy = new (Vector<double>[m]);

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << ex.what() << endl;
		stop( _("memory allocation failed") );
	}

	set_y();		// calls 'set_sy'

	const double th_0 = xs[1];
	th0 = th_0 + 1;
	th0MC = xs[ns-1] + 1;
	set_theta0(th_0, INIT);		// initializes 'z', 'w'

	const double a0 = (*py)[ is[1] ];
	alpha0 = a0 + 1;
	set_alpha0(a0, INIT);

	prev_SL= -1; 
	set_SL();


	set_acc();


	try{
		C = new double[3];

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << ex.what() << endl;
		stop( _("memory allocation failed") );
	}
	C[0]= get_C(m-2); C[1]= get_C(m-1); C[2]= get_C(m); 


	old_th = prev_th =  Inf; 
	ah = a_low = a_high =  0; 


	return;
}

