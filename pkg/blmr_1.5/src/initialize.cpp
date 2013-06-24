//
//

#include "blmr.h"



void Cblmr::initialize( void )
// initialize variables, printout 
// precalculate working variables in the object
{
	int i,j;


// set  'cov_matrix_I'  and  'cov_matrix_diagonal'  flags

	cov_matrix_I = true;
	cov_matrix_diagonal = true;
	for (i=0;i<n;i++) for (j=0;j<n;j++) {
		const double sij= *(S_in+i*n+j);
		if (i==j && sij!=1.) cov_matrix_I = false;
		if (i!=j && sij!=0.) {cov_matrix_I = false; cov_matrix_diagonal = false;}
	}


// echo input

	Rcout << endl << _( "initializing:" ) << endl;
	Rcout << _( "model " ) << model_in;
	if (variance_unknown) Rcout << _("  with variance unknown");  else  Rcout << _("  with variance known");
	Rcout << endl << "y = "; for (i=0;i<n;i++) Rcout << " " << Y_in[i];
	Rcout << endl << "x = "; for (i=0;i<n;i++) Rcout << " " << X_in[i];
	Rcout << endl;
	if ( cov_matrix_diagonal  &&  !cov_matrix_I ) {
		if(variance_unknown) Rcout << _("weights:  ");
		if(!variance_unknown) Rcout << _("variances:  ");
		for (i=0;i<n;i++) Rcout << " " << *(S_in+i*n+i);
		Rcout << endl;
	} else  {
		if(!cov_matrix_I ) {
			if(variance_unknown) Rcout << _("weights matrix:  "); else  Rcout << _("variance-covariance matrix:  ");
			Matrix<double> S(n,n,0.);
			for (i=0;i<n;i++) for (j=0;j<n;j++) S[i][j] = *(S_in+i*n+j);
			Rcout << S;
		}
		Rcout << endl;
	}

	if(model_in==1)  Model= M1;  else  {
		if(model_in==2 || model_in== -2)  Model= M2;		// treat Model -2 as Model 2 internally
			else  if(model_in== 3 || model_in== -3)  Model= M3;
	}

	if(Model==M1)  { m= n-2; k1= 1; } 
	if(Model==M2)  { m= n-1; k1= 0; }
	if(Model==M3)  { m= n; k1= -1; }

Rcout << "i here1" << endl;
	try {
		prS = new (Matrix<double>[n*n]);
		pirS = new (Matrix<double>[n*n]);

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 3 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}

	set_Sigma();

Rcout << "i here2" << endl;


	try{
		px = new (Vector<double>[n]);
		pv1h = new (Vector<double>[n]);
		pxh = new (Vector<double>[n]);
		psig1 = new (Vector<double>[n]);
		psigx = new (Vector<double>[n]);
		pQ = new (Matrix<double>[m*n]);
		is = new int[n];

		nan_n = new (Vector<double>[n]);
		nan_m = new (Vector<double>[m]);
		pnse1 = new (Vector<double>[n]);
		pnuse1 = new (Vector<double>[n]);
		pusen = new (Vector<double>[n]);
		puqe1 = new (Vector<double>[m]);
		puqen = new (Vector<double>[m]);
		puqx = new (Vector<double>[m]);
		pm1h = new (Vector<double>[m]);

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 4 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


	set_x();		// allocates memory for arrays with size 'ns' or 'ns*m'

Rcout << "i here3" << endl;

	try{
		py = new (Vector<double>[n]);
		psy = new (Vector<double>[n]);
		pqy = new (Vector<double>[m]);

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 6 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}

	set_y();		// calls 'set_sy'
Rcout << "i here4" << endl;

	const double th_0 = xs[1];
	set_theta0(th_0, INIT);	// initializes 'z', 'w'
Rcout << "i here5" << endl;

	const double a0 = (*py)[ is[1] ];
	set_alpha0(a0, INIT);
Rcout << "i here6" << endl;

	prev_SL= -1; 
	set_SL();
Rcout << "i here7" << endl;

	const int  digits= 6;
	Rcout << setprecision( digits );
	rel_print_eps =  pow( 10., -(digits-1) );

	set_acc();
Rcout << "i here8" << endl;


	try{
		C = new double[4];

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 7 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}
	C[0]= get_C(n-1); C[1]= get_C(n-2); C[2]= get_C(n-3); C[3]= get_C(n-4);


	old_th = prev_th =  Inf; 
	ah = a_low = a_high =  0; 

Rcout << "i here9" << endl;

	return;
}

