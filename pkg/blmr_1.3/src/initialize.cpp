//
//

#include "blmr.h"



void Cblmr::initialize( const double *Xdata, const double *Ydata, const double *Sigma)
// initialize variables, printout 
// precalculate working variables in the object
{
	C = NULL; is = NULL;
	Xs = Xd = Yd = Sig = NULL;
	q11 = qx1 = qxx = ck = qff = NULL;
	q10 = qx0 = a0 = b0 = NULL;
	px = py = psy = pqy = NULL;
	pmq1 = pq1 = pqx = ps1 = psx = ps1c = psxc = NULL;
	psig1 = psigx = pv1h = pxh = NULL;
	vund = vund2 = pnse1 = pnuse1 = puqe1 = puqen = puqx = pusen = NULL;
	prS = pirS = pQ = NULL; 


	if ( (n<4 && model==M1) || (n<3 && model==M2) ) {
		stop( _("invalid number of data observations, must be at least 4 in M1 or 3 in M2") );
	}


	Rcout << endl << _( "initializing:" ) << endl;


	if(model==M1) { m= n-2; k1= 1; } else { m= n-1; k1= 0; }

	C = new double[4];
	C[0]= get_C(n-1); C[1]= get_C(n-2); C[2]= get_C(n-3); C[3]= get_C(n-4);


	Sig = new double[n*n]; 
	prS = new (Matrix<double>[n*n]);
	pirS = new (Matrix<double>[n*n]);

	set_Sigma(Sigma);


	vund = new (Vector<double>[m]);
	vund2 = new (Vector<double>[n]);
	pnse1 = new (Vector<double>[n]);
	pnuse1 = new (Vector<double>[n]);
	pusen = new (Vector<double>[n]);
	puqe1 = new (Vector<double>[m]);
	puqen = new (Vector<double>[m]);
	puqx = new (Vector<double>[m]);

	Xd = new double[n];
	px = new (Vector<double>[n]);
	pv1h = new (Vector<double>[n]);
	pxh = new (Vector<double>[n]);
	psig1 = new (Vector<double>[n]);
	psigx = new (Vector<double>[n]);
	pQ = new (Matrix<double>[m*n]);
	is = new int[n];

	set_x(Xdata);

	set_xs();


	Yd = new double[n];
	py = new (Vector<double>[n]);
	psy = new (Vector<double>[n]);
	pqy = new (Vector<double>[m]);

	th0 = Xs[1];
	alpha0 = Yd[1];
	z = w = 0.;

	set_y(Ydata);


	ah = 0; old_th = inf;
	prev_th=inf; a_low=0; a_high=0;

	prev_SL=-1; 
	set_SL();

	srand( (unsigned)time( NULL ) );


//printout
	Rcout << setprecision(6);
	rel_print_eps = 0.00001;

	int i,j,n_model;
	if(model==M1) n_model=1; else n_model=2;
	Rcout << _( "model " ) << n_model;
	if (variance_unknown) Rcout << _("  with variance unknown");  else  Rcout << _("  with variance known");
	Rcout << endl << "x = "; for (i=0;i<n;i++) Rcout << " " << Xdata[i];
	Rcout << endl << "y = "; for (i=0;i<n;i++) Rcout << " " << Ydata[i];
	Rcout << endl;
	if ( cov_matrix_diagonal  &&  !cov_matrix_I ) {
		if(variance_unknown) Rcout << _("weights = ");
		if(!variance_unknown) Rcout << _("variances = ");
		for (i=0;i<n;i++) Rcout << " " << *(Sig+i*n+i);
		Rcout << endl;
	} else  {
		if(!cov_matrix_I ) {
			Matrix<double> S(n,n,0.);
			for (i=0;i<n;i++) for (j=0;j<n;j++) S[i][j] = *(Sig+i*n+j);
			if(variance_unknown) Rcout << _("weights matrix = "); else  Rcout << _("variance-covariance matrix = ");
			Rcout << S;
		}
		Rcout << endl;
	}



	return;
}

