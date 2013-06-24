//
// routines to get various matrices used in Cblmr

#include "blmr.h"



void Cblmr::get_Q(void) const
//  Generate orthogonal matrix Q that zeroes the 1-vector and the x-vector.
//  Use 1-hat and x-hat for rows 1 and 2 and identity matrix for the rest, 
//  then apply Gram-Schmidt.
{
	Matrix<double> Q(m,n,0.);


	if(Model== M3) {
		for(int i=0;i<n;i++) Q[i][i]= 1.;
		*pQ = Q;
		return;
	}


	Vector<double>  *beta= NULL,  *alpha= NULL;
	try{
		beta = new (Vector<double>[n*n]);
		alpha = new (Vector<double>[n*n]);
	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 9 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


	const Vector<double>  dummy( n, 0. );

	beta[0] = dummy;
	beta[1] = dummy;


//  Basis vectors as initial 'rest of Q' has two zero columns, say i and j.
//  The rows will be independent if  x[i]  is not equal to  x[j],
//  so choose  i = 0  and  j = n-1  for the zero columns.

	int i,j;
	for ( i = 2; i < n; i++ ) {
		beta[i] = dummy;
		beta[i][i-1] = 1.;
	}


// confirm that rows are independent
	Array2D<double> rows2D(n,n);
	for(j=0;j<n;j++) rows2D[0][j] = (*pv1h)[j];
	for(j=0;j<n;j++) rows2D[1][j] = (*pxh)[j];
	for(i=2;i<n;i++) for(j=0;j<n;j++) rows2D[i][j] = beta[i][j];
	LU<double> rows(rows2D);
	if(!rows.isNonsingular()) stop( _("'get_Q' failed") );


// Gram-Schmidt
	alpha[0] = *pv1h;
	alpha[1] = *pxh;

	for (i=2;i<n;i++) {
		alpha[i] = beta[i];
		for (j=0;j<i;j++) alpha[i] -= (beta[i]*alpha[j])*alpha[j];
		for (j=0;j<n;j++) if (fabs(alpha[i][j]) < zero_eq) alpha[i][j] = 0.;
		alpha[i] = 1./sqrt(alpha[i]*alpha[i])*alpha[i];
	}


	if (Model==M1) for (i=0;i<n-2;i++)  for (j=0;j<n;j++)  Q[i][j] = alpha[i+2][j];
	if (Model==M2) for (i=0;i<n-1;i++)  for (j=0;j<n;j++)  Q[i][j] = alpha[i+1][j];

	delete[] beta;  delete[] alpha;

	*pQ = Q;
	return;
}




void Cblmr::get_rS_irS(void) const
//  get matrix rS that is the square-root of  Sigma = cov matrix of Y 
//  get matrix irS that is the inverse square-root of Sigma  
{
	int i, j;

	Matrix<double>  rSd(n,n,0.), irSd(n,n,0.);

	if (cov_matrix_diagonal) {

		for (i=0;i<n;i++) {
			double Sii= *(S_in+i*n+i);
			if( model_in < 0 )  Sii= *( S_in + (n-1-i)*n + (n-1-i) );
			rSd[i][i] = sqrt( Sii );
			irSd[i][i] = 1./rSd[i][i];
		}
		*prS = rSd;
		*pirS = irSd;

	}  else  {

		Matrix<double>  Se(n,n);
		Array2D<double>  S2D(n,n), Se2D(n,n), Sd2D(n,n);

//  get var-cov matrix 

		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			double Sij= *(S_in+i*n+j);
			if( model_in < 0 )  Sij= *( S_in + (n-1-i)*n + (n-1-j) );
			if( fabs(Sij) < zero_eq )  Sij= 0.;
			S2D[i][j] = Sij;
		}

		Eigenvalue<double> Seig(S2D);

		Seig.getV(Se2D);

		Seig.getD(Sd2D);

		for (i=0;i<n;i++) for (j=0;j<n;j++)  Se[i][j] = Se2D[i][j];

		double  maxD =0., minD =Inf;
		for (i=0;i<n;i++) {
			double  Di= Sd2D[i][i];
			rSd[i][i] = sqrt(Di);
			irSd[i][i] = 1./rSd[i][i];
			if( Di > maxD )  maxD= Di;
			if( Di < minD )  minD= Di;
		}

		if( minD/maxD < 1.e-7 )  Rf_warning( _("weights matrix ill-conditioned for 'clr' method") );


		*prS = Se*rSd*transpose(Se);	// Se is orthogonal so its transpose is its inverse
		*pirS = Se*irSd*transpose(Se);

		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			if(fabs( (*prS)[i][j] ) < zero_eq) (*prS)[i][j]=0.;
			if(fabs( (*pirS)[i][j] ) < zero_eq) (*pirS)[i][j]=0.;
		}

	}

	return;
}






void Cblmr::get_M( Matrix<double> *const pM ) const
//  generate orthogonal matrix M with first row gam(theta0)
{
	Vector<double>  *beta= NULL,  *alpha= NULL;
	try{
		beta = new (Vector<double>[m*m]);
		alpha = new (Vector<double>[m*m]);
	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 10 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


	int i,j;

	alpha[0] = gam(th0,k0);
	for (j=0;j<m;j++)  if (fabs(alpha[0][j]) < zero_eq)  alpha[0][j] = 0.;
//Rcout << "y  " << *py << endl;
//Rcout << "th0 k0 alph0  " << th0 << " " << k0 << "   " << alpha[0] << endl;

	const Vector<double>  dummy(m,0.);
	beta[0] = dummy;

// At start of Gram-Schmidt,  first row of  'M'  is  gam(th0),  other rows are basis vectors, 
// leaving one column 'j' zero,  so need   gam(th0)[j] != 0   to have independent rows.
// Choose  'j'  such that  gam(th0)[j]  is the element with maximum absolute value.
	int maxj =0;
	double  max =0.;
	for( j=0; j < m; j++ )  if( fabs(alpha[0][j]) > max )  max= fabs(alpha[0][j]),  maxj= j; 
	if( max < zero_eq )  stop( _( "gam(th0) is the zero vector, so cannot get matrix M" )  );

	for ( i= 1; i < maxj+1; i++ ) {
		beta[i] = dummy;
		beta[i][i-1] = 1.;
	}
	for ( i= maxj+1; i < m; i++ ) {
		beta[i] = dummy;
		beta[i][i] = 1.;
	}


// confirm rows are independent
	Array2D<double> M2D(m,m);
	for (j=0;j<m;j++) M2D[0][j]= alpha[0][j];
	for (i=1;i<m;i++) for (j=0;j<m;j++) M2D[i][j] = beta[i][j];
//Rcout << "M2D  " << M2D << endl;
	LU<double> lum(M2D);
	if(!lum.isNonsingular()) stop( _("'get_M' failed") );


	for (i=1;i<m;i++) {
		alpha[i] = beta[i];
		for (j=0;j<i;j++) alpha[i] -= (beta[i]*alpha[j])*alpha[j];
		for (j=0;j<m;j++) if (fabs(alpha[i][j]) < zero_eq) alpha[i][j] = 0.;
		alpha[i] = 1./sqrt(alpha[i]*alpha[i])*alpha[i];
	}


// confirm 'M' is nonsingular = invertible
	for (i=0;i<m;i++) for (j=0;j<m;j++) M2D[i][j] = alpha[i][j];
//Rcout << "M2D  " << M2D << endl;
	LU<double> lum2(M2D);
	if(!lum2.isNonsingular()) stop( _("'get_M' failed") );


// confirm orthogonal, eased tolerance
	Matrix<double> Mtest(m,m), Mtest2(m,m);
	for (i=0;i<m;i++) for (j=0;j<m;j++) Mtest[i][j] = alpha[i][j];
	Mtest2= Mtest*transpose(Mtest);
//Rcout << "Mtest2  " << Mtest2 << endl;
	for (i=0;i<m;i++) for (j=0;j<m;j++) {
		if(i!=j) if(fabs(Mtest2[i][j])>zero_eq*1000)  stop( _("'get_M' failed") );
		if(i==j) if(fabs(Mtest2[i][j] - 1.)>zero_eq*1000)  stop( _("'get_M' failed") );
	}

	for (i=0;i<m;i++) for (j=0;j<m;j++)  (*pM)[i][j] = alpha[i][j];

	delete[] beta; delete[] alpha;

	return;
}






