//
// routines to get various matrices used in Clmbr

#include "lmbr.h"



void  Clmbr::get_Q(void)  const
//  Generate orthogonal matrix 'Q' that annihilates the x-matrix,
//  and  'Q1'  that annihilates the x-matrix except for the 1- and x1- vectors.
{
	Matrix<double>  Q(m,n,0.),  Q1(m1,n,0.);
	int  i, j, k;

	Vector<double>  *beta= NULL,  *alpha= NULL;
	try{
		beta = new (Vector<double>[n*n]);
		alpha = new (Vector<double>[n*n]);
	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 9 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}
	const Vector<double>  dummy( n, 0. );
	for(i=0;i<n;i++)  { beta[i]= dummy;  alpha[i]= dummy; }

	int xcol;
	if(Model==M3) xcol=0; else xcol=1;
	Matrix<double>  X( n, xrank );
	for (i=0;i<n;i++)  for (j=0;j<xrank;j++)  
		if( model_in > 0 )  X[i][j] = *(X_in+i*xrank+j);  else  X[i][j] = *(X_in+(n-1-i)*xrank+j);
	if( model_in < 0 )  for (i=0;i<n;i++)  X[i][xcol] *= -1.;
	X= *pirS * X;

// extend 'X' to a set of independent 'beta' vectors for Gram-Schmidt
// first, reduce columns of X to identify significant rows
	double maxrows[xrank];
	for(j=0;j<xrank;j++) {
		int  maxi, maxk;
		double  maxX= 0.;
		for(k=j;k<xrank;k++)  for(i=0;i<n;i++) 
				if( fabs(X[i][k]) > maxX ) { maxX= fabs(X[i][k]); maxi= i; maxk= k; }
		if(maxk!=j) for(i=0;i<n;i++) 
			{ const double tmp= X[i][j]; X[i][j]= X[i][maxk]; X[i][maxk]= tmp; }
		maxX= X[maxi][j];
		for(i=0;i<n;i++) X[i][j]= X[i][j]/maxX;
		for(k=j+1;k<xrank;k++) {
			const double  Xmi= X[maxi][k];
			for(i=0;i<n;i++)  X[i][k]= X[i][k] - X[i][j]*Xmi;
		}
		maxrows[j] = maxi;
	}

// reset X
	for (i=0;i<n;i++)  for (j=0;j<xrank;j++)  
		if( model_in > 0 )  X[i][j] = *(X_in+i*xrank+j);  else  X[i][j] = *(X_in+(n-1-i)*xrank+j);
	if( model_in < 0 )  for (i=0;i<n;i++)  X[i][xcol] *= -1.;
	X= *pirS * X;

// the first 'beta' vectors are the X-columns, but
// place the  1- and x1- vectors as the final two
	int nm;
	if(Model==M3) nm=1; else nm=2;
	for(i=0;i<xrank-nm;i++)  for(j=0;j<n;j++)  beta[i][j]= X[j][i+nm];
	for(j=0;j<n;j++) beta[i][j]= X[j][0];  
	if(Model!=M3) { i++;  for(j=0;j<n;j++) beta[i][j]= X[j][1]; }

//  fill-out with basis vectors
	int bcol= 0;
	for(i=xrank;i<n;i++)  {
		bool mrow=true;
		while( mrow ) {
			mrow= false;
			for(j=0;j<xrank;j++)  if(bcol==maxrows[j]) bcol++, mrow=true;
		}
		beta[i][bcol++]= 1.;
	}


// confirm that rows are independent
	Array2D<double> M2D(n,n);
	for(i=0;i<n;i++) for(j=0;j<n;j++)  M2D[i][j] = beta[i][j];
	LU<double> rows(M2D);
	if(!rows.isNonsingular()) stop( _("'get_Q' failed") );


// modified Gram-Schmidt
	for (i=0;i<n;i++) {
		alpha[i] = 1./sqrt(beta[i]*beta[i]) * beta[i];
		for (j=0;j<n;j++)  if (fabs(alpha[i][j]) < zero_eq) alpha[i][j] = 0.;
		for (j=i+1;j<n;j++)  beta[j] -= (beta[j]*alpha[i])*alpha[i];
	}


// confirm invertible
	for (i=0;i<n;i++) for (j=0;j<n;j++) M2D[i][j] = alpha[i][j];
	LU<double> lum2(M2D);
	if(!lum2.isNonsingular()) stop( _("'get_Q' failed") );

// confirm orthogonal, eased tolerance
	Matrix<double> Mtest(n,n), Mtest2(n,n);
	for (i=0;i<n;i++) for (j=0;j<n;j++) Mtest[i][j] = M2D[i][j];
	Mtest2= Mtest*transpose(Mtest);
	for (i=0;i<n;i++) for (j=0;j<n;j++) {
		if(i!=j) if(fabs(Mtest2[i][j])>zero_eq*1000)  stop( _("'get_Q' failed") );
		if(i==j) if(fabs(Mtest2[i][j] - 1.)>zero_eq*1000)  stop( _("'get_Q' failed") );
	}


// Q and Q1 = the final 'm' and 'm1' alpha 
	for (i=0;i<m;i++)  for (j=0;j<n;j++)  Q[i][j] = alpha[n-m+i][j];
	for (i=0;i<m1;i++)  for (j=0;j<n;j++)  Q1[i][j] = alpha[n-m1+i][j];

	*pQ = Q;
	*pQ1 = Q1;


	delete[] beta;  delete[] alpha;

	return;
}




void  Clmbr::get_rS_irS(void)  const
//  get matrix rS that is the square-root of  Sigma = cov. matrix of Y 
//  get matrix irS that is the inverse square-root of Sigma  
{
	int i, j;

	Matrix<double>  rSd(n,n,0.),  irSd(n,n,0.);

	double  maxD =0.,  minD =Inf;

	if (cov_matrix_diagonal) {

		for (i=0;i<n;i++) {
			double Di= *(S_in+i*n+i);
			if( model_in < 0 )  Di= *( S_in + (n-1-i)*n + (n-1-i) );
			if( inverse )
				rSd[i][i] = sqrt( Di ),  irSd[i][i] = 1./rSd[i][i];
			else
				irSd[i][i] = sqrt( Di ),  rSd[i][i] = 1./irSd[i][i];
			if( Di > maxD )  maxD= Di;
			if( Di < minD )  minD= Di;
		}
		*prS = rSd;
		*pirS = irSd;

	}  else  {

		Matrix<double>  Se(n,n);
		Array2D<double>  S2D(n,n), Se2D(n,n), Sd2D(n,n);

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

		for (i=0;i<n;i++) {
			double  Di= Sd2D[i][i];
			if( inverse )
				rSd[i][i] = sqrt(Di),  irSd[i][i] = 1./rSd[i][i];
			else
				irSd[i][i] = sqrt(Di),  rSd[i][i] = 1./irSd[i][i];
			if( Di > maxD )  maxD= Di;
			if( Di < minD )  minD= Di;
		}


		*prS = Se*rSd*transpose(Se);	
		*pirS = Se*irSd*transpose(Se);

		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			if(fabs( (*prS)[i][j] ) < zero_eq) (*prS)[i][j]=0.;
			if(fabs( (*pirS)[i][j] ) < zero_eq) (*pirS)[i][j]=0.;
		}

	}

	if( minD/maxD < 1.e-7 )  Rf_warning( _("weights matrix ill-conditioned for 'clr' method") );

	return;
}






void  Clmbr::get_M( Matrix<double> *const pM )  const
//  generate orthogonal matrix 'M' with first row 'gam(theta0)'
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

	Vector<double>  g0 = gam(th0,k0);
	for (j=0;j<m;j++)  if (fabs(g0[j]) < zero_eq)  g0[j] = 0.;


// Pivot on the element of 'g0' with maximum absolute value.
	int  maxj =0;
	double  max =0.;
	for( j=0; j < m; j++ )  if( fabs(g0[j]) > max )  max= fabs(g0[j]),  maxj= j; 
	if( max < zero_eq )  stop( _( "'get_M' failed:  'gam(th0)' is the zero vector" )  );

	beta[0] = g0;
	const Vector<double>  dummy(m,0.);
	for ( i= 1; i < maxj+1; i++ ) {
		beta[i] = dummy;
		beta[i][i-1] = 1.;
	}
	for ( i= maxj+1; i < m; i++ ) {
		beta[i] = dummy;
		beta[i][i] = 1.;
	}


	Array2D<double> M2D(m,m);
	for (j=0;j<m;j++) M2D[0][j]= beta[0][j];
	for (i=1;i<m;i++) for (j=0;j<m;j++) M2D[i][j] = beta[i][j];

// confirm starting rows are independent
	LU<double> lum(M2D);
	if(!lum.isNonsingular()) stop( _("'get_M' failed") );


// modified Gram-Schmidt
	for (i=0;i<m;i++) {
		alpha[i] = 1./sqrt(beta[i]*beta[i]) * beta[i];
		for (j=0;j<m;j++)  if (fabs(alpha[i][j]) < zero_eq) alpha[i][j] = 0.;
		for (j=i+1;j<m;j++) beta[j] -= (beta[j]*alpha[i])*alpha[i];
	}
	alpha[0] = g0;


// confirm 'M' is invertible
	for (i=0;i<m;i++) for (j=0;j<m;j++) M2D[i][j] = alpha[i][j];
	LU<double> lum2(M2D);
	if(!lum2.isNonsingular()) stop( _("'get_M' failed") );

// confirm orthogonal, eased tolerance
	Matrix<double> Mtest(m,m), Mtest2(m,m);
	for (i=0;i<m;i++) for (j=0;j<m;j++) Mtest[i][j] = M2D[i][j];
	Mtest2= Mtest*transpose(Mtest);
	for (i=0;i<m;i++) for (j=0;j<m;j++) {
		if(i!=j) if(fabs(Mtest2[i][j])>zero_eq*1000)  stop( _("'get_M' failed") );
		if(i==j) if(fabs(Mtest2[i][j] - 1.)>zero_eq*1000)  stop( _("'get_M' failed") );
	}


	for (i=0;i<m;i++) for (j=0;j<m;j++)  (*pM)[i][j] = alpha[i][j];

	delete[] beta; delete[] alpha;

	return;
}






