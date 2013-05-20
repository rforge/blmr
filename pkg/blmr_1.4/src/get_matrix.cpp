//
// routines to get various matrices used in Cblmr

#include "blmr.h"



void Cblmr::get_Q(void) const
//  generate orthogonal matrix Q that zeroes the 1-vector and the x-vector
//  use 1-hat and x-hat for rows 1 and 2 and identity matrix for the rest, 
//  then apply Gram-Schmidt
{
	Vector<double>  *beta= NULL,  *alpha= NULL;
	try{
		beta = new (Vector<double>[n*n]);
		alpha = new (Vector<double>[n*n]);
	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 10 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


	const Vector<double>  dummy(n,0.);

	beta[0] = dummy;
	beta[1] = dummy;


//  Basis vectors as initial 'rest of Q' has two zero columns, say i and j.
//  The rows will be independent if  x[i]  is not equal to. x[j].
//  Therefore, make  i_d[n_d-2]  and  i_d[n_d-1] = n-1  the zero columns.

	int i,j;
	for (i=2;i<i_d[n_d-2]+2;i++) {
		beta[i] = dummy;
		beta[i][i-2] = 1.;
	}
	for (i=i_d[n_d-2]+2;i<i_d[n_d-1]+1;i++) {
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


	Matrix<double> Q(m,n);
	if (Model==M1) for (i=0;i<n-2;i++)  for (j=0;j<n;j++)  Q[i][j] = alpha[i+2][j];
	if (Model==M2) for (i=0;i<n-1;i++)  for (j=0;j<n;j++)  Q[i][j] = alpha[i+1][j];

	*pQ = Q;

	delete[] beta;  delete[] alpha;

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
			if( model_in == -2 )  Sii= *( S_in + (n-1-i)*n + (n-1-i) );
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
			if( model_in == -2 )  Sij= *( S_in + (n-1-i)*n + (n-1-j) );
			if( fabs(Sij) < zero_eq )  Sij= 0.;
			S2D[i][j] = Sij;
		}

		Eigenvalue<double> Seig(S2D);

		Seig.getV(Se2D);

		Seig.getD(Sd2D);

		for (i=0;i<n;i++) for (j=0;j<n;j++)  Se[i][j] = Se2D[i][j];

		for (i=0;i<n;i++) {
			rSd[i][i] = sqrt(Sd2D[i][i]);
			irSd[i][i] = 1./rSd[i][i];
		}


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
		Rcout << _("message: ") << 11 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}


// first row of  'M'  is  gam(th0)
// starting rows for Gram-Schmidt are basis vectors, leaving one column 'j' zero, 
// so need   gam(th0)[j] != 0   to have independent rows.

	int i,j;
	j=0;
	while( j < m  &&  fabs(gam(th0,k0)[j]) < zero_eq )  j++;
	if(j==m) stop( _( "gam(th0) is the zero vector, so cannot get matrix M" )  );

	const Vector<double> dummy(m,0.);

	beta[0] = dummy;
	for (i=1;i<j+1;i++) {
		beta[i] = dummy;
		beta[i][i-1] = 1.;
	}
	for (i=j+1;i<m;i++) {
		beta[i] = dummy;
		beta[i][i] = 1.;
	}

	alpha[0] = gam(th0,k0);
	for (j=0;j<m;j++) if (fabs(alpha[0][j]) < zero_eq) alpha[0][j] = 0.;


// confirm rows are independent
	Array2D<double> M2D(m,m);
	for (j=0;j<m;j++) M2D[0][j]= alpha[0][j];
	for (i=1;i<m;i++) for (j=0;j<m;j++) M2D[i][j] = beta[i][j];
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
	LU<double> lum2(M2D);
	if(!lum2.isNonsingular()) stop( _("'get_M' failed") );


// confirm orthogonal, eased tolerance
	Matrix<double> Mtest(m,m), Mtest2(m,m);
	for (i=0;i<m;i++) for (j=0;j<m;j++) Mtest[i][j] = alpha[i][j];
	Mtest2= Mtest*transpose(Mtest);
	for (i=0;i<m;i++) for (j=0;j<m;j++) {
		if(i!=j) if(fabs(Mtest2[i][j])>zero_eq*100)  stop( _("'get_M' failed") );
		if(i==j) if(fabs(Mtest2[i][j] - 1.)>zero_eq*100)  stop( _("'get_M' failed") );
	}

	for (i=0;i<m;i++) for (j=0;j<m;j++) (*pM)[i][j] = alpha[i][j];

	delete[] beta; delete[] alpha;

	return;
}






