//
// routines to get various matrices used in Cblmr

#include "blmr.h"



void Cblmr::get_Q(void) const
{
//  generate orthogonal matrix Q that zeroes the 1-vector and the x-vector
//  use 1-hat and x-hat for rows 1 and 2 and identity matrix for the rest, 
//  then apply Gram-Schmidt


//  Identity matrix 'for rest of Q' has two zero rows.
//  If the zero rows are i and j then need to have xi not equal to. xj
//  to have independent columns:  So make is[n-2] and is[n-1] = n-1 the zero rows.


	Vector<double> dummy(n,0.), *beta, *alpha;
	beta = new (Vector<double>[n*n]);
	alpha = new (Vector<double>[n*n]);

	beta[0] = dummy;
	beta[1] = dummy;

	int i,j;
	for (i=2;i<is[ns-2];i++) {
		beta[i] = dummy;
		beta[i][i-2] = 1.;
	}
	for (i=is[ns-2];i<n;i++) {
		beta[i] = dummy;
		beta[i][i-1] = 1.;
	}

	alpha[0] = *pv1h;
	alpha[1] = *pxh;
	for (i=2;i<n;i++) {
		alpha[i] = beta[i];
		for (j=0;j<i;j++) alpha[i] -= (beta[i]*alpha[j])*alpha[j];
		for (j=0;j<n;j++) if (fabs(alpha[i][j]) < zero_eq) alpha[i][j] = 0.;
		alpha[i] = 1./sqrt(alpha[i]*alpha[i])*alpha[i];
	}


	Matrix<double> Q(m,n);
	if (model==M1) for (i=0;i<n-2;i++)  for (j=0;j<n;j++)  Q[i][j] = alpha[i+2][j];
	if (model==M2) for (i=0;i<n-1;i++)  for (j=0;j<n;j++)  Q[i][j] = alpha[i+1][j];

	*pQ = Q;

	delete[] beta;
	delete[] alpha;

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
			rSd[i][i] = sqrt( *(Sig+i*n+i) );
			irSd[i][i] = 1./rSd[i][i];
		}
		*prS = rSd;
		*pirS = irSd;

	}  else  {

		Matrix<double> Se(n,n);
		Array2D<double> S2D(n,n),Se2D(n,n),Sd2D(n,n);

//  get var-cov matrix 

		for (i=0;i<n;i++) for (j=0;j<n;j++) S2D[i][j] = *(Sig+i*n+j);

		Eigenvalue<double> veig(S2D);

		veig.getV(Se2D);

		veig.getD(Sd2D);

		for (i=0;i<n;i++) for (j=0;j<n;j++) Se[i][j] = Se2D[i][j];

		for (i=0;i<n;i++) {
			rSd[i][i] = sqrt(Sd2D[i][i]);
			irSd[i][i] = 1./rSd[i][i];
		}


// the eigenvector matrix is orthogonal so its transpose is its inverse

		*prS = Se*rSd*transpose(Se);
		*pirS = Se*irSd*transpose(Se);

		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			if(fabs( (*prS)[i][j] ) < zero_eq) (*prS)[i][j]=0.;
			if(fabs( (*pirS)[i][j] ) < zero_eq) (*pirS)[i][j]=0.;
		}

	}

	return;
}






void Cblmr::get_M(Matrix<double> *pM) const
//  generate orthogonal matrix M with first row gam(theta0)
{
	int i,j;

// first row of  M  is  gam(th0)
// initial matrix for rest of M has basis vectors and one row zero
// so need   gam(th0)[j] != 0   to have independent columns.
	j=m-1;
	if ( fabs(gam(th0,gk0)[m-1]) < zero_eq ) {
		j=0;
		while(fabs(gam(th0,gk0)[j]) < zero_eq  &&  j<m-1 ) j++;
		if(j==m-1) stop( _( "gam(th0)=0, so cannot get matrix M" )  );
	}

	Vector<double> dummy(m,0.),*beta,*alpha;
	beta = new (Vector<double>[m*m]);
	alpha = new (Vector<double>[m*m]);

	beta[0] = dummy;
	for (i=1;i<j+1;i++) {
		beta[i] = dummy;
		beta[i][i-1] = 1.;
	}
	for (i=j+1;i<m;i++) {
		beta[i] = dummy;
		beta[i][i] = 1.;
	}

	alpha[0] = gam(th0,gk0);
	for (j=0;j<m;j++) if (fabs(alpha[0][j]) < zero_eq) alpha[0][j] = 0.;

	for (i=1;i<m;i++) {
		alpha[i] = beta[i];
		for (j=0;j<i;j++) alpha[i] -= (beta[i]*alpha[j])*alpha[j];
		for (j=0;j<m;j++) if (fabs(alpha[i][j]) < zero_eq) alpha[i][j] = 0.;
		alpha[i] = 1./sqrt(alpha[i]*alpha[i])*alpha[i];
	}

	Array2D<double> M2D(m,m);
	for (i=0;i<m;i++) for (j=0;j<m;j++) M2D[i][j] = alpha[i][j];
	LU<double> lum(M2D);
	if(!lum.isNonsingular()) stop( _("'get_M' failed") );


	for (i=0;i<m;i++) for (j=0;j<m;j++) (*pM)[i][j] = alpha[i][j];

	delete[] beta;
	delete[] alpha;

	return;
}






