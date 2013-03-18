//


#include "blmr.h"




void Cblmr::set_Sigma(const double *Sigma)
// precalculate working variables that depend only on x values
// and Sigma matrix
// if pointer to Sigma = NULL or 0, assign Sigma = identity matrix
{
	int i,j;

	if (Sigma==0) {

		cov_matrix_I = true;
		cov_matrix_diagonal = true;
		for(i=0;i<n;i++) for(j=0;j<n;j++) *(Sig+i*n+j)=0.;
		for(i=0;i<n;i++) *(Sig+i*n+i) = 1.;
		Matrix<double> Ip(n,n,0.);
		for (i=0;i<n;i++) Ip[i][i] = 1.;		
		*prS = Ip;
		*pirS = Ip;

	}

	else

	{
		Array2D<double> S2D(n,n);
		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			double sij = *(Sigma+i*n+j);
			if ( isinf(sij) || isnan(sij) ) stop( _("weights or covariance matrix has invalid entries") );
			S2D[i][j] = sij;
		}
		Cholesky<double> pdtest(S2D);
		if ( !pdtest.is_spd() )  stop( _("weights or covariance matrix is not positive-definite") );

		cov_matrix_I = true;
		cov_matrix_diagonal = true;
		for (i=0;i<n;i++)  for (j=0;j<n;j++) {
			*(Sig+i*n+j) = *(Sigma+i*n+j);
			if (i==j && *(Sig+i*n+j)!=1.) cov_matrix_I = false;
			if (i!=j && *(Sig+i*n+j)!=0.) {cov_matrix_I = false; cov_matrix_diagonal = false;}
		}
		get_rS_irS();
	}


	if (Xd != NULL) set_x(Xd);

	return;
}



