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
		bool Sbad = false;
		Array2D<double> S2D(n,n);
		for (i=0;i<n;i++) for (j=0;j<n;j++) S2D[i][j] = *(Sigma+i*n+j);
		Cholesky<double> pdtest(S2D);
		if ( !pdtest.is_spd() ) {
			Rcout << "Covariance matrix has invalid entries or is not positive-definite:" << endl;
			Rcout << S2D << endl;
			Rcout << "Model's covariance matrix ";
			if (Sig==NULL) Rcout << "defaults to identity matrix."; else Rcout << "not changed.";
			Rcout << endl;
			if (Sig==NULL) set_Sigma();
			Sbad = true;
		}  else  {
			cov_matrix_I = true;
			cov_matrix_diagonal = true;
			for (i=0;i<n;i++)  for (j=0;j<n;j++) {
				*(Sig+i*n+j) = *(Sigma+i*n+j);
				if (i==j && *(Sig+i*n+j)!=1.) cov_matrix_I = false;
				if (i!=j && *(Sig+i*n+j)!=0.) {cov_matrix_I = false; cov_matrix_diagonal = false;}
			}
			get_rS_irS();
		}
		if (Sbad) return;
	}


	if (Xd != NULL) set_x(Xd);

	return;
}



