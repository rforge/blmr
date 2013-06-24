//


#include "blmr.h"




void Cblmr::set_Sigma( void )
// check and assign flags for covariance matrix, then get its root and inverse root
{
	int i,j;

	Array2D<double> S2D(n,n);
	for (i=0;i<n;i++) for (j=0;j<n;j++) {
		double sij = *(S_in+i*n+j);
		if ( isinf(sij) || isnan(sij) ) stop( _("weights or covariance matrix has invalid entries") );
		if( fabs(sij) < zero_eq )  sij= 0.; 
		S2D[i][j] = sij;
	}
	Cholesky<double> pdtest(S2D);
	if ( !pdtest.is_spd() )  stop( _("weights or covariance matrix is not positive-definite") );

	get_rS_irS();

	if ( px != NULL )  set_x();

	return;
}



