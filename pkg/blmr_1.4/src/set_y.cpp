//
// during initialization, set y values

#include "blmr.h"


void Cblmr::set_y( void )
{
	int i;
	for(i=0;i<n;i++) 
		if ( isinf( *(Y_in+i) ) || isnan( *(Y_in+i) ) ) stop( _("invalid y value") );


	Vector<double>  Y(n,0.), sy(n,0.);

	for (i=0;i<n;i++) {
		Y[i] = Y_in[i];
		if( model_in == -2 )  Y[i]= Y_in[n-1-i];
	}

	sy = (*pirS)*Y;

	double sY[n];
	for (i=0;i<n;i++) {
		sY[i] = sy[i];
		if( model_in == -2 )  sY[i]= sy[n-1-i];		// 'set_sy'  takes sy-values in original order
	}

	set_sy(sY, INIT);

	return;
}

