//
// during initialization, set y values

#include "lmbr.h"


void Clmbr::set_y( void )
{
	int i;
	for(i=0;i<n;i++) 
		if ( isinf( *(Y_in+i) ) || isnan( *(Y_in+i) ) ) stop( _("invalid y value") );


	Vector<double>  Y(n,0.),  irsy(n,0.);

	for (i=0;i<n;i++)  if( model_in > 0 )  Y[i] = Y_in[i];  else  Y[i]= Y_in[n-1-i];

	irsy = *pirS * Y;

	double irsY[n];
	for (i=0;i<n;i++)
		if( model_in > 0 )  irsY[i] = irsy[i];  else  irsY[i]= irsy[n-1-i];		// 'set_sy'  takes sy-values in original order

	set_sy( irsY, INIT );

	return;
}

