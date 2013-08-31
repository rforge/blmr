//
// during initialization, set y values

#include "lmbr.h"


void Clmbr::set_y( void )
{
	int i;
	for(i=0;i<n;i++) 
		if ( isinf( *(y_in+i) ) || isnan( *(y_in+i) ) )  stop( _("invalid y value") );

	Vector<double>  y(n,0.),  irSy(n,0.);

	for (i=0;i<n;i++)  if( model_in > 0 )  y[i] = y_in[i];  else  y[i]= y_in[n-1-i];

	irSy = y;
	if( vectorS )  for(i=0;i<n;i++)  irSy[i] =  *( irS + i ) * y[i];
	if( matrixS )  for(i=0;i<n;i++) {
		irSy[i] = 0.;
		for( int k=0; k<n; k++ )  irSy[i] += *( irS + k*n + i ) * y[k];
	} 

	double irsy[n];
	for (i=0;i<n;i++)
		if( model_in > 0 )  irsy[i] = irSy[i];  else  irsy[i]= irSy[n-1-i];		// 'set_sy'  takes sy-values in original order
																				//  i.e. not reversed for LT or LT0 models

	set_sy( irsy, INIT );

	return;
}

