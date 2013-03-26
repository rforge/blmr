//

#include "blmr.h"


void Cblmr::set_y(const double *Ydata)
{
	int i;
	for(i=0;i<n;i++) 
		if ( isinf( *(Ydata+i) ) || isnan( *(Ydata+i) ) ) stop( _("invalid y value") );

	Vector<double> Y(n);
	for (i=0;i<n;i++) Y[i] = Ydata[i];
	Y = (*pirS)*Y;

	double sY[n];
	for (i=0;i<n;i++) sY[i] = Y[i];

	set_sy(sY);

	return;
}
