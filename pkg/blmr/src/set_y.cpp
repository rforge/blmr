//

#include "blmr.h"


void Cblmr::set_y(const double *Ydata)
{
	int i;
	bool Ybad = false;
	for (i=0;i<n;i++) {
		if ( isinf( *(Ydata+i) ) || isnan( *(Ydata+i) ) ) { 
			Rcout << "Invalid y[" << i << "] value " << *(Ydata+i) << endl;
			Rcout << "Y values not ";
			if (Yd==NULL) Rcout << "assigned."; else Rcout << "changed.";
			Rcout << endl;
			Ybad = true;
			break;
		}
	}
	if (Ybad) return;

	Vector<double> Y(n);
	for (i=0;i<n;i++) Y[i] = Ydata[i];
	Y = (*pirS)*Y;

	double sY[n];
	for (i=0;i<n;i++) sY[i] = Y[i];

	set_sy(sY);

	return;
}
