//

#include "blmr.h"



void  Cblmr::set_sy(double *sy, METHOD met)
// precalculate numbers and vectors based on  irS*y  values 
{
	bool sybad = false;
	for (int i=0;i<n;i++) {
		if ( isinf(sy[i]) || isnan(sy[i]) ) { 
			Rcout << "Invalid (irS*y)[i]  value: " << endl;
			Rcout << sy[i] << endl;
			Rcout << "irSy values not changed." << endl;
			sybad = true;
			break;
		}
	}
	if (sybad) return;

	Vector<double> vsy(n);
	for (int i=0;i<n;i++) vsy[i] = sy[i];
	*psy = vsy;


	*py = *prS*(*psy);
	for (int i=0;i<n;i++) Yd[i] = (*py)[i];

	y1 = *psy*(*psig1);
	yx = *psy*(*psigx);
	sysq = *psy*(*psy);

	*pqy = (*pQ)*(*psy);
	qysq = *pqy*(*pqy);

	omega =  max( qysq - max_mk() , 0. );

	set_theta0(th0,met);

	if (met==INIT || met==GEO2 || met==AF2 || met==MC)  set_alpha0(alpha0,met);

	set_acc();

	return;
}

