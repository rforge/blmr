//

#include "blmr.h"



void  Cblmr::set_sy(double *const sy, const METHOD met)
// precalculate numbers and vectors based on  irS*y  values 
{
	for (int i=0;i<n;i++)
		if ( isinf(sy[i]) || isnan(sy[i]) )  stop( _("invalid value in 'irSy' vector") );



	Vector<double> vsy(n);
	for (int i=0;i<n;i++) vsy[i] = sy[i];
	if( model_in == -2 )  for (int i=0;i<n;i++) vsy[i] = sy[n-1-i];
	*psy = vsy;


	*py = *prS*(*psy);

	y1 = *psy*(*psig1);
	yx = *psy*(*psigx);
	sysq = *psy*(*psy);

	*pqy = (*pQ)*(*psy);

	qysq = *pqy*(*pqy);

	omega =  max( qysq - max_mk() , 0. );
	if(omega < zero_eq) omega= 0.;

	if ( met != INIT ) {
		set_theta0(th0,met);
		set_alpha0(alpha0,met);
		set_SL();
		set_acc();
	}

	return;
}

