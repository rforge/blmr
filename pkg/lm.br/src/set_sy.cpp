//

#include "lmbr.h"



void  Clmbr::set_sy(double *const irsy, const METHOD met)
// precalculate numbers and vectors based on  irS*y  values 
{
	for (int i=0;i<n;i++)
		if ( isinf(irsy[i]) || isnan(irsy[i]) )  stop( _("invalid value in 'irSy' vector") );



	Vector<double> virsy(n);
	for (int i=0;i<n;i++)  if( model_in > 0 )  virsy[i] = irsy[i];  else  virsy[i] = irsy[n-1-i];
	*psy = *pQ1*virsy;


	*py = *prS*virsy;

	y1 = *psy*(*psig1);
	yx = *psy*(*psigx);
	sysq = *psy*(*psy);  

	*pqy = (*pQ)*virsy;

	qysq = *pqy*(*pqy);

	double  max_gy;
	mle( false, &max_gy );
	omega =  qysq - max_gy ;
	if(omega < zero_eq)  omega= 0.;

	if ( met != INIT ) {
		set_theta0( th0, met );
		set_alpha0( alpha0, met );
		set_SL();
		set_acc();
	}

	return;
}

