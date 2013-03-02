//

#include "blmr.h"




double Cblmr::F(int k, double arg)
// calculate the F() function defined in K,S&Z eq.(14)
{

	double Fx;
	if (arg <= -1. + zero_eq)
		Fx = 0.;
	else
		if (arg >= 1. - zero_eq)
			Fx = 1.;
		else
			Fx = Rf_pt( arg*sqrt(k/(1-arg*arg)), k, 1, 0) ;
	return Fx;
}




double Cblmr::get_C(const int k) const
{
	int i = 0;
	double c = 1.;

	if ( k % 2 ) {
		i = (k-1)/2;
		while (i>0) { c *= i / (i - 0.5); i--;}
		c /= PI;
	} else {
		i = k/2 - 1;
		while (i>0) { c *= (i + 0.5) / i; i--;}
		c /= 2;
	}

	return c;
}




double Cblmr::fk(int k, double arg)
// calculate the f() function defined in K,S&Z eq.(13)
{
	double gv = 0.;
	if ( fabs(arg) < 1. )  gv = C[n-1-k]*pow( 1.- arg*arg, k/2. - 1.);
	return gv;
}



