//


#include "blmr.h"






double Cblmr::geo_ex(double *err)
// calculate significance level for theta0 outside data range (xs[0], xs[ns-1])
// using Knowles and Siegmund's geometric formulae
{
	if(err!=0) *err=0.;
	if (variance_unknown) return geo_vu_ex();
					else  return geo_vk_ex(err);
}



double Cblmr::geo_vu_ex(void)
// case variance unknown
{
	return  2*F(m-1,-w) + pow(1.-w*w, m*0.5-1)*Lgamma/PI;
}



double Cblmr::geo_vk_ex(double *err)
// case variance known
{
	if(err!=0) *err = 0;

	return  2*( Rf_pnorm5(-w, 0,1,1,0) + Lgamma/sqrt(2*PI) * Rf_dnorm4(w, 0,1,0) );
}



