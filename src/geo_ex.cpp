//


#include "blmr.h"






double Cblmr::geo_ex(double *err)
// calculate significance level for theta0 outside data range (Xd[0], Xd[n-1])
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

	double sl, er=0, *per = &er;
	sl = 2*integrate( w, max(aex,2*w), &Cblmr::ipr_ex, 1, per);

	if(err!=0) *err = acc_sl_abs/2 + *per;

	return sl;
}



double Cblmr::ipr_ex(double s, int k)
// integrand for significance level for theta0 outside data range, variance known
{
	if (s < w+zero_eq) return 0.;
	return  ( 2*F(m-1,-w/s) + pow(1.-w*w/s/s,m*0.5-1)*Lgamma/PI ) * Rf_dchisq(s*s,m ,0) * s ;
}


