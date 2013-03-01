//


#include "blmr.h"




double Cblmr::sl_af(int mode)
// estimate significance level for broken line regression by Approximate F method
{
	double errors0sq = qysq - star_z*star_z;

	int q;  if (th0ex) q = 2; else q = 1;  if (mode==2) q++;

	double sl;
	if (variance_unknown) {
		double Fstat = (m-2)*1./q*fabs(errors0sq/omega - 1.);
		sl =  1. -  Rf_pf(Fstat,q,m-2,1,0) ;
	}  else  {
		double CHIstat = fabs(errors0sq - omega);
		sl =  1. -  Rf_pchisq( CHIstat, q ,1,0) ;
	}

	return sl;
}





double Cblmr::sl_af2(void) 
{

	double sl;
	int q;
	if (th0ex) q = 3; else q = 2;

	if (variance_unknown) {
		double Fstat = (m-2)*1./q*fabs(lambdasq/omega - 1.);
		sl =  1. -  Rf_pf(Fstat,q,m-2,1,0) ;
	}  else  {
		double CHIstat = fabs(lambdasq - omega);
		sl =  1. -  Rf_pchisq( CHIstat, q ,1,0) ;
	}

	return sl;
}
