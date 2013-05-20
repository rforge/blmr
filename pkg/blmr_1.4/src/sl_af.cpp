//


#include "blmr.h"




double Cblmr::sl_af(int mode)  const
// estimate significance level for broken line regression by Approximate F method
{
	const double errors0sq = qysq - prime_z*prime_z;

	int q;  if (th0ex) q = 2; else q = 1;  if (mode==2) q++;

	double sL;
	if (variance_unknown) {
		const double Fstat = (m-2)*1./q*fabs(errors0sq/omega - 1.);
		sL =  1. -  Rf_pf(Fstat,q,m-2,1,0) ;
	}  else  {
		const double CHIstat = fabs(errors0sq - omega);
		sL =  1. -  Rf_pchisq( CHIstat, q ,1,0) ;
	}

	return sL;
}





double Cblmr::sl_af2(void)  const
{

	double sL;
	int q;
	if (th0ex) q = 3; else q = 2;

	if (variance_unknown) {
		const double Fstat = (m-2)*1./q*fabs(lambdasq/omega - 1.);
		sL =  1. -  Rf_pf(Fstat,q,m-2,1,0) ;
	}  else  {
		const double CHIstat = fabs(lambdasq - omega);
		sL =  1. -  Rf_pchisq( CHIstat, q ,1,0) ;
	}

	return sL;
}
