//


#include "blmr.h"




void  Cblmr::set_alpha0(const double a_0, const METHOD met)
// precalculate numbers and vectors based on alpha0, used in routines
{
	if ( isinf(a_0) || isnan(a_0) )  stop( _("invalid 'alpha0' value") );

	alpha0 = a_0;

	Vector<double>  star_y(n);
	star_y = *psy - alpha0*(*psig1);

	if (Model==M1) {

		if (th0ex) {

			Vector<double>  pf0(n);
			pf0 = *psigx - th0*(*psig1);
			const double  yf0 = star_y*pf0;

			lambdasq = star_y*star_y - yf0*yf0/(pf0*pf0);
			lambda = sqrt( max(lambdasq,0.) );

		}  else  {

			Vector<double>  gsm0(n),  gfr0(n),  gbar0(n);
			gsm0 = gsm(th0,k0);  for(int i=0; i<m; i++)  if( fabs(gsm0[i]) < zero_eq )  gsm0[i] = 0.;
			gfr0 = gfr(th0,k0);  for(int i=0; i<m; i++)  if( fabs(gfr0[i]) < zero_eq )  gfr0[i] = 0.;
			gbar0 = gbar(th0,k0);  for(int i=0; i<m; i++)  if( fabs(gbar0[i]) < zero_eq )  gbar0[i] = 0.;

			const double ysm0 = star_y*gsm0;
			const double yfr0 = star_y*gfr0;
			lambdasq = star_y*star_y - ysm0*ysm0 - yfr0*yfr0;
			lambda = sqrt( max(lambdasq,0.) );

			if ( ! (met==AF || met==AF2) ) {
				const double ybar0 = star_y*gbar0;
				const double c3 = (*pv1h*gfr0)*(*pxh*gsm0) - (*pv1h*gsm0)*(*pxh*gfr0);
				c1 = -lambda*c3;
				c2 = prime_z + ybar0*c3;
			}
		}
	}


	if (Model==M2) {

		if (th0ex) {

			lambdasq = star_y*star_y;
			lambda = sqrt( max( lambdasq, 0.) );

		}  else  {

			Vector<double>  gfr0(n),  gbar0(n);
			gfr0 = gfr(th0,k0);  for(int i=0; i<m; i++)  if( fabs(gfr0[i]) < zero_eq )  gfr0[i] = 0.;
			gbar0 = gbar_prime(th0,k0);  for(int i=0; i<m; i++)  if( fabs(gbar0[i]) < zero_eq )  gbar0[i] = 0.;

			const double yfr0 = star_y*gfr0;
			lambdasq = star_y*star_y - yfr0*yfr0;
			lambda = sqrt( max( lambdasq, 0.) );

			if ( ! (met==AF || met==AF2) ) {
				const double ybar0 = star_y*gbar0;
				const double c3 = *pv1h*gfr0;
				c1 = -lambda*c3;
				c2 = prime_z + ybar0*c3;
			}
		}
	}


	if(lambdasq < zero_eq)  lambdasq= 0.; 
	if(lambda < zero_eq)  lambda= 0.;


	if(omega==0) c= 1;  else  c= sqrt(  max( 1 - omega/lambdasq , 0. )  );


	return;
}
