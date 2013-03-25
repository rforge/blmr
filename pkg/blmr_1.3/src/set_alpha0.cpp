//


#include "blmr.h"




void  Cblmr::set_alpha0(double a_0, METHOD met)
// precalculate numbers and vectors based on alpha0, used in routines
{
	if ( isinf(a_0) || isnan(a_0) )  stop( _("invalid 'alpha0' value") );

	alpha0 = a_0;

	Vector<double> star_y(n);
	star_y = *psy - alpha0*(*psig1);

	if (model==M1) {

		if (th0ex) {

			Vector<double> pf0(n);
			pf0 = *psigx - th0*(*psig1);
			double yf0 = star_y*pf0;

			lambdasq = star_y*star_y - yf0*yf0/(pf0*pf0);
			if (!variance_unknown) lambda = sqrt(lambdasq);

		}  else  {

			Vector<double> gfr0(n), gsm0(n), gbar0(n);
			gsm0 = gsm(th0,gk0);
			gfr0 = gfr(th0,gk0);
			gbar0 = gbar(th0,gk0);

			double ysm0 = star_y*gsm0;
			double yfr0 = star_y*gfr0;
			lambdasq = star_y*star_y - ysm0*ysm0 - yfr0*yfr0;
			lambda = sqrt(lambdasq);

			if ( ! (met==AF || met==AF2) ) {
				double star_z2 = star_y*gbar0;
				double star_z3 = (*pv1h*gfr0)*(*pxh*gsm0) - (*pv1h*gsm0)*(*pxh*gfr0);
				c1 = -lambda*star_z3;
				c2 = star_z + star_z2*star_z3;
			}
		}
	}


	if (model==M2) {

		if (th0ex) z_tilde_M2 = 0.; else  z_tilde_M2 = star_y*gfr(th0,gk0);
		lambdasq = star_y*star_y - z_tilde_M2*z_tilde_M2;
		lambda = sqrt(lambdasq);

	}

	c = sqrt( max( 1 - omega/lambdasq , 0. )  );


	return;
}
