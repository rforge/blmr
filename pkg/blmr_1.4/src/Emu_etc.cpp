//


#include "blmr.h"



double Cblmr::Emupr(double th, int k)
// calculate (E-mu+)*pr, in Knowles,Siegmund,Zhang's formula
{
	if (k==gk0 || fabs(th-th0)<zero_eq) return 0.;

	double ro = rho(th,k), drosq = drhosq(th,k), r=1-ro*ro, zz = 1.-z*z,
		wzr = w-z*ro, zwr = z-w*ro;

	if (r < zero_eq) return 0.;
	
	double musq = drosq*zwr*zwr/r/r;

	double pisq = zz - wzr*wzr/r;

	if ( pisq < zero_eq )  return  0.;

	double OMsq = dgsq(th,k) - drosq/r;

	double tausq =  OMsq * pisq ;

	if( musq/tausq >1-zero_eq )  return 0.;

	double amu=sqrt(musq), tau=sqrt(tausq), ambt=amu/tau;

	double rr= sqrt(r*zz), g= wzr/rr, pr= fk(m-2,g)/rr;

	return  pr*( -amu*F(m-3,-ambt) + tau*fk(m-1,ambt)/(m-2) );
}



double Cblmr::Emupr_vk(double th, int k)
// calculate (E-mu+)*pr, in geometric-expectation formula for variance known
{
	if (k==gk0 || fabs(th-th0)<zero_eq) return 0.;

	double ro = rho(th,k), dro = drho(th,k), r=1-ro*ro, rr = sqrt( r );

	if (r < zero_eq) return 0.;
	
	double namu = -fabs( dro*(z-w*ro)/r );

	double OMarg = dgsq(th,k) - dro*dro/r;
	if (OMarg < zero_eq) return 0.; 
	double OM = sqrt( OMarg );

	double  mbO = namu/OM;

	double pr = Rf_dnorm4( (w-z*ro)/rr ,0,1,0)/rr;

	double pn_mbO =  Rf_pnorm5(mbO ,0,1,1,0) , dn_mbO = Rf_dnorm4(mbO ,0,1,0) ;

	return  ( namu*pn_mbO + OM*dn_mbO )*pr;
}







