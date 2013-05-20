//


#include "blmr.h"



double Cblmr::mu_by_tau(const double th, const int k)  const
{
	const double  ro = rho(th,k), dro = drho(th,k), r=1-ro*ro, wzr = w-z*ro;
	
	double  mu = dro*(z-w*ro)/r;
	if(fabs(mu)<zero_eq) mu= 0.;

	if (fabs(th-th0) < zero_eq || k==k0)  return  copysign( inf, mu );

	const double  pisq = 1. - z*z - wzr*wzr/r;
	if ( pisq < zero_eq )  return  copysign( inf, mu );

	const double  OMsq = dgsq(th,k) - dro*dro/r;
	if ( OMsq < zero_eq )  return  copysign( inf, mu );

	const double  tau = sqrt( OMsq * pisq ) ;

	return mu/tau;
}




double Cblmr::amu_by_Omega(const double th, const int k)  const
{
	if (k==k0 || fabs(th-th0)<zero_eq)  return inf;

	const double  ro = rho(th,k), drosq = drhosq(th,k), r=1-ro*ro, zwr = z-w*ro;

	const double  musq = drosq*zwr*zwr/r/r;

	const double  OMsq = dgsq(th,k) - drosq/r;

	return  sqrt(musq/OMsq);
}




double Cblmr::Emupr(const double th, const int k)  const
// calculate (E-mu+)*pr, in Knowles,Siegmund,Zhang's formula
{
	if (k==k0 || fabs(th-th0)<zero_eq) return 0.;

	const double  ro = rho(th,k), drosq = drhosq(th,k), r=1-ro*ro, zz = 1.-z*z,
		wzr = w-z*ro, zwr = z-w*ro;

	if (r < zero_eq) return 0.;
	
	double  musq = drosq*zwr*zwr/r/r;
	if(fabs(musq)<zero_eq) musq= 0.;

	const double  pisq = zz - wzr*wzr/r;
	if ( pisq < zero_eq )  return  0.;

	const double  OMsq = dgsq(th,k) - drosq/r;
	if( OMsq < zero_eq )  return 0.;

	const double tausq =  OMsq * pisq ;

	if( musq/tausq >1-zero_eq )  return 0.;

	const double  tau=sqrt(tausq), ambt= sqrt(musq/tausq);

	const double  rr= sqrt(r*zz), g= wzr/rr, pr= fk(m-2,g)/rr;

	return  pr*( tau*sF(m-3,-ambt) );
}



double Cblmr::Emupr_vk(const double th, const int k)  const
// calculate (E-mu+)*pr, in geometric-expectation formula for variance known
{
	if (k==k0 || fabs(th-th0)<zero_eq)  return 0.;

	const double  ro = rho(th,k), dro = drho(th,k), r=1-ro*ro, rr = sqrt( r );

	if (r < zero_eq) return 0.;
	
	const double  namu = -fabs( dro*(z-w*ro)/r );

	const double  OMsq = dgsq(th,k) - dro*dro/r;
	if (OMsq < zero_eq) return 0.; 
	const double  OM = sqrt( OMsq );

	const double  mbO = namu/OM;

	const double  pr = Rf_dnorm4( (w-z*ro)/rr ,0,1,0)/rr;

	const double  pn_mbO =  Rf_pnorm5(mbO ,0,1,1,0) ,  dn_mbO = Rf_dnorm4(mbO ,0,1,0) ;

	return  ( namu*pn_mbO + OM*dn_mbO )*pr;
}







