
//   integrands in  Knowles, Siegmund, Zhang's  geometric-expectation formulae for conditional
//   likelihood-ratio significance levels 


#include "blmr.h"





double Cblmr::amu_by_Omega(const double th, const int k)  const
{
	if (k==k0 || fabs(th-th0)<zero_eq)  return Inf;

	const double  ro = rho(th,k),  rosq = rhosq(th,k),  r = 1 - rosq;
	double  zwr = fabs(z-w*ro);
	if( zwr < zero_eq )  return 0.;

	if( isinf(th) )  {
		if( B[k]-rosq < 0 )  return 0.;
		return  zwr * sqrt( (B[k]-rosq) / (1-B[k]) / r );

	}  else  {

		const double  drosq = drhosq(th,k),  OMsq = dgsq(th,k) - drosq/r;

		if( OMsq <= 0 )  return Inf;

		return  sqrt(drosq/OMsq)*zwr/r;
	}
}



double Cblmr::Emupr(const double th, const int k)  const
// calculate (E-mu+)*pr, in Knowles,Siegmund,Zhang's geometric-expectation formula for variance unknown
// as a function of theta
{
	if ( k==k0 || fabs(th-th0) < zero_eq || isinf(th) )  return 0.;

	const double  rosq = rhosq(th,k),  ro = rho(th,k),  r=1-rosq,  zz = 1.-z*z,  wzr = w-z*ro;
	double  drosq = drhosq(th,k),  zwr = fabs(z-w*ro);
	if ( drosq < zero_eq*zero_eq )  drosq= 0.;
	if ( zwr < zero_eq)  zwr= 0.;
	
	const double  pisq = zz - wzr*wzr/r;
	if ( pisq < zero_eq*zero_eq )  return  0.;

	const double  OMsq = dgsq(th,k) - drosq/r;
	if( OMsq < zero_eq*zero_eq )  return 0.;

	const double  tau =  sqrt( OMsq * pisq ),  ambt = sqrt(drosq) * zwr / r / tau ;

	if( ambt > 1 - zero_eq )  return 0.;

	const double  rr= sqrt(r*zz), g= wzr/rr, pr= fk(m-2,g)/rr;

	return  pr * tau * sF(m-3,-ambt) ;
}



double  Cblmr::Emupr_vk(const double th, const int k)  const
// calculate (E-mu+)*pr, the integrand in geometric-expectation formula for variance known
// as a function of theta
{
	if (k==k0 || fabs(th-th0)<zero_eq || isinf(th) )  return 0.;

	const double  rosq = rhosq(th,k),  r=1-rosq,  rr = sqrt(r),  ro= rho(th,k);
	double  zwr= fabs(z-w*ro);
	if( zwr < zero_eq )  zwr= 0.;

	const double  drosq = drhosq(th,k),  namu = -zwr*sqrt(drosq)/r ;

	const double  OMsq = dgsq(th,k) - drosq/r;
	if (OMsq < zero_eq*zero_eq)  return 0.; 

	const double  OM = sqrt( OMsq ),  mbO = -zwr*sqrt(drosq/OMsq)/r;

	const double  pr = Rf_dnorm4( (w-z*ro)/rr ,0,1,0)/rr;

	const double  pn_mbO =  Rf_pnorm5(mbO ,0,1,1,0) ,   dn_mbO = Rf_dnorm4(mbO ,0,1,0) ;

	return  ( namu*pn_mbO + OM*dn_mbO )*pr;
}







