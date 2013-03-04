//


#include "blmr.h"



double Cblmr::Emupr(double th, int k)
// calculate (E-mu+)*pr, in Knowles,Siegmund,Zhang's formula
// arranged to minimize sqrt's
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

	double mbt_invsq = tausq/musq;

	if ( mbt_invsq < 1.+zero_eq) return 0.;

// arrange to minimize sqrt's
	double  g_invsq = r*zz/wzr/wzr,  gsq = wzr*wzr/r/zz;

	double  mbd = musq/r/zz,  tbd = tausq/r/zz;

	double f1= fk_arg_invsq(m-2,g_invsq), F1= F_arg_neginvsq(m-3,mbt_invsq), fF= f1*F1;
	double f2= fk_arg_invsq(m-2,g_invsq), f22= fk_arg_invsq(m-1,mbt_invsq), ffd= f2*f22/(m-2.);

	double result;
	if(m%2) result= ffd*sqrt(tbd*(1.-gsq)) - fF*sqrt(mbd*(1.- gsq)); 
		else  result= ffd*sqrt(tbd-mbd) - fF*sqrt(mbd);

	return result;
}



double Cblmr::Emupr_vk(double th, int k)
// calculate (E-mu+)*pr, in Knowles,Siegmund,Zhang's formula for variance known
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



double Cblmr::fk_arg_invsq(int k, double arg)
// calculate the f() function defined in K,S&Z eq.(13)
// arg inverse squared, and if k odd take k-1, i.e. leave out the sqrt
{
	if ( arg < 1.+ zero_eq) return 0.;
	int pk = k;
	if ( k%2 ) pk = k-1;
	double term = 1.-1./arg, prod = term;
	for (int i=1;i<pk/2-1;i++) prod *= term;
	return  C[n-1-k]*prod;
}



double Cblmr::F_arg_neginvsq(int k, double arg)
// calculate the F() function defined in K,S&Z eq.(14)
{
	if (arg <= 1.) return 0.;
	return 1.- Rf_pt( sqrt(k*1./(arg-1.)) , k , 1, 0);
}






