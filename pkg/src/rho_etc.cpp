//

#include "blmr.h"





double Cblmr::mu_by_tau(double th, int k)
{
	double ro = rho(th,k), dro = drho(th,k), r=1-ro*ro, wzr = w-z*ro;
	
	double mu = dro*(z-w*ro)/r;

	if (fabs(th-th0) < zero_eq || k==gk0)  return  copysign( inf, mu );

	double pisq = 1. - z*z - wzr*wzr/r;
	if ( pisq < zero_eq )  return  copysign( inf, mu );

	double OMsq = dgsq(th,k) - dro*dro/r;

	if ( OMsq < zero_eq )  return  copysign( inf, mu );

	double tau = sqrt( OMsq * pisq ) ;

	return mu/tau;
}




double Cblmr::mu_by_tau_sq(double th, int k)
{
	if (fabs(th-th0) < zero_eq) return inf;

	double ro = rho(th,k), drosq = drhosq(th,k), r=1-ro*ro, zz = 1.-z*z,
		wzr = w-z*ro, zwr = z-w*ro;

	double  sign_dro = a0[k] - b0[k]*th;
	if (th < th0)  sign_dro = -sign_dro;

	double sign_mu = sign_dro*zwr;
	
	if (Xs[k]<th0+zero_eq && th0<Xs[k+1]+zero_eq) return  copysign( inf, sign_mu );

	double musq = drosq*zwr*zwr/r/r;

	double pisq = zz - wzr*wzr/r;

	if ( pisq < zero_eq )  return  copysign( inf, sign_mu );

	double OMsq = dgsq(th,k) - drosq/r;

	if ( OMsq < zero_eq )  return  copysign( inf, sign_mu );

	double tausq =  OMsq * pisq ;

	return  copysign( musq/tausq, sign_mu ) ;
}



double Cblmr::rho_inv(double s, int k)
// returns th such that rho(th) = s  and th is in (x[k], x[k+1])
// rhok_inv gives the same result for s and for -s 
// if both possible values of th are in the data interval then 
// returns the lower value if s<0, the greater value if s>=0
{
	double s2 = s*s;
	double a = s2*q11[k]-q10[k]*q10[k], b = 2*(q10[k]*qx0[k]-s2*qx1[k]), c = s2*qxx[k]-qx0[k]*qx0[k];
	double r = b*b-4*a*c;
	if (r < 0.) {
		Rcout << "no inverse-rho in interval " << k << "  for rho= " << s << ",  rad^2= " << r << endl;
		return und;
	}
	double rad = sqrt(r);
	double r1 = (-b-rad)/2/a,  r2 = (-b+rad)/2/a;

	bool i1 = false, i2 = false;
	if (Xs[k] <= r1  &&  r1 <= Xs[k+1]) i1 = true;
	if (Xs[k] <= r2  &&  r2 <= Xs[k+1]) i2 = true;
	if (i1 && i2) { if (s < 0) return min(r1,r2); else return max(r1,r2); }
	if (i1) return r1;
	if (i2) return r2;
	Rcout << "no inverse-rho in interval " << k << " for rho= " << s << endl;
	return und;
}





double Cblmr::rho(double th)
{
	int k=0;
	while (Xs[k+1]<th) k++;
	return rho(th,k);
}


double Cblmr::rho(double th, int k)
{
	return (qx0[k] - q10[k]*th)/sqrt(ff(th,k));
}



double Cblmr::rhosq(double th, int k)
{
	double num = qx0[k] - q10[k]*th;
	return num*num/ff(th,k);
}



double Cblmr::drho(double th, int k)
{
	double  fsq = ff(th,k), dro = (a0[k] - b0[k]*th)/sqrt(fsq)/fsq;
	if (th < th0)  dro = -dro;
	return dro;
}


double Cblmr::drhosq(double th, int k)
{
	double  ab = a0[k] - b0[k]*th, fsq = ff(th,k);
	return ab*ab/(fsq*fsq*fsq);
}



double Cblmr::dgsq(double th, int k)
{
	double  fsq = ff(th,k);
	return ck[k]/fsq/fsq;
}

