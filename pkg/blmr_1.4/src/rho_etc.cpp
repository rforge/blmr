//

#include "blmr.h"



double Cblmr::rho(const double th)  const
{
	int k=0;
	while (x_d[k]<th && k<n_d) k++;
	return rho(th,k);
}



double Cblmr::rho(const double th, const int k)  const
{
	if(th>=x_d[n_d-1])  return und;  else
		if(Model==M1  &&  th<x_d[0])  return und;  else  {
			const double fsq= ff(th,k);
			return (qx0[k] - q10[k]*th)/sqrt(fsq);
		}
}



double Cblmr::rhosq(const double th, const int k)  const
{
	if(th>=x_d[n_d-1])  return und;  else
		if(Model==M1  &&  th<x_d[0])  return und;  else  {
			const double fsq= ff(th,k);
			const double num= qx0[k] - q10[k]*th;
			return num*num/fsq;
		} 
}



double Cblmr::drho(const double th, const int k)  const
{
	if(th>=x_d[n_d-1])  return und;  else
		if(Model==M1  &&  th<x_d[0])  return und;  else  {
			const double fsq= ff(th,k);
			double dro = (a0[k] - b0[k]*th)/sqrt(fsq)/fsq; 
			if (th < th0)  dro = -dro;
			return dro;
		}
}



double Cblmr::drhosq(const double th, const int k)  const
{
	if(th>=x_d[n_d-1])  return und;  else
		if(Model==M1  &&  th<x_d[0])  return und;  else  {
			const double  fsq = ff(th,k);
			const double ab= a0[k] - b0[k]*th; 
			return ab*ab/(fsq*fsq*fsq);
		}
}



double Cblmr::dgsq(const double th, const int k)  const
{
	if(th>=x_d[n_d-1])  return und;  else
		if(Model==M1  &&  th<x_d[0])  return und;  else  {
			const double  fsq = ff(th,k);
			return ck[k]/fsq/fsq;
		}
}



double Cblmr::rho_inv(const double s, const int k, const int hi_lo)  const
// Returns 'th' such that rho(th) = s  and 'th' is in (x[k-1], x[k]);
// 'rho_inv' gives the same result for s and for -s, so it checks the sign of rho(th); 
// result 'th' is a quadratic root, and if both roots are in the data interval then 
// returns the lower value if 'hi_lo' < 0, the greater value if 'hi_lo' > 0 .
{
	if( k >= n_d-1 )  return und;     // 'rho' is constant on final data interval, in all models
	if( fabs( rho(x_d[k-1],k) - s ) < zero_eq )  return x_d[k-1];
	if(k<n_d-1) if( fabs( rho(x_d[k],k) - s ) < zero_eq )  return x_d[k];
	const double s2 = s*s;
	const double a = s2*q11[k]-q10[k]*q10[k], b = q10[k]*qx0[k]-s2*qx1[k], 
		c = s2*qxx[k]-qx0[k]*qx0[k];
	const double r = b*b-a*c;
	if( r<0  &&  fabs(s) < 1.e-4) {    //computation unreliable for small 's'
		const double r1 = -b/a;   
		if (x_d[k-1] <= r1  &&  r1 <= x_d[k])  return r1;
	}
	if (r >= 0.) {
		const double rad = sqrt(r);
		const double r1 = (-b-rad)/a,  r2 = (-b+rad)/a;

		bool i1 = false, i2 = false;
		if (x_d[k-1] <= r1  &&  r1 <= x_d[k]  &&  s*rho(r1,k) >=0 ) i1 = true;
		if (x_d[k-1] <= r2  &&  r2 <= x_d[k]  &&  s*rho(r2,k) >=0 ) i2 = true;
		if (i1 && i2) { if (hi_lo < 0) return min(r1,r2); else return max(r1,r2); }
		if (i1) return r1;
		if (i2) return r2;
	}
	Rcout << _("'rho_inv' routine called for  interval, rho, radicand = ") << endl;
	Rcout << k << " " << s << " " << r << endl;
	stop ( _("no inverse-rho for given 'rho' value in given data interval") );
	return und;
}


