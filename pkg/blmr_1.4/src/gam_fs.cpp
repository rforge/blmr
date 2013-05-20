//

#include "blmr.h"


 
Vector<double> Cblmr::gam(const double th, const int k) const
// calculate gamma
{
	if ( (Model==M1 && th <=x_d[0]) || x_d[n_d-1] <= th) return *und_m; else
		if (x_d[n_d-2] <= th  &&  th < x_d[n_d-1]) return *puqen; else
			if (Model==M1  &&  (x_d[0] < th  &&  th <= x_d[1]) ) return *puqe1; else
				if (Model==M2  &&  th <= x_d[0]) return *puqx; else

	return  1./sqrt( ff(th,k) ) * q_f(th,k);
}


 
Vector<double> Cblmr::gfr(const double th, const int k) const
// gamma frown
{
	if (th >= x_d[n_d-1]) return *und_n; else
		if (x_d[n_d-2] <= th  &&  th < x_d[n_d-1]) return *pusen; else

	return  1./sqrt(sf(th,k)*sf(th,k)) * sf(th,k);
}


 
Vector<double> Cblmr::gsm(const double th, const int k) const
// this function "gamma-smile" used in Model M1 only
{
	if (th <= x_d[0])  return *und_n;  else
		if (cov_matrix_diagonal) {
			if (x_d[0] < th  &&  th < x_d[1])  return *pnuse1;  else
				return  1./sqrt(sfc(th,k)*sfc(th,k)) * sfc(th,k);
		} else {
			if (x_d[0] < th  &&  th <= x_d[1]) {
				const double e1fr = *pnse1*gfr(th,k);
				return  1./sqrt(se1sq - e1fr*e1fr) * (*pnse1 - e1fr*gfr(th,k));
			} else {
				Vector<double> f1(n);
				f1 = sfc(th,k) - (sfc(th,k)*sf(th,k))/(sf(th,k)*sf(th,k)) * sf(th,k);
				return  1./sqrt(f1*f1) * f1;
			}
		}
}


 
Vector<double> Cblmr::gbar(const double th, const int k) const
// gamma bar
{
	if ( (Model==M1 && th <= x_d[0])  ||  x_d[n_d-1] <= th )  return *und_n;  else  {

		Vector<double> fbar(n);
		fbar = *pv1h - (*pv1h*gsm(th,k))*gsm(th,k) - (*pv1h*gfr(th,k))*gfr(th,k);

		return  1./sqrt(fbar*fbar) * fbar;
	}
}


 
Vector<double> Cblmr::gbar_prime(const double th, const int k) const
// gamma bar prime
{
	if ( (Model==M1 && th <= x_d[0])  ||  x_d[n_d-1] <= th )  return *und_n;  else  {

		Vector<double> fbar(n);
		fbar = *pv1h - (*pv1h*gfr(th,k))*gfr(th,k);

		return  1./sqrt(fbar*fbar) * fbar;
	}
}


 
Vector<double> Cblmr::q_f(const double th, const int k) const
{
	return  pqx[k] - th*pq1[k];
}


 
Vector<double> Cblmr::sf(const double th, const int k) const
{
	return  psx[k] - th*ps1[k];
}


 
Vector<double> Cblmr::sfc(const double th, const int k) const
{
	return  psxc[k] - th*ps1c[k];
}



double Cblmr::ff(const double th, const int k) const
{
	return  qxx[k] + (q11[k]*th - 2*qx1[k])*th;
}


