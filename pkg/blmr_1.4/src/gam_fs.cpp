//

#include "blmr.h"



Vector<double> Cblmr::gam(double th, int k) const
// calculate gamma
{
	if ( (Model==M1 && th <=xs[0]) || xs[ns-1] <= th) return *vund; else
		if (xs[ns-2] <= th  &&  th < xs[ns-1]) return *puqen; else
			if (Model==M1  &&  (xs[0] < th  &&  th <= xs[1]) ) return *puqe1; else
				if (Model==M2  &&  th <= xs[0]) return *puqx; else

	return 1./sqrt( ff(th,k) ) * q_f(th,k);
}


Vector<double> Cblmr::gfr(double th, int k) const
{
	if (th >= xs[ns-1]) return *vund2; else
		if (xs[ns-2] <= th  &&  th < xs[ns-1]) return *pusen; else

	return 1./sqrt(sf(th,k)*sf(th,k)) * sf(th,k);
}


Vector<double> Cblmr::gsm(double th, int k) const
{
	if (th <= xs[0]) return *vund2; else

	if (cov_matrix_diagonal) {
		if (xs[0] < th  &&  th <= xs[1]) return *pnuse1; else
			return  1./sqrt(sfc(th,k)*sfc(th,k)) * sfc(th,k);
	} else {
		if (xs[0] < th  &&  th <= xs[1]) {
			double e1fr = *pnse1*gfr(th,0);
			return 1./sqrt(se1sq-e1fr*e1fr) * (*pnse1 - e1fr*gfr(th,0));
		} else {
			Vector<double> f1(n);
			f1 = sfc(th,k) - (sfc(th,k)*sf(th,k))/(sf(th,k)*sf(th,k)) * sf(th,k);
			return  1./sqrt(f1*f1) * f1;
		}
	}
}


Vector<double> Cblmr::gbar(double th, int k) const
{
	if (th <= xs[0] ||  xs[ns-1] <= th ) return *vund2; else {

	Vector<double> fbar(n);
	fbar = *pv1h - (*pv1h*gsm(th,k))*gsm(th,k) - (*pv1h*gfr(th,k))*gfr(th,k);
	return  1./sqrt(fbar*fbar) * fbar;

	}
}



Vector<double> Cblmr::q_f(double th, int k) const
{
	return pqx[k] - th*pq1[k];
}



Vector<double> Cblmr::sf(double th, int k) const
{
	return psx[k] - th*ps1[k];
}



Vector<double> Cblmr::sfc(double th, int k) const
{
	return psxc[k] - th*ps1c[k];
}



double Cblmr::ff(double th, int k) const
{
	return qxx[k] + (q11[k]*th - 2*qx1[k])*th;
}

