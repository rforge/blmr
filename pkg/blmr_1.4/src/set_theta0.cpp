//

#include "blmr.h"



void  Cblmr::set_theta0(double th_0, METHOD met)
// precalculate numbers and vectors based on theta0, used in routines
{
	if ( isinf(th_0) || isnan(th_0) )  stop( _("invalid theta0 value") );


	int k0=0;
	while ( xs[k0+1] - th_0 < zero_eq  &&  k0 < ns-1 ) k0++;
	if ( fabs( xs[k0] - th_0 ) < zero_eq) th_0 = xs[k0];
	gk0 = k0;
	th0 = th_0;


// if th0 exterior to data set   th0ex = true  
	if ( (Model==M1 && th0<=xs[0]) || xs[ns-1]<=th0 ) th0ex = true; else th0ex = false;

	int i,k;

	if (th0ex) {

		star_z = 0.;
		z = 0.;

	} else {

		Vector<double> g0(m);
		g0 = gam(th0,k0);

		star_z = *pqy*g0;
		if(variance_unknown)  z = star_z/sqrt(qysq); else z = star_z;

		if (met==GEO || met==GEO2 || met==INIT)  for (k=0;k<ns;k++)  {
			q10[k] = pq1[k]*g0;  if (fabs(q10[k]) < zero_eq) q10[k] = 0.;
			qx0[k] = pqx[k]*g0;  if (fabs(qx0[k]) < zero_eq) qx0[k] = 0.;
			a0[k] = qx1[k]*qx0[k] - qxx[k]*q10[k];  if (fabs(a0[k]) < zero_eq) a0[k] = 0.;
			b0[k] = q11[k]*qx0[k] - qx1[k]*q10[k];  if (fabs(b0[k]) < zero_eq) b0[k] = 0.;
		}

	}


	if (met==MC || met==INIT) {

		Matrix<double> M(m,m,0.), *pM = NULL;
		pM = new (Matrix<double>[m*m]);
		*pM = M;
		if (th0ex) { for(i=0;i<m;i++) (*pM)[i][i] =1.; } else  get_M(pM);

		for (k=0;k<ns;k++)  {
			pmq1[k] = (*pM)*pq1[k];  
			for(i=0;i<m;i++)  if( fabs(pmq1[k][i]) < zero_eq )  pmq1[k][i] = 0.;
		}

		delete[] pM;
		pM = NULL;
	}


	return;
}

