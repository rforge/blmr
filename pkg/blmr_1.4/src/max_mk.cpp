//

#include "blmr.h"


double Cblmr::max_mk()  const
// get the maximum of (gamma(theta)*qy )^2 .  
// Algorithm in paper "Computing Confidence Intervals in Broken Line Regression".
{
	double  vf= (x_d[n_d-1]-x_d[n_d-2])*(*pqy*pq1[n_d-1]),  max= vf*vf/qff[n_d-1];

	for (int k=n_d-3;k>=k1;k--) {

		const double  v1 = (*pqy)*pq1[k+1],  vx = vf + v1*x_d[k+1];  
		vf = vx - v1*x_d[k];

		const double va = vx*qx1[k+1] - v1*qxx[k+1], vb = vx*q11[k+1] - v1*qx1[k+1], thk = va/vb;

		double mk;
		if(x_d[k]<thk && thk<x_d[k+1]) mk =(vx*vb - v1*va)/ck[k+1]; else  mk =vf*vf/qff[k+1];
			
		if (mk > max)  max = mk;
	}

	return max;
}



bool Cblmr::m_gt_w(const double wsq, const Vector<double> &s)  const
// Version of max_mk for Monte Carlo evaluation of SL.  Check whether  max(<gam(theta).u>^2)  >  wsq .
// Use pre-calculated  M*"stump of 1"  vectors for the dot products of 'u' 
// by gamma*u = gamma*(M-transpose*s) = (M*gamma)*s .  
{
	double  sf= (x_d[n_d-1]-x_d[n_d-2])*(s*pmq1[n_d-1]),  mk= sf*sf/qff[n_d-1];
	if (mk > wsq)  return true; 

	for (int k=n_d-3;k>=k1;k--) 
	{
		const double  s1 = s*pmq1[k+1],  sx = sf + s1*x_d[k+1];  
		sf = sx - s1*x_d[k];

		const double sa = sx*qx1[k+1] - s1*qxx[k+1], sb = sx*q11[k+1] - s1*qx1[k+1], thk = sa/sb;

		if(x_d[k]<thk && thk<x_d[k+1]) mk =(sx*sb - s1*sa)/ck[k+1];  else  mk =sf*sf/qff[k+1];

		if (mk > wsq)  return true;
	}

	return false;
}


