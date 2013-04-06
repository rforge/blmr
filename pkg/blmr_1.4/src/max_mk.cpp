//

#include "blmr.h"


double Cblmr::max_mk() const
// get the maximum of (gamma(theta)*qy )^2 .  
// Algorithm in paper "Computing Confidence Intervals in Broken Line Regression".
{
	double  vf= (xs[ns-1]-xs[ns-2])*(*pqy*pq1[ns-2]), max= vf*vf/qff[ns-2];

	for (int k=ns-3;k>=k1;k--) {

		double  mk,  v1 = (*pqy)*pq1[k],  vx = vf + v1*xs[k+1];  vf = vx - v1*xs[k];

		double va = vx*qx1[k] - v1*qxx[k], vb = vx*q11[k] - v1*qx1[k], thk = va/vb;

		if(xs[k]<thk && thk<xs[k+1]) mk =(vx*vb - v1*va)/ck[k]; else  mk =vf*vf/qff[k];
			
		if (mk > max)  max = mk;
	}

	return max;
}



bool Cblmr::m_gt_w(double wsq, Vector<double> &s)
// Version of max_mk for Monte Carlo evaluation.  Check whether  <gam(theta).u>^2  >  wsq .
// Pre-calculate M*(1-tilde vector) for the dot products of u 
// by gamma*u = gamma*(M-transpose*s) = (M*gamma)*s .  
{
	double  sf= (xs[ns-1]-xs[ns-2])*(s*pmq1[ns-2]), mk= sf*sf/qff[ns-2];
	if (mk > wsq)  return true; 

	for (int k=ns-3;k>=k1;k--) 
	{
		double  s1 = s*pmq1[k],  sx = sf + s1*xs[k+1];  sf = sx - s1*xs[k];

		double sa = sx*qx1[k] - s1*qxx[k], sb = sx*q11[k] - s1*qx1[k], thk = sa/sb;

		if(xs[k]<thk && thk<xs[k+1]) mk =(sx*sb - s1*sa)/ck[k];  else  mk =sf*sf/qff[k];

		if (mk > wsq)  return true;
	}

	return false;
}


