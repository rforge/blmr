//

#include "blmr.h"



void  Cblmr::set_theta0(const double th_0, const METHOD met)
// precalculate numbers and vectors that depend on theta0
{
	if ( isinf(th_0) || isnan(th_0) )  stop( _("invalid 'theta0' value") );


	int k_0 = 0;
	while ( x_d[k_0] - th_0 < zero_eq  &&  k_0 < n_d ) k_0++;
	th0 = th_0;
	if( fabs(th_0) < zero_eq )  th0= 0.;
	if( k_0 > 0 )   if ( fabs( x_d[k_0-1] - th_0 ) < zero_eq) th0 = x_d[k_0-1];
	if( k_0 < n_d )   if ( fabs( x_d[k_0] - th_0 ) < zero_eq)  { th0 = x_d[k_0];  k_0++; }
	k0 = k_0;


//  'th0ex' = true  if 'th0' exterior to  'x'  values 
	if ( (Model==M1 && th0<=x_d[0]) || x_d[n_d-1]<=th0 )  th0ex = true;  else  th0ex = false;

	int i,k;

	if (th0ex) {

		prime_z = 0.;
		z = 0.;

	} else {

		Vector<double>  g0(m);
		g0 = gam(th0,k0);

		prime_z = *pqy*g0;

		if (met==GEO || met==GEO2 || met==INIT)  {
			for (k=0;k<n_d+1;k++)  {
				q10[k] = pq1[k]*g0;  if (fabs(q10[k]) < zero_eq) q10[k] = 0.;
				qx0[k] = pqx[k]*g0;  if (fabs(qx0[k]) < zero_eq) qx0[k] = 0.;
				a0[k] = qx1[k]*qx0[k] - qxx[k]*q10[k];  if (fabs(a0[k]) < zero_eq) a0[k] = 0.;
				b0[k] = q11[k]*qx0[k] - qx1[k]*q10[k];  if (fabs(b0[k]) < zero_eq) b0[k] = 0.;
			}
		}
	}


	if (met==MC || met==INIT) {	

		Matrix<double>   M(m,m,0.),  *pM = NULL;
		try{ pM = new (Matrix<double>[m*m]); }  catch( bad_alloc &ex ) {
			Rcout << _("message: ") << 9 << " " << ex.what() << endl;
			stop( _("memory allocation failed") );
		}

		*pM = M;
		if (th0ex)  { for(i=0;i<m;i++) (*pM)[i][i] =1.; }  else  get_M(pM);

		for (k=0;k<n_d+1;k++)  {
			pmq1[k] = (*pM)*pq1[k];  
			for(i=0;i<m;i++)  if( fabs(pmq1[k][i]) < zero_eq )  pmq1[k][i] = 0.;
		}

		delete[] pM;
	}


	return;
}

