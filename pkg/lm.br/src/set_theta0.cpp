//

#include "lmbr.h"



void  Clmbr::set_theta0( const double th_0, const METHOD met)
// precalculate numbers and vectors that depend on theta0, including 'w' and 'z'
{
	if ( !isinf(th_0) && isnan(th_0) )  stop( _("invalid 'theta0' value") );

	int k_0 = 0;
	while ( xs[k_0] - th_0 < zero_eq  &&  k_0 < ns ) k_0++;
	th0 = th_0;
	if( fabs(th_0) < zero_eq )  th0= 0.;
	if( k_0 > 0 )   if ( fabs( xs[k_0-1] - th_0 ) < zero_eq) th0 = xs[k_0-1];
	if( k_0 < ns )   if ( fabs( xs[k_0] - th_0 ) < zero_eq)  { th0 = xs[k_0];  k_0++; }
	k0 = k_0;


//  'th0ex' = true  if 'th0' is exterior to  'x'  values 

	if ( (Model==M1 && th0<=xs[0]) || xs[ns-1]<=th0 )  th0ex = true;  else  th0ex = false;



	double  wsq,  max_gy;
	mle( false, &max_gy );
	if (variance_unknown)  wsq = max_gy/qysq;  else  wsq = max_gy;
	w = sqrt( max( 0., wsq ) );


	if (th0ex) {

		prime_z = 0.;
		z = 0.;

	} else {

		Vector<double>  g0(m), qf0(m);
		g0 = gam(th0,k0);
		qf0 = q_f(th0,k0);
		for(int i=0; i<m; i++)  if( fabs(g0[i]) < zero_eq )  g0[i] = 0.;
		for(int i=0; i<m; i++)  if( fabs(qf0[i]) < zero_eq )  qf0[i] = 0.;

		prime_z = *pqy*g0;

		if (variance_unknown)  z = prime_z/sqrt(qysq);  else  z = prime_z;


		if (met==GEO || met==GEO2 || met==INIT)  {

			g0u2 = *puqen*g0;   if( fabs(g0u2) < zero_eq )  g0u2= 0.;
			if(Model==M1)  g0u1 = *puqe1*g0;  else  
				if(Model==M2)  g0u1 = *puqx*g0;  else
					if(Model==M3)  g0u1 = *pv1h*g0; 
			if( fabs(g0u1) < zero_eq )  g0u1= 0.;

			for (int k=0;k<ns+1;k++)  {
				q10[k] = pq1[k]*g0;  if (fabs(q10[k]) < zero_eq) q10[k] = 0.;
				qx0[k] = pqx[k]*g0;  if (fabs(qx0[k]) < zero_eq) qx0[k] = 0.;
				a0[k] = qx1[k]*qx0[k] - qxx[k]*q10[k];  if (fabs(a0[k]) < zero_eq) a0[k] = 0.;
				b0[k] = q11[k]*qx0[k] - qx1[k]*q10[k];  if (fabs(b0[k]) < zero_eq) b0[k] = 0.;
				f01[k] = qf0*pq1[k];  if (fabs(f01[k]) < zero_eq) f01[k] = 0.;
				f0x[k] = qf0*pqx[k];  if (fabs(f0x[k]) < zero_eq) f0x[k] = 0.;

// calculate B[k]
				if( ck[k]==0. )  B[k]= 1.;  else  {

					if( ( xs[ns-2] <= th0 && th0 < xs[ns-1] )  ||
							( Model==M1 && xs[0]<th0 && th0<=xs[1] )  ||  
								( Model==M2 && th0<=xs[0] ) || (Model==M3 && isinf(th0))  ) 	// 'th0' on an end-interval

						B[k] =  ( qxx[k]*q10[k]*q10[k] - 2*qx1[k]*q10[k]*qx0[k] + q11[k]*qx0[k]*qx0[k] ) / ck[k] ;

					else						
						B[k] =  ( qxx[k]*f01[k]*f01[k] - 2*qx1[k]*f01[k]*f0x[k] + q11[k]*f0x[k]*f0x[k] ) / ck[k] / ff(th0,k0);

			  		if (fabs(B[k]) < zero_eq)  B[k] = 0.;
					if ( B[k] > 1 - zero_eq )  B[k] = 1.;
				}
			}
		}
	}

	if( fabs( w - fabs(z) ) < zero_eq )  w = fabs(z);


	if (met==MC || met==INIT) {	

		Matrix<double>   M(m,m,0.);

		if (th0ex)  { for(int i=0;i<m;i++) M[i][i] =1.; }  else  get_M( &M );

		for (int k=0; k < ns + 1; k++)  {
			pmq1[k] = M*pq1[k];  
			for(int i=0;i<m;i++)  if( fabs(pmq1[k][i]) < zero_eq )  pmq1[k][i] = 0.;
		}

		if(Model==M3) {
			*pm1h= M*(*pv1h);  
			for(int i=0;i<m;i++)  if( fabs( (*pm1h)[i] ) < zero_eq )  (*pm1h)[i] = 0.;
		}
	}


	return;
}

