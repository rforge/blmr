//

#include "lmbr.h"



void  Clmbr::set_theta0( const double th_0, const METHOD met )
// precalculate numbers and vectors that depend on theta0, including 'w' and 'z'
{
	if ( !isinf(th_0) && isnan(th_0) )  stop( _("invalid 'theta0' value") );

	if( th_0 != th0 )  { //1

		int k_0 = 0;
		while ( xs[k_0] - th_0 < zero_eq  &&  k_0 < ns )  k_0++;
		th0 = th_0;
		if( fabs(th_0) < zero_eq )  th0= 0.;
		if( k_0 > 0 )   if ( fabs( xs[k_0-1] - th_0 ) < zero_eq )  th0 = xs[k_0-1];
		if( k_0 < ns )   if ( fabs( xs[k_0] - th_0 ) < zero_eq )  { th0 = xs[k_0];  k_0++; }
		k0 = k_0;


//  'th0ex' = true  if 'th0' is exterior to  'x'  values 

		if ( (Model==M1 && th0<=xs[0]) || xs[ns-1]<=th0 )  th0ex = true;  else  th0ex = false;


		double  wsq,  max_gy;
		mle( false, &max_gy );
		if (variance_unknown)  wsq = max_gy/qysq;  else  wsq = max_gy;
		w = sqrt( max( 0., wsq ) );


		if (th0ex) { //2

			prime_z = 0.;
			z = 0.;

		} else { //2

			Vector<double>  g0(m), qf0(m);
			g0 = gam(th0,k0);
			qf0 = q_f(th0,k0);

			prime_z = *pqy*g0;

			if (variance_unknown)  z = prime_z/sqrt(qysq);  else  z = prime_z;


			if (met==GEO || met==GEO2 || met==INIT)  { //3

				g0u2 = *puqen*g0;  
				if(Model==M1)  g0u1 = *puqe1*g0;  else  
					if(Model==M2)  g0u1 = *puqx*g0;  else
						if(Model==M3)  g0u1 = *pv1h*g0; 

				for (int k=0;k<ns+1;k++)  { //4
					q10[k] = pq1[k]*g0;  
					qx0[k] = pqx[k]*g0;  
					a0[k] = qx1[k]*qx0[k] - qxx[k]*q10[k];  
					b0[k] = q11[k]*qx0[k] - qx1[k]*q10[k];  
					f01[k] = qf0*pq1[k];  
					f0x[k] = qf0*pqx[k];  

// calculate B[k]
					if( ck[k] == 0 )  B[k]= 1.;  else  { //5

						if( ( xs[ns-2] <= th0 && th0 < xs[ns-1] )  ||
								( Model==M1 && xs[0]<th0 && th0<=xs[1] )  ||  
									( Model==M2 && th0<=xs[0] ) || (Model==M3 && isinf(th0))  ) 	// 'th0' on an end-interval

							B[k] =  ( qxx[k]*q10[k]*q10[k] - 2*qx1[k]*q10[k]*qx0[k] + q11[k]*qx0[k]*qx0[k] ) / ck[k] ;

						else						
							B[k] =  ( qxx[k]*f01[k]*f01[k] - 2*qx1[k]*f01[k]*f0x[k] + q11[k]*f0x[k]*f0x[k] ) / ck[k] / ff(th0,k0);

				  		if ( B[k] < 0. )  B[k] = 0.;
						if ( B[k] > 1. )  B[k] = 1.;
					} //5
				} //4
			} //3
		} //2

		if( fabs( w - fabs(z) ) < zero_eq )  w = fabs(z);

	} //1


	if ( met==MC  &&  !th0ex  &&  th_0 != th0MC )  {	

// to speed-up Monte Carlo evaluation, pre-multiply vectors by  
// orthogonal matrix 'M' which has first row = 'gamma(th0)' .
// Use LAPACK routines DGEQRF to generate 'M' and DORMQR to multiply by 'M' .
 
		int  i,  j,  ng =1;
		Vector<double>  g0(m);
		g0 = gam(th0,k0);
		double  M[m],  _tau[ng];
		for(i=0;i<n;i++)  M[i] = g0[i];

		int  lwork = -1,  info;
		double  tmp[1];
		{
			F77_CALL(dgeqrf)( &m, &ng, M, &m, _tau, tmp, &lwork, &info );
			if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );  else  lwork= *tmp; 
			double  work[lwork];
			F77_CALL(dgeqrf)( &m, &ng, M, &m, _tau, work, &lwork, &info );
			if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );
		}

		const int  nCm = ns+2;
		double  Cm[m*nCm];
		Vector<double>  cq(m);
		for(j=0;j<nCm-1;j++)  { cq = pq1[j];  for(i=0;i<m;i++)  *(Cm+j*m+i) = cq[i]; }
		cq = *pv1h;  for(i=0;i<m;i++)  *(Cm+(nCm-1)*m+i) = cq[i]; 

		{
			const char  side = 'L',  tp = 'N';
			lwork= -1;
			F77_CALL(dormqr)( &side, &tp, &m, &nCm, &ng, M, &m, _tau, Cm, &m, tmp, &lwork, &info );
			if( info )  stop( _("LAPACK routine 'dormqr' failed") );  else  lwork= *tmp; 
			double  work[lwork];
			F77_CALL(dormqr)( &side, &tp, &m, &nCm, &ng, M, &m, _tau, Cm, &m, work, &lwork, &info );
			if( info )  stop( _("LAPACK routine 'dormqr' failed") );
		}

//			if( pmq1 != NULL )  {  delete[] pmq1;  delete[] pm1h;  pmq1 = pm1h = NULL;  }

		if( pmq1 == NULL )  {
			try{
				pmq1 = new (Vector<double>[(ns+1)*m]);
				if(Model==M3)  pm1h = new (Vector<double>[m]);

			} catch( bad_alloc &ex ) {
				Rcout << _("message: ") << ex.what() << endl;
				stop( _("memory allocation failed") );
			}
		}

		for(j=0;j<nCm-1;j++) { for(i=0;i<m;i++) cq[i]= *(Cm+j*m+i);  pmq1[j]= cq; } 
		if(Model==M3)  { for(i=0;i<m;i++) cq[i]= *(Cm+(nCm-1)*m+i);  *pm1h= cq; }

		th0MC = th0;
	}


	return;
}

