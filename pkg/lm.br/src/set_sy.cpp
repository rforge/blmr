//

#include "lmbr.h"



void  Clmbr::set_sy(double *const irsy, const METHOD met)
// recalculate numbers and vectors based on  rW*y = irS*y  values 
{
	int i;
	for (i=0;i<n;i++)
		if ( isinf(irsy[i]) || isnan(irsy[i]) )  stop( _("invalid value in 'rWy' vector") );

	double  virsy[n];
	for (i=0;i<n;i++)  if( model_in > 0 )  virsy[i] = irsy[i];  else  virsy[i] = irsy[n-1-i];

	Vector<double>  vy(n),  sy(m1),  qy(m);

	for(i=0;i<n;i++)  vy[i] =  virsy[i];
	if( vectorS )  for(i=0;i<n;i++)  vy[i] =  *( rS + i ) * virsy[i];
	if( matrixS )  for(i=0;i<n;i++)  {
		vy[i] = 0.;
		for( int k=0; k<n; k++)  vy[i] += *( rS + k*n + i ) * virsy[k];
	} 
	*py = vy;


// multiply by Q
	{
		const char  side = 'L',  tp = 'T';
		int  ny =1,  lwork= -1, info;
		double tmp[1];

		F77_CALL(dormqr)( &side, &tp, &n, &ny, &xrank, Q, &n, tau, virsy, &n, tmp, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dormqr' failed") );  else  lwork= *tmp;
		double  work[lwork]; 

		F77_CALL(dormqr)( &side, &tp, &n, &ny, &xrank, Q, &n, tau, virsy, &n, work, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dormqr' failed") );
	}


	for(i=0;i<m1;i++) sy[i]= *(virsy+i+n-m1);
	for(i=0;i<m;i++)  qy[i]= *(virsy+i+n-m);

	*psy = sy;
	*pqy = qy;

	y1 = *psy*(*psig1);
	yx = *psy*(*psigx);
	sysq = *psy*(*psy);  
	qysq = *pqy*(*pqy);

	double  max_gysq;
	mle( false, &max_gysq );
	omega =  qysq - max_gysq ;
	if(omega < 0)  omega= 0.;

	if ( met != INIT ) {
		const double  th0i = th0,  a0i = alpha0;
		th0 += 1;
		alpha0 += 1;
		set_theta0( th0i, met );
		set_alpha0( a0i, met );
		set_SL();
		set_acc();
	}

	return;
}

