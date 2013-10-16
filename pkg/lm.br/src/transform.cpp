//

#include "lmbr.h"



void  Clmbr::transform( void )
// transform to a canoncial model,  precalculate working variables and vectors 
// that depend only on 'x' values and 'Sigma' matrix
{
	int i,j;

	j= 0;
	for(i=1;i<n;i++)  if( (*px)[i] != (*px)[i-1] )  is[j++]= i-1;
	is[j]= n-1;

	for(i=0;i<ns;i++)  xs[i] = (*px)[ is[i] ];


	double *  X= Calloc( n*xrank, double );
	for (i=0;i<n;i++)  for (j=0;j<xrank;j++)  
		if( model_in > 0 )  *(X+j*n+i) = *(x_in+j*n+i);  else  *(X+j*n+i) = *(x_in+j*n+(n-1-i));

	int xcol;
	if(Model==M3) xcol=0; else xcol=1;
	if( model_in < 0 )  for (i=0;i<n;i++)  *(X+xcol*n+i) *= -1.;


//  Canonically reduce the model, using an orthogonal matrix Q for change-of-basis.
//  Use LAPACK routines DGEQRF to generate Q and DORMQR to multiply by Q .


//	compute  X =  irS * X :
	double *  temp= Calloc( n, double );
	if( vectorS )  for(j=0;j<xrank;j++)  for(i=0;i<n;i++)  *(X+j*n+i) *=  *( irS + i );
	if( matrixS )  for(j=0;j<xrank;j++)  {
		for(i=0;i<n;i++)  temp[i] = *(X+j*n+i);
		for(i=0;i<n;i++) {
			*(X+j*n+i) = 0.;
			for( int k=0; k<n; k++) *(X+j*n+i) += *( irS + k*n + i ) * temp[k];
		}
	} 


// 'Q' and 'tau' are global arrays, declared in "lmbr.h"
// set the pre-transform columns of 'Q' to be the columns of irS*X ,  
// but with the  1- and x1- vectors as the final two cols
	for(i=0;i<n;i++)  for(j=xcol+1;j<xrank;j++)  *(Q + (j-xcol-1)*n + i ) =  *(X+j*n+i);
	for(i=0;i<n;i++)  for(j=0;j<xcol+1;j++)  *(Q + (j+xrank-xcol-1)*n + i ) =  *(X+j*n+i);
	Free( X );

	int  lwork,  info;
	double  tmp[1];
	{
		lwork = -1;

		F77_CALL(dgeqrf)( &n, &xrank, Q, &n, tau, tmp, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );  else  lwork= *tmp; 
		double *  work= Calloc( lwork, double );

		F77_CALL(dgeqrf)( &n, &xrank, Q, &n, tau, work, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );
		Free( work );
	}


// first, use DORMQR to setup some vectors and constants

	{
		Vector<double>  e1(n,0.),  en(n,0.);
		for(i=0;i<=is[0];i++)  e1[i] = 1.;
		for(i=is[ns-2]+1;i<=is[ns-1];i++)  en[i] = 1.;

		const int  nC0 = 4;
		double *  C0= Calloc( n*nC0, double );
		for(i=0;i<n;i++) {
			*(C0+i) = 1;
			*(C0+n+i) = (*px)[i];
			*(C0+2*n+i) = e1[i];
			*(C0+3*n+i) = en[i];
		}

		if( vectorS )  for(j=0;j<nC0;j++)  for(i=0;i<n;i++)  *(C0+n*j+i) *=  *( irS + i );
		if( matrixS )  for(j=0;j<nC0;j++)  {
			for(i=0;i<n;i++)  temp[i] = *(C0+n*j+i);
			for(i=0;i<n;i++) {
				*(C0+n*j+i) = 0.;
				for( int k=0; k<n; k++) *(C0+n*j+i) += *( irS + k*n + i ) * temp[k];
			}
		} 


		const char  side = 'L',  tp = 'T';
		{
			lwork= -1;

			F77_CALL(dormqr)( &side, &tp, &n, &nC0, &xrank, Q, &n, tau, C0, &n, tmp, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );  else  lwork= *tmp;
			double *  work= Calloc( lwork, double );

			F77_CALL(dormqr)( &side, &tp, &n, &nC0, &xrank, Q, &n, tau, C0, &n, work, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );
			Free( work );
		}
 

		Vector<double>  s1(m1), sx(m1), se1(m1), sen(m1), qsx(m), qe1(m), qen(m);

		for(i=0;i<m1;i++) {
			const int  im1 = i+n-m1;
			s1[i] = *(C0+im1);
			sx[i] = *(C0+1*n+im1);
			se1[i]= *(C0+2*n+im1);
			sen[i]= *(C0+3*n+im1);
		}

		for(i=0;i<m;i++) {
			const int  im = i+n-m;
			qsx[i]= *(C0+1*n+im);
			qe1[i]= *(C0+2*n+im);
			qen[i]= *(C0+3*n+im);
		}
		Free( C0 );


		const Vector<double>  dummy_m1(m1,0.),  dummy_m(m,0.); 

		*psig1 = s1;
		*psigx = sx;

		s11 = s1*s1;
		sx1 = sx*s1;
		sxx = sx*sx;
		n1 = sqrt(s11);
		*pv1h = 1./n1 * s1;

		const double  x1 = sx*(*pv1h);
		*pxh = 1./sqrt( sx*sx - x1*x1 ) * ( sx - x1*(*pv1h) );


		*nan_m1 = dummy_m1;
		(*nan_m1)[0] = NaN;
		*pnse1 = -1.*se1;		
		se1sq = *pnse1*(*pnse1);
		*pnuse1 = 1./sqrt(se1sq) * (*pnse1);
		*pusen = 1./sqrt(sen*sen) * sen;

		*nan_m = dummy_m;
		(*nan_m)[0] = NaN;
		*puqe1 = 1./sqrt( qe1*qe1 ) * qe1;
		*puqen = 1./sqrt( qen*qen ) * qen;

		if(Model==M1) *puqx = *nan_m;  else  *puqx = 1./sqrt( qsx*qsx ) * qsx;

// initialize vectors for method="MC"
		for(j=0;j<ns+1;j++)  pmq1[j]= dummy_m; 
		if(Model==M3)  *pm1h= dummy_m; 
	}



// pre-calculate  Q* "stub 1-"  and  Q* "stub x-"  vectors and scalars

	{
		const int  nC1 = ns+1;
		double *  C1= Calloc( n*nC1, double );
		double *  Cx= Calloc( n*nC1, double );

		Vector<double>  vk1(n,1.), kx(n);
		kx = *px;

		for(j=0;j<nC1;j++) {

			for(i=0;i<n;i++) *(C1+j*n+i) = vk1[i];
			for(i=0;i<n;i++) *(Cx+j*n+i) = kx[i];

			if(j==0) {
				for(i=0;i<=is[j];i++) vk1[i] = 0.;
				for(i=0;i<=is[j];i++) kx[i] = 0.;
			} else  if(j<ns) {
				for(i=is[j-1]+1;i<=is[j];i++) vk1[i] = 0.;
				for(i=is[j-1]+1;i<=is[j];i++) kx[i] = 0.;
			}
		}

		if( vectorS )  for(j=0;j<nC1;j++)  for(i=0;i<n;i++)  *(C1+n*j+i) *=  *( irS + i );
		if( matrixS )  for(j=0;j<nC1;j++)  {
			for(i=0;i<n;i++)  temp[i] = *(C1+n*j+i);
			for(i=0;i<n;i++) {
				*(C1+n*j+i) = 0.;
				for( int k=0; k<n; k++) *(C1+n*j+i) += *( irS + k*n + i ) * temp[k];
			}
		} 

		const char  side = 'L',  tp = 'T';
		{
			lwork= -1;

			F77_CALL(dormqr)( &side, &tp, &n, &nC1, &xrank, Q, &n, tau, C1, &n, tmp, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );  else  lwork= *tmp; 
			double *  work= Calloc( lwork, double );

			F77_CALL(dormqr)( &side, &tp, &n, &nC1, &xrank, Q, &n, tau, C1, &n, work, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );
			Free( work );
		}
	

		if( vectorS )  for(j=0;j<nC1;j++)  for(i=0;i<n;i++)  *(Cx+n*j+i) *=  *( irS + i );
		if( matrixS )  for(j=0;j<nC1;j++)  {
			for(i=0;i<n;i++)  temp[i] = *(Cx+n*j+i);
			for(i=0;i<n;i++) {
				*(Cx+n*j+i) = 0.;
				for( int k=0; k<n; k++) *(Cx+n*j+i) += *( irS + k*n + i ) * temp[k];
			}
		} 
		Free( temp );


		{
			lwork= -1;

			F77_CALL(dormqr)( &side, &tp, &n, &nC1, &xrank, Q, &n, tau, Cx, &n, tmp, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );  else  lwork= *tmp; 
			double *  work= Calloc( lwork, double );

			F77_CALL(dormqr)( &side, &tp, &n, &nC1, &xrank, Q, &n, tau, Cx, &n, work, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );
			Free( work );
		}
	

		Vector<double>  s1(m1), sx(m1), q1(m), qx(m);


		for (j=0;j<ns+1;j++) {

			for(i=0;i<m1;i++) s1[i]= *(C1+j*n+i+n-m1);
			for(i=0;i<m1;i++) sx[i]= *(Cx+j*n+i+n-m1);

			ps1[j] = s1;
			psx[j] = sx;

			for(i=0;i<m;i++) q1[i]= *(C1+j*n+i+n-m);
			for(i=0;i<m;i++) qx[i]= *(Cx+j*n+i+n-m);

			pq1[j] = q1;
			pqx[j] = qx;

			q11[j] = q1*q1;							
			qx1[j] = qx*q1;							
			qxx[j] = qx*qx;							
			ck[j] = qxx[j]*q11[j] - qx1[j]*qx1[j];  
			if(j==0)  qff[j]= NaN;  else  qff[j]= (qx-xs[j-1]*q1)*(qx-xs[j-1]*q1);  
		}
		
		Free( C1 );  Free( Cx );
	}
	
	Lgamma = 0.;
	double gg;
	if(k1== -1)  {
		gg = gam(-Inf,0)*gam(xs[0],0);
		gg = min( 1., gg);  gg = max( -1., gg);
		Lgamma += acos( gg );
	}
	for (i=max(k1,0);i<ns-2;i++)  {
		gg = gam(xs[i],i)*gam(xs[i+1],i+1);
		gg = min( 1., gg);  gg = max( -1., gg);
		Lgamma += acos( gg );
	}


	if ( py != NULL )   set_y();

	return;
}

