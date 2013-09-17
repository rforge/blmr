//

#include "lmbr.h"




void Clmbr::set_x( void )
// precalculate working variables and vectors 
// that depend only on 'x' values and 'Sigma' matrix
{
	int i,j;


	double  X[n][xrank];
	for (i=0;i<n;i++)  for (j=0;j<xrank;j++)  
		if( model_in > 0 )  X[i][j] = *(x_in+j*n+i);  else  X[i][j] = *(x_in+j*n+(n-1-i));

	int xcol;
	if(Model==M3) xcol=0; else xcol=1;
	if( model_in < 0 )  for (i=0;i<n;i++)  X[i][xcol] *= -1.;

	Vector<double>  x(n);
	for (i=0;i<n;i++) x[i] = X[i][xcol];
	const double  min_xdiff= (x[n-1]-x[0])*0.001;
	double  xi= x[0] - 1 - min_xdiff,  xib;

	for(i=0;i<n;i++) {
		xib= xi;
		xi= x[i];
		if ( isinf( xi ) || isnan( xi ) )  stop( _("invalid 'x' value") );
		if( xib > xi )  stop( _("'x' values must be non-decreasing") );
		const double xdiff= xi - xib;
		if ( 0 < xdiff  &&  xdiff < min_xdiff  ) {
			Rcout << _("consider a repeat predictor value instead of the values") << endl;
			Rcout << xib << ",  " << xi << endl;
			Rf_warning( _("predictor values might be too close for reliable computations") );
		}
	}


	ns= 0;
	for(i=1;i<n;i++)  if( x[i] != x[i-1] )  is[ns++]= i-1;
	is[ns++]= n-1;


	bool  lack= false;
	if(Model==M1)  if( ns < 4 )  lack= true;
	if(Model==M2)  if( ns < 3 )  lack= true;
	if(Model==M3)  if( ns < 2 )  lack= true;
	if(variance_unknown)  if( m < 3 )  lack= true;
	if(lack)  stop( _("number of seperate 'x' values below minimum for changepoint inference") );



//  Canonically reduce the model, using an orthogonal matrix Q for change-of-basis.
//  Use LAPACK routines DGEQRF to generate Q and DORMQR to multiply by Q .


//	compute  X =  irS * X :
	if( vectorS )  for(j=0;j<xrank;j++)  for(i=0;i<n;i++)  X[i][j] *=  *( irS + i );
	if( matrixS )  for(j=0;j<xrank;j++)  {
		double  temp[ n ];
		for(i=0;i<n;i++)  temp[i] = X[i][j];
		for(i=0;i<n;i++) {
			X[i][j] = 0.;
			for( int k=0; k<n; k++) X[i][j] += *( irS + k*n + i ) * temp[k];
		}
	} 


// 'Q' and 'tau' are global arrays, declared in "lmbr.h"
// set the starting columns of 'Q' to be the columns of irS*X ,  
// but with the  1- and x1- vectors as the final two cols
	for(i=0;i<n;i++) for(j=xcol+1;j<xrank;j++)  *(Q + (j-xcol-1)*n + i ) =  X[i][j];
	for(i=0;i<n;i++) for(j=0;j<xcol+1;j++)  *(Q + (j+xrank-xcol-1)*n + i ) =  X[i][j];


	int  lwork,  info;
	double  tmp[1];
	{
		lwork = -1;

		F77_CALL(dgeqrf)( &n, &xrank, Q, &n, tau, tmp, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );  else  lwork= *tmp; 
		double  work[lwork];

		F77_CALL(dgeqrf)( &n, &xrank, Q, &n, tau, work, &lwork, &info );

		if( info )  stop( _("LAPACK routine 'dgeqrf' failed") );
	}


// first, use DORMQR to setup some vectors and constants

	{
		Vector<double>  e1(n,0.),  en(n,0.);
		for(i=0;i<=is[0];i++)  e1[i] = 1.;
		for(i=is[ns-2]+1;i<=is[ns-1];i++)  en[i] = 1.;

		const int  nC0 = 4;
		double  C0[n*nC0];
		for(i=0;i<n;i++) {
			*(C0+i) = 1;
			*(C0+n+i) = x[i];
			*(C0+2*n+i) = e1[i];
			*(C0+3*n+i) = en[i];
		}

		if( vectorS )  for(j=0;j<nC0;j++)  for(i=0;i<n;i++)  *(C0+n*j+i) *=  *( irS + i );
		if( matrixS )  for(j=0;j<nC0;j++)  {
			double  temp[ n ];
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
			double  work[lwork]; 

			F77_CALL(dormqr)( &side, &tp, &n, &nC0, &xrank, Q, &n, tau, C0, &n, work, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );
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


		const Vector<double>  dummy_m1(m1,0.),  dummy_m(m,0.); 

		*px = x;
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
	}


// array sizes that depend on 'ns'

	if( xs != NULL )  {

		delete[] xs;
		delete[] ps1;  delete[] psx;  
		delete[] q11;  delete[] qx1;  delete[] qxx;  delete[] ck;  delete[] qff;
		delete[] q10;  delete[] qx0;  delete[] a0;  delete[] b0;
		delete[] f01;  delete[] f0x;  delete[] B;
		delete[] pq1;  delete[] pqx;  

		xs = NULL;
		ps1 = psx = NULL;
		q11 = qx1 = qxx = ck = qff = NULL;
		q10 = qx0 = a0 = b0 = NULL;
		f01 = f0x = NULL;  B = NULL;
		pq1 = pqx = NULL;
	}


	try{
		xs = new double[ns];
		ps1 = new (Vector<double>[(ns+1)*m1]);
		psx = new (Vector<double>[(ns+1)*m1]);
		pq1 = new (Vector<double>[(ns+1)*m]);
		pqx = new (Vector<double>[(ns+1)*m]);
		q11 = new double[ns+1];
		qx1 = new double[ns+1];
		qxx = new double[ns+1];
		ck = new double[ns+1];
		qff = new double[ns+1];
		q10 = new double[ns+1];
		qx0 = new double[ns+1];
		a0 = new double[ns+1];
		b0 = new double[ns+1];
		f01 = new double[ns+1];
		f0x = new double[ns+1];
		B = new double[ns+1];

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << ex.what() << endl;
		stop( _("memory allocation failed") );
	}



	for(i=0;i<ns;i++)  xs[i] = x[ is[i] ];


// pre-calculate Q*"stub 1-" and Q*"stub x-" vectors and scalars

	{
		const int  nC1 = ns+1;
		double  C1[n*nC1],  Cx[n*nC1];

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
			double  temp[ n ];
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
			double  work[lwork];

			F77_CALL(dormqr)( &side, &tp, &n, &nC1, &xrank, Q, &n, tau, C1, &n, work, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );
		}
	

		if( vectorS )  for(j=0;j<nC1;j++)  for(i=0;i<n;i++)  *(Cx+n*j+i) *=  *( irS + i );
		if( matrixS )  for(j=0;j<nC1;j++)  {
			double  temp[ n ];
			for(i=0;i<n;i++)  temp[i] = *(Cx+n*j+i);
			for(i=0;i<n;i++) {
				*(Cx+n*j+i) = 0.;
				for( int k=0; k<n; k++) *(Cx+n*j+i) += *( irS + k*n + i ) * temp[k];
			}
		} 


		{
			lwork= -1;

			F77_CALL(dormqr)( &side, &tp, &n, &nC1, &xrank, Q, &n, tau, Cx, &n, tmp, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );  else  lwork= *tmp; 
			double  work[lwork];

			F77_CALL(dormqr)( &side, &tp, &n, &nC1, &xrank, Q, &n, tau, Cx, &n, work, &lwork, &info );

			if( info )  stop( _("LAPACK routine 'dormqr' failed") );
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


