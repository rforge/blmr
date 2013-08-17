//

#include "lmbr.h"




void Clmbr::set_x( void )
// precalculate working variables that depend only on x values
// and Sigma matrix
{
	int i,j;


	Vector<double>  x(n);
	int xcol;
	if(Model==M3) xcol=0; else xcol=1;
	for (i=0;i<n;i++) {
		if( model_in > 0 )  x[i]= *( X_in + i*xrank + xcol );  else  x[i]= *( X_in + (n-1-i)*xrank + xcol )*(-1.); 
	}

	const double  min_xdiff= (x[n-1]-x[0])*0.001;
	double  xi= x[0] - 1 - min_xdiff,  xib;

	for(i=0;i<n;i++) {
		xib= xi;
		xi= x[i];
		if ( isinf( xi ) || isnan( xi ) )  stop( _("invalid 'x' value") );
		if( xib > xi )  stop( _("'x' values must be non-decreasing") );
		const double xdiff= xi - xib;
		if ( 0 < xdiff  &&  xdiff < min_xdiff  ) {
			Rcout << _("suggest combining these data into two observations at a repeat 'x' value:") << endl;
			Rcout << "x[" << i-1 << "] =" << xib << " ,  x[" << i << "] =" << xi << endl;
			Rf_warning( _("'x' values too close for reliable computations") );
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


	get_Q();
	Matrix<double>  QS(m,n),  Q1S(m1,n);
	QS = (*pQ)*(*pirS);
	Q1S = (*pQ1)*(*pirS);
	for (i=0;i<m;i++) for (j=0;j<n;j++) if (fabs(QS[i][j])<zero_eq) QS[i][j]=0.;
	for (i=0;i<m1;i++) for (j=0;j<n;j++) if (fabs(Q1S[i][j])<zero_eq) Q1S[i][j]=0.;


	const Vector<double>  v1(n,1.), dummy_m1(m1,0.), dummy_m(m,0.); 
	Vector<double>  s1(m1), sx(m1);

	*px = x;
	s1 = Q1S*v1;	for(i=0;i<m1;i++) if(fabs(s1[i])<zero_eq) s1[i]=0.;
	sx = Q1S*x;		for(i=0;i<m1;i++) if(fabs(sx[i])<zero_eq) sx[i]=0.;
	*psig1 = s1;
	*psigx = sx;

	s11 = s1*s1;
	sx1 = sx*s1;
	sxx = sx*sx;
	n1 = sqrt(s11);
	*pv1h= 1./n1 * s1;

	const double x1 = sx*(*pv1h);
	*pxh = 1./sqrt( sx*sx - x1*x1 ) * ( sx - x1*(*pv1h) );


	Vector<double>  e1(n,0.),  en(n,0.),  sen(m1);

	for(i=0;i<=is[0];i++)  e1[i] = 1.;
	for(i=is[ns-2]+1;i<=is[ns-1];i++)  en[i] = 1.;
	*nan_m1 = dummy_m1;
	(*nan_m1)[0] = NaN;
	*pnse1 = -1.*(Q1S*e1);		for(i=0;i<m1;i++) if(fabs((*pnse1)[i])<zero_eq) (*pnse1)[i]=0.;
	se1sq = *pnse1*(*pnse1);
	*pnuse1 = 1./sqrt(se1sq) * (*pnse1);
	sen = Q1S*en;				for(i=0;i<m1;i++) if(fabs(sen[i])<zero_eq) sen[i]=0.;
	*pusen = 1./sqrt(sen*sen) * sen;



	*nan_m = dummy_m;
	(*nan_m)[0] = NaN;
	*puqe1 = 1./sqrt( (QS*e1)*(QS*e1) ) * (QS*e1);
	*puqen = 1./sqrt( (QS*en)*(QS*en) ) * (QS*en);
	for(j=0;j<m;j++) if( fabs((*puqe1)[j]) < zero_eq )  (*puqe1)[j] = 0.;
	for(j=0;j<m;j++) if( fabs((*puqen)[j]) < zero_eq )  (*puqen)[j] = 0.;

	if(Model==M1) *puqx = *nan_m;  else {
		*puqx = 1./sqrt( (QS*x)*(QS*x) ) * (QS*x);
		for(j=0;j<m;j++) if( fabs((*puqx)[j]) < zero_eq )  (*puqx)[j] = 0.;
	}


// array sizes that depend on 'ns'

	if( xs != NULL ) {

		delete[] xs;
		delete[] ps1;  delete[] psx;  
		delete[] ps1c;  delete[] psxc;
		delete[] q11;  delete[] qx1;  delete[] qxx;  delete[] ck;  delete[] qff;
		delete[] q10;  delete[] qx0;  delete[] a0;  delete[] b0;
		delete[] f01;  delete[] f0x;  delete[] B;
		delete[] pmq1;  delete[] pq1;  delete[] pqx;  

		xs = NULL;
		ps1 = psx = ps1c = psxc = NULL;
		q11 = qx1 = qxx = ck = qff = NULL;
		q10 = qx0 = a0 = b0 = NULL;
		f01 = f0x = NULL;
		B = NULL;
		pmq1 = pq1 = pqx = NULL;
	}


	try{
		xs = new double[ns];
		ps1 = new (Vector<double>[(ns+1)*m1]);
		psx = new (Vector<double>[(ns+1)*m1]);
		ps1c = new (Vector<double>[(ns+1)*m1]);
		psxc = new (Vector<double>[(ns+1)*m1]);
		pmq1 = new (Vector<double>[(ns+1)*m]);
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
		Rcout << _("message: ") << 5 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}



	for(i=0;i<ns;i++)  xs[i] = x[ is[i] ];



	Vector<double>  vk1(n,1.), kx(n), q1(m), qx(m);

	kx = *px;


	for (i=0;i<ns+1;i++) {

		s1 = Q1S*vk1;
		sx = Q1S*kx;

		for(j=0;j<m1;j++) {
			if (fabs(s1[j]) < zero_eq) s1[j] = 0.;
			if (fabs(sx[j]) < zero_eq) sx[j] = 0.;
		}
		ps1[i] = s1;
		psx[i] = sx;

		s1 = Q1S*(v1-vk1);
		sx = Q1S*(x-kx);
		for(j=0;j<m1;j++) {
			if (fabs(s1[j]) < zero_eq) s1[j] = 0.;
			if (fabs(sx[j]) < zero_eq) sx[j] = 0.;
		}
		ps1c[i] = s1;
		psxc[i] = sx;


		q1 = QS*vk1;
		qx = QS*kx;

		for(j=0;j<m;j++) {
			if (fabs(q1[j]) < zero_eq) q1[j] = 0.;
			if (fabs(qx[j]) < zero_eq) qx[j] = 0.;
		}

		pq1[i] = q1;
		pqx[i] = qx;

		q11[i] = q1*q1;							if (fabs(q11[i]) < zero_eq) q11[i] = 0.;
		qx1[i] = qx*q1;							if (fabs(qx1[i]) < zero_eq) qx1[i] = 0.;
		qxx[i] = qx*qx;							if (fabs(qxx[i]) < zero_eq) qxx[i] = 0.;
		ck[i] = qxx[i]*q11[i] - qx1[i]*qx1[i];  if (fabs(ck[i]) < zero_eq) ck[i] = 0.;
		if(i==0)  qff[i]= NaN;  else  { qff[i]= (qx-xs[i-1]*q1)*(qx-xs[i-1]*q1);  
		if (fabs(qff[i]) < zero_eq) qff[i]= 0.; }


		if(i==0) {
			for(j=0;j<=is[i];j++) vk1[j] = 0.;
			for(j=0;j<=is[i];j++) kx[j] = 0.;
		} else  if(i<ns) {
			for(j=is[i-1]+1;j<=is[i];j++) vk1[j] = 0.;
			for(j=is[i-1]+1;j<=is[i];j++) kx[j] = 0.;
		}

	}

	
	Lgamma = 0.;
	double gg;
	if(k1== -1)  {
		gg = gam(-Inf,0)*gam(xs[0],0);
		gg = min( 1., gg);
		gg = max( -1., gg);
		Lgamma += acos( gg );
	}
	for (i=max(k1,0);i<ns-2;i++)  {
		gg = gam(xs[i],i)*gam(xs[i+1],i+1);
		gg = min( 1., gg);
		gg = max( -1., gg);
		Lgamma += acos( gg );
	}

	if ( py != NULL )   set_y();


	return;
}


