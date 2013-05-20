//

#include "blmr.h"




void Cblmr::set_x( void )
// precalculate working variables that depend only on x values
// and Sigma matrix
{
	int i,j;


	const double  min_xdiff= (X_in[n-1]-X_in[0])*0.001;
	double  xi= X_in[0] - 1 - min_xdiff,  xib;

	for(i=0;i<n;i++) {
		xib= xi;
		xi= X_in[i];
		if ( isinf( xi ) || isnan( xi ) )  stop( _("invalid 'x' value") );
		if( xib > xi )  stop( _("'x' values must be non-decreasing") );
		const double xdiff= xi - xib;
		if ( 0 < xdiff  &&  xdiff < min_xdiff  ) {
			Rcout << _("suggest combining these data into two observations at a repeat 'x' value:") << endl;
			Rcout << "x[" << i-1 << "] =" << xib << " ,  x[" << i << "] =" << xi << endl;
			Rf_warning( _("'x' values too close for reliable computations") );
		}
	}



	Vector<double> x(n);

	for (i=0;i<n;i++) {
		if( model_in != -2 )  x[i]= X_in[i];  else  x[i]= -X_in[n-1-i];
	}

	n_d= 0;
	for(i=1;i<n;i++)  if( x[i] != x[i-1] )  i_d[n_d++]= i-1;
	i_d[n_d++]= n-1;
	for(i=n_d;i<n;i++) i_d[i]= 0.;	//fill out 


	bool lack= false;
	if(Model==M1)  {
		if( n_d < 4 )  lack= true;
		if(variance_unknown)  if( n < 5 )  lack= true;
	}
	if(Model==M2)  {
		if( n_d < 3 )  lack= true;
		if(variance_unknown)  if( n < 4 )  lack= true;
	}
	if(lack)  stop( _("number of 'x' values less than minimum") );



	const Vector<double>  v1(n,1.), dummy_n(n,0.), dummy_m(m,0.); 
	Vector<double>  s1(n), sx(n);

	*px = x;
	s1 = *pirS*v1;
	sx = *pirS*x;
	*psig1 = s1;
	*psigx = sx;

	s11 = s1*s1;
	sx1 = sx*s1;
	sxx = sx*sx;
	n1 = sqrt(s11);
	*pv1h= 1./n1 * s1;

	const double x1 = sx*(*pv1h);
	*pxh = 1./sqrt( sx*sx - x1*x1 ) * ( sx - x1*(*pv1h) );


	Vector<double>  e1(n,0.),  en(n,0.),  sen(n);

	for(i=0;i<=i_d[0];i++)  e1[i] = 1.;
	for(i=i_d[n_d-2]+1;i<=i_d[n_d-1];i++)  en[i] = 1.;
	*und_n = dummy_n;
	(*und_n)[0] = und;
	*pnse1 = -1.*(*pirS*e1);
	se1sq = *pnse1*(*pnse1);
	*pnuse1 = 1./sqrt(se1sq) * (*pnse1);
	sen = *pirS*en;
	*pusen = 1./sqrt(sen*sen) * sen;


	get_Q();
	Matrix<double> QS(m,n);
	QS = (*pQ)*(*pirS);
	for (i=0;i<m;i++) for (j=0;j<n;j++) if (fabs(QS[i][j])<zero_eq) QS[i][j]=0.;

	*und_m = dummy_m;
	(*und_m)[0] = und;
	*puqe1 = 1./sqrt( (QS*e1)*(QS*e1) ) * (QS*e1);
	*puqen = 1./sqrt( (QS*en)*(QS*en) ) * (QS*en);

	if(Model==M1) *puqx = *und_m;  else {
		*puqx = 1./sqrt( (QS*x)*(QS*x) ) * (QS*x);
		for(j=0;j<m;j++) if( fabs((*puqx)[j]) < zero_eq )  (*puqx)[j] = 0.;
	}


// array sizes that depend on 'n_d'

	if( x_d != NULL ) {

		delete[] x_d;
		delete[] ps1;  delete[] psx;  
		delete[] ps1c;  delete[] psxc;
		delete[] q11;  delete[] qx1;  delete[] qxx;  delete[] ck;  delete[] qff;
		delete[] q10;  delete[] qx0;  delete[] a0;  delete[] b0;
		delete[] pmq1;  delete[] pq1;  delete[] pqx;  

		x_d = NULL;
		ps1 = psx = ps1c = psxc = NULL;
		q11 = qx1 = qxx = ck = qff = NULL;
		q10 = qx0 = a0 = b0 = NULL;
		pmq1 = pq1 = pqx = NULL;
	}


	try{
		x_d = new double[n_d];
		ps1 = new (Vector<double>[(n_d+1)*n]);
		psx = new (Vector<double>[(n_d+1)*n]);
		ps1c = new (Vector<double>[(n_d+1)*n]);
		psxc = new (Vector<double>[(n_d+1)*n]);
		pmq1 = new (Vector<double>[(n_d+1)*m]);
		pq1 = new (Vector<double>[(n_d+1)*m]);
		pqx = new (Vector<double>[(n_d+1)*m]);
		q11 = new double[n_d+1];
		qx1 = new double[n_d+1];
		qxx = new double[n_d+1];
		ck = new double[n_d+1];
		qff = new double[n_d+1];
		q10 = new double[n_d+1];
		qx0 = new double[n_d+1];
		a0 = new double[n_d+1];
		b0 = new double[n_d+1];

	} catch( bad_alloc &ex ) {
		Rcout << _("message: ") << 5 << " " << ex.what() << endl;
		stop( _("memory allocation failed") );
	}



	for(i=0;i<n_d;i++)  x_d[i] = x[ i_d[i] ];



	Vector<double>  vk1(n,1.), kx(n), q1(m), qx(m);

	kx = *px;


	for (i=0;i<n_d+1;i++) {

		s1 = *pirS*vk1;
		sx = *pirS*kx;

		for(j=0;j<n;j++) {
			if (fabs(s1[j]) < zero_eq) s1[j] = 0.;
			if (fabs(sx[j]) < zero_eq) sx[j] = 0.;
		}
		ps1[i] = s1;
		psx[i] = sx;

		s1 = *pirS*(v1-vk1);
		sx = *pirS*(x-kx);
		for(j=0;j<n;j++) {
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
		if(i==0)  qff[i]= und;  else  { qff[i]= (qx-x_d[i-1]*q1)*(qx-x_d[i-1]*q1);  if (fabs(qff[i]) < zero_eq) qff[i]= 0.; }


		if(i==0) {
			for(j=0;j<=i_d[i];j++) vk1[j] = 0.;
			for(j=0;j<=i_d[i];j++) kx[j] = 0.;
		} else  if(i<n_d) {
			for(j=i_d[i-1]+1;j<=i_d[i];j++) vk1[j] = 0.;
			for(j=i_d[i-1]+1;j<=i_d[i];j++) kx[j] = 0.;
		}

	}
	

	Lgamma = 0.;
	for (i=k1;i<n_d-2;i++)  Lgamma += acos( gam(x_d[i],i)*gam(x_d[i+1],i+1) );


	if ( py != NULL )   set_y();


	return;
}


