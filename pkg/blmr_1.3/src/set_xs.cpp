//

#include "blmr.h"




void Cblmr::set_xs(void)
// precalculate working variables that depend only on x values
// and Sigma matrix
{
	int i,j;
	const double zero_eq = ldexp(1.,-40);


	if(Xs != NULL) {
		delete[] Xs;
		delete[] pmq1;  delete[] pq1;  delete[] pqx;  delete[] ps1;  delete[] psx;  
		delete[] ps1c;  delete[] psxc;
		delete[] q11;  delete[] qx1;  delete[] qxx;  delete[] ck;  delete[] qff;
		delete[] q10;  delete[] qx0;  delete[] a0;  delete[] b0;
	}


	Xs = new double[ns];
	pmq1 = new (Vector<double>[ns*m]);
	pq1 = new (Vector<double>[ns*m]);
	pqx = new (Vector<double>[ns*m]);
	ps1 = new (Vector<double>[ns*n]);
	psx = new (Vector<double>[ns*n]);
	ps1c = new (Vector<double>[ns*n]);
	psxc = new (Vector<double>[ns*n]);
	q11 = new double[ns];
	qx1 = new double[ns];
	qxx = new double[ns];
	ck = new double[ns];
	qff = new double[ns];
	q10 = new double[ns];
	qx0 = new double[ns];
	a0 = new double[ns];
	b0 = new double[ns];


	for(i=0;i<ns;i++) Xs[i] = Xd[ is[i] ];


	Matrix<double> QS(m,n);
	QS = (*pQ)*(*pirS);
	for (i=0;i<m;i++) for (j=0;j<n;j++) if (fabs(QS[i][j])<zero_eq) QS[i][j]=0.;


	Vector<double> v1(n,1.), x(n), vk1(n,1.), kx(n), s1(n), sx(n), q1(m), qx(m);

	x = *px;
	kx = *px;


	for (i=0;i<ns;i++) {


		if(i==0) {
			for(j=0;j<=is[i];j++) vk1[j] = 0.;
			for(j=0;j<=is[i];j++) kx[j] = 0.;
		} else {
			for(j=is[i-1]+1;j<=is[i];j++) vk1[j] = 0.;
			for(j=is[i-1]+1;j<=is[i];j++) kx[j] = 0.;
		}

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
		qff[i] = (qx-Xs[i]*q1)*(qx-Xs[i]*q1);	if (fabs(qff[i]) < zero_eq) qff[i] = 0.;

	}

	Lgamma = 0.;
	for (i=k1;i<ns-2;i++) Lgamma += acos( gam(Xs[i],i)*gam(Xs[i+1],i+1) );



	if (Yd != NULL) set_y(Yd);

	return;
}


