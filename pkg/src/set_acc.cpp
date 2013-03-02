//


#include "blmr.h"




void  Cblmr::set_acc( double acc )
// set the accuracies, absolute and relative for significance level estimates,
// for x- and y- boundaries of confidence intervals and regions,
// and increments for grid searches and integration limits
{
	if ( isnan(acc) || acc<0 || acc>0.1 ) {
		Rcout << "Invalid acc value: " << acc << endl;
		Rcout << "New value not assigned. acc remains  " << acc_sl_abs << endl;
		return;
	}

	subints = 10;	// minimum number of subintervals per data interval 
					// for grid searches and int. limits

	acc_sl_abs = acc;		// maximum absolute error in significance level estimates
	acc_sl_rel = min(10*acc,0.01);	// maximum relative error in sl estimates

	int i, j;

// maximum error in x- boundaries
	acc_xb = (Xd[n-1] - Xd[0])*acc_sl_rel/64.;	
	j=1; while(acc_xb < ldexp(1.,-j))j++; acc_xb = ldexp(1.,-j);

// maximum error in y- boundaries
	double maxY= -inf, minY= inf;
	for(i=0;i<n;i++) { if(Yd[i]>maxY) maxY=Yd[i]; if(Yd[i]<minY) minY=Yd[i]; }
	acc_yb = (maxY - minY)*acc_sl_rel/64.;	
	j=1; while(acc_yb < ldexp(1.,-j))j++; acc_yb = ldexp(1.,-j);

//accuracy for integration limits
	inc_x = (Xd[n-1] - Xd[0])/128;
	double xinc;
	for(i=k1;i<ns-2;i++) { xinc=(Xs[i+1]-Xs[i])/subints; if(xinc<inc_x) inc_x= xinc; }
	j=1; while(inc_x < ldexp(1.,-j))j++; inc_x = ldexp(1.,-j);

//starting increment for y- grid searches
	inc_y = (maxY - minY)/128;
	j=1; while(inc_y < ldexp(1.,-j))j++; inc_y = ldexp(1.,-j);


	aex = sqrt( Rf_qchisq(1-acc/2, m ,1,0) );	//used in  geo_vk_ex
	

	return;
}

