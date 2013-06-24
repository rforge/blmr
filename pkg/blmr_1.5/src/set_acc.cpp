//


#include "blmr.h"




void  Cblmr::set_acc( const double acc )
// set the accuracies, absolute and relative for significance level estimates,
// for x- and y- boundaries of confidence intervals and regions,
// and increments for grid searches and integration limits
{
	if ( isnan(acc) || acc<0 || acc>1 )  stop( _("invalid 'acc' value") );
	
	subints = 10;	// minimum number of subintervals per data interval for grid searches

	acc_rho = 0.0001;		//  tolerance for finding 'rho' values by 'bisect' or 'rho_inv' 

	acc_sl_abs = acc;					// maximum absolute error in significance level estimates
	acc_sl_rel = min( 10*acc, 0.01 );	// maximum relative error in significance level estimates


	int i;

// maximum error in x- boundaries
	acc_xb = (xs[ns-1] - xs[0])*acc_sl_rel/64.;	
	i=1;  while( acc_xb < ldexp(1.,-i) )  i++;  acc_xb = ldexp(1.,-i);

// maximum error in y- boundaries
	double maxY= -Inf, minY= Inf;
	for(i=0;i<n;i++)  {  if( (*py)[i] > maxY )  maxY= (*py)[i];  if( (*py)[i] < minY )  minY= (*py)[i];  }
	const double dY= maxY - minY;
	acc_yb = dY*acc_sl_rel/64.;	
	i=1;  while( acc_yb < ldexp(1.,-i) )  i++;  acc_yb = ldexp(1.,-i);

//accuracy for integration limits
	inc_x = acc_xb;
	double  xinc;
	for( i= max( k1, 0 );  i < ns-2;  i++ )  { 
		xinc= ( xs[i+1] - xs[i] )/subints;  
		if( xinc < inc_x )  inc_x= xinc; 
	}
	i=1;  while( inc_x < ldexp(1.,-i) )  i++;  inc_x = ldexp(1.,-i);

//starting increment for y- grid searches
	inc_y = dY/128;
	i=1;  while( inc_y < ldexp(1.,-i) )  i++;  inc_y = ldexp(1.,-i);

//check for trivial case
	trivial= false;
	if ( variance_unknown  &&  omega/m < zero_eq )  trivial= true;


	return;
}


