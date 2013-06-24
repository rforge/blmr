//

#include "blmr.h"




const int  bis_it_limit= 50;	// maximum number of iterations in bisection routines



double Cblmr::bisect(const double a, const double b, double (Cblmr::*fn)(double,int), const int k, const double value, const double crit)
// find  x  such that  value < fn(x) < value + crit   if  crit > 0 ,   or   value - crit < fn(x) < value   if   crit < 0
{
	double  x1= a, x2= b, f1 = (this->*fn)(x1,k) - value,  f2 = (this->*fn)(x2,k) - value;
//Rcout << "bisect: x1 x2 f1 f2  " << x2 << " " << x2 << " " << f1 << " " << f2 << endl;
	if (f1*f2>0 || f1==f2 || (isnan(f1*f2) && !isinf(f1*f2)) ) {
		Rcout << _("'bisect' to find x such that f(x)=0 starting from input values:") << endl;
		Rcout << "x1, x2, f(x1), f(x2) =  " << x1 << " " << x2 << " " << f1 << " " << f2 << endl;
		stop( _("'bisect' cannot find interim point from starting values") );
	}
	int iteration=0;
	while (  fabs(x1-x2) > fabs(crit)  && (iteration<bis_it_limit) ) {
		const double  xmean = (x1+x2)/2,  fx = (this->*fn)(xmean,k)-value;
		if(f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; } 
//Rcout << "x1 x2 f1 f2  " << x2 << " " << x2 << " " << f1 << " " << f2 << endl;
		iteration++;
	}
	if(iteration==bis_it_limit)  Rf_warning( _("'bisect' failed to reach tolerance after maximum number of iterations") );
if(iteration==bis_it_limit)  stop("debug");
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}



double Cblmr::bisect(const double a, const double b, double (Cblmr::*fn)(double,int) const, const int k, const double value, const double crit)  const
// "const" version of above routine
{
	double  x1= a, x2= b, f1 = (this->*fn)(x1,k) - value,  f2 = (this->*fn)(x2,k) - value;
//Rcout << "bisect const: x1 x2 f1 f2  " << x2 << " " << x2 << " " << f1 << " " << f2 << endl;
	if (f1*f2>0 || f1==f2 || (isnan(f1*f2) && !isinf(f1*f2)) ) {
		Rcout << "const " << _("'bisect' to find x such that f(x)=0 starting from input values:") << endl;
		Rcout << "x1, x2, f(x1), f(x2) =  " << x1 << " " << x2 << " " << f1 << " " << f2 << endl;
		Rcout << "y " << *py << endl;
		stop( _("'bisect' const  cannot find interim point from starting values") );
	}
	int iteration=0;
	while (  fabs(x1-x2) > fabs(crit)  && (iteration<bis_it_limit) ) {
		const double  xmean = (x1+x2)/2,  fx = (this->*fn)(xmean,k)-value;
		if(f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; } 
//Rcout << "x1 x2 f1 f2  " << x2 << " " << x2 << " " << f1 << " " << f2 << endl;
		iteration++;
	}
	if(iteration==bis_it_limit)  Rf_warning( _("'bisect' const  failed to reach tolerance after maximum number of iterations") );
if(iteration==bis_it_limit)  Rcout << "y " << *py << endl;
if(iteration==bis_it_limit)  stop("debug");
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}




double Cblmr::bisect_sl(const double a, const double b, const METHOD met, double crit)
// use bisection to find x such that   SL < sl(x) < SL + crit   if  crit > 0 ,
// or    SL - crit < sl(x) < SL   if  crit < 0 ,   used in 'ci' routine
{
	double  x1= a, x2= b, f1= sl(x1,met,false)-SL, f2= sl(x2,met,false)-SL;
	if( fabs(f1)<zero_eq && fabs(f1-f2)<zero_eq) return (x1+x2)/2.;
	const double p12 = f1*f2;
	if (p12>0 || f1==f2 || fabs(p12)>1 || isnan(p12)) {
		Rcout << _("'bisect_sl' to find x such that sl(x) = critical SL, starting from input values:") << endl;
		Rcout << "x1, x2, sl(x1), sl(x2) =  " << x1 << " " << x2 << " " << f1+SL << " " << f2+SL << endl;
		stop( _("'bisect_sl' cannot find interim point from starting values") );
	}
	int iteration=0;
	while (  fabs(x1-x2) > fabs(crit)  && (iteration<bis_it_limit) ) {
		const double  xmean= (x1+x2)/2,  fx= sl(xmean,met,false) -SL;
		if (f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; }
		iteration++;
	}
	if(iteration==bis_it_limit)  Rf_warning( _("'bisect_sl' failed to reach tolerance after maximum number of iterations") );
if(iteration==bis_it_limit) stop("debug");
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}


