//

#include "blmr.h"


bool args2real(SEXP* args, int nargs) { 
	if( nargs != 2) return false;
	return ( TYPEOF( args[1] ) == REALSXP ); }

bool args3real(SEXP* args, int nargs) { 
	if( nargs != 3) return false;
	return ( TYPEOF( args[1] ) == REALSXP ); }

bool args4real(SEXP* args, int nargs) { 
	if( nargs != 4) return false;
	return ( TYPEOF( args[1] ) == REALSXP ); }



RCPP_MODULE(Cblmr){

	using namespace Rcpp ;


	void (Cblmr::*sl1) (double) = &Cblmr::slR ;
	void (Cblmr::*sl2) (double, double) = &Cblmr::slR ;
	void (Cblmr::*sl3) (double, string) = &Cblmr::slR ;
	void (Cblmr::*sl4) (double, double, string) = &Cblmr::slR ;
	void (Cblmr::*sl5) (double, string, double) = &Cblmr::slR ;
	void (Cblmr::*sl6) (double, double, string, double) = &Cblmr::slR ;
	double (Cblmr::*sl7) (double, string, double, bool) = &Cblmr::slR ;
	double (Cblmr::*sl8) (double, double, string, double, bool) = &Cblmr::slR ;

	void (Cblmr::*ci1) (double) = &Cblmr::ciR ;
	void (Cblmr::*ci2) (double, string) = &Cblmr::ciR ;
	void (Cblmr::*ci3) (void) = &Cblmr::ciR ;

	void (Cblmr::*cr1) (double) = &Cblmr::crR ;
	void (Cblmr::*cr2) (double, string) = &Cblmr::crR ;
	void (Cblmr::*cr3) (double, string, double) = &Cblmr::crR ;
	void (Cblmr::*cr4) (void) = &Cblmr::crR ;



	class_<Cblmr>( "Cblmr" )
	
	    // expose the constructor
	    .constructor< NumericVector, NumericVector, int, bool, NumericMatrix >()    

		// expose methods
		.method( "sl", sl1 , "significance level for theta0 by default method CLR-GEO" )
		.method( "sl", sl2 , "SL for (theta0, alpha0) by default method", &args2real )
		.method( "sl", sl3 , "SL for th0 by specified method" )
		.method( "sl", sl4 , "SL for (th0,a0) by specified method", &args3real )
		.method( "sl", sl5 , "SL for th0 by specified method and accuracy" )
		.method( "sl", sl6 , "SL for (th0,a0) by specified method and accuracy", &args4real )
		.method( "sl", sl7 , "SL for th0 by specified method, accuracy and output flag" )
		.method( "sl", sl8 , "SL for (th0,a0) by specified method, accuracy and output flag" )
		.method( "ci", ci1 , "printout confidence interval for theta" )
		.method( "ci", ci2 , "confidence interval for theta by specified method" )
		.method( "ci", ci3 , "confidence interval for theta at default confidence level and method" )
		.method( "cr", cr1 , "confidence region for (theta,alpha)" )
		.method( "cr", cr2 , "confidence region for (theta,alpha) by specified method" )
		.method( "cr", cr3 , "confidence region for (theta,alpha) by specified method and increment" )
		.method( "cr", cr4 , "confidence region for (theta,alpha) at default confidence level" )
		.method( "mle", &Cblmr::MLE , "maximum likelihood estimates of parameters" )
		.method( "sety", &Cblmr::SET_irSy , 
			"set new values for inverse-root of Sigma matrix times y-vector" )
	;

}                     


