//

#include "lmbr.h"


bool args2real(SEXP* args, int nargs) { 
	if( nargs != 2) return false;
	return ( TYPEOF( args[1] ) == REALSXP ); }

bool args3real(SEXP* args, int nargs) { 
	if( nargs != 3) return false;
	return ( TYPEOF( args[1] ) == REALSXP ); }

bool args4real(SEXP* args, int nargs) { 
	if( nargs != 4) return false;
	return ( TYPEOF( args[1] ) == REALSXP ); }



RCPP_MODULE(Clmbr){

	using namespace Rcpp ;


	void (Clmbr::*sl1) (double) = &Clmbr::slR ;
	void (Clmbr::*sl2) (double, double) = &Clmbr::slR ;
	void (Clmbr::*sl3) (double, string) = &Clmbr::slR ;
	void (Clmbr::*sl4) (double, double, string) = &Clmbr::slR ;
	void (Clmbr::*sl5) (double, string, double) = &Clmbr::slR ;
	void (Clmbr::*sl6) (double, double, string, double) = &Clmbr::slR ;
	double (Clmbr::*sl7) (double, string, double, bool) = &Clmbr::slR ;
	double (Clmbr::*sl8) (double, double, string, double, bool) = &Clmbr::slR ;

	void (Clmbr::*ci1) (double) = &Clmbr::ciR ;
	void (Clmbr::*ci2) (double, string) = &Clmbr::ciR ;
	void (Clmbr::*ci3) (void) = &Clmbr::ciR ;

	void (Clmbr::*cr1) (double) = &Clmbr::crR ;
	void (Clmbr::*cr2) (double, string) = &Clmbr::crR ;
	void (Clmbr::*cr3) (double, string, double) = &Clmbr::crR ;
	void (Clmbr::*cr4) (void) = &Clmbr::crR ;


	class_<Clmbr>( "Cpp_Clmbr" )
	
	    // expose the constructors
	    .constructor< NumericVector, NumericMatrix, int, NumericMatrix, bool, bool >()    
	    .constructor< NumericVector, NumericMatrix, int, bool >()

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
		.method( "mle", &Clmbr::MLE , "printout maximum likelihood estimates of parameters" )
		.method( "param", &Clmbr::PARAM , "return maximum likelihood estimates of parameters" )
		.method( "sety", &Clmbr::SET_rWy , 
			"set new values for inverse-root of Sigma matrix times y-vector" )
	;

}                     


