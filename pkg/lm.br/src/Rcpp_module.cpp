//

#include "lmbr.h"



RCPP_MODULE(Clmbr){

	using namespace Rcpp ;


	void (Clmbr::*sl1) (int, double, double) = &Clmbr::sl1R ;
	void (Clmbr::*sl2) (int, double, double, double) = &Clmbr::sl2R ;
	NumericVector (Clmbr::*sl3) (int, int, int, double, double) = &Clmbr::sl3R ;
	NumericVector (Clmbr::*sl4) (int, int, int, double, double, double) = &Clmbr::sl4R ;

	void (Clmbr::*cr1) (double, int, double) = &Clmbr::cr1R ;
	NumericMatrix (Clmbr::*cr2) (double, int, double, int) = &Clmbr::cr2R ;


	class_<Clmbr>( "Cpp_Clmbr" )

// expose the constructor
	    .constructor< NumericVector, NumericMatrix, NumericMatrix, int, int, int >()

// expose methods
		.method( "sl", sl1 , "SL for th0 by specified method, accuracy" )
		.method( "sl", sl2 , "SL for (th0,a0) by specified method, accuracy" )
		.method( "sl", sl3 , "SL for th0 by specified method, accuracy and output flag" )
		.method( "sl", sl4 , "SL for (th0,a0) by specified method, accuracy and output flag" )
		.method( "ci", &Clmbr::ciR , "confidence interval for theta by specified method" )
		.method( "cr", cr1 ,
			"confidence region for (theta,alpha) by specified method and increment" )
		.method( "cr", cr2 , 
			"confidence region for (theta,alpha) by specified method and increment, return matrix" )
		.method( "mle", &Clmbr::MLE , "printout maximum likelihood estimates of parameters" )
		.method( "param", &Clmbr::PARAM , "return maximum likelihood estimates of parameters" )
		.method( "sety", &Clmbr::SET_rWy , 
			"reset values for square-root of weights times y-vector" )
	;

}

