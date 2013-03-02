//

#include "blmr.h"




const double rmax = static_cast<double>(RAND_MAX + 2);


double Cblmr::unif_rand(void) const
//returns pseudo-random number between 0 and 1, excluding exactly 0 or 1
{
	return static_cast<double>(rand()+1)/rmax;
}




double Cblmr::rnorm(void) const
{
	#define BIG 134217728	// 2^27  unif_rand() alone is not of high enough precision
	double u1 = unif_rand();
	u1 = static_cast<int>(BIG*u1) + unif_rand();
	return Rf_qnorm5(u1/BIG ,0,1,1,0) ;
}


