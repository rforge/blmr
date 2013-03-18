//
//  includes and global variable definitions for class Cblmr


#if !defined Cblmr_g		//prevents compiler from repeating this code in other files
#define Cblmr_g

#define R_NO_REMAP

#include "jama125/jama_eig.h"
#include "jama125/jama_cholesky.h"
#include "jama125/jama_lu.h"
#include "jama125/tnt.h"

#include <Rcpp.h>

#include <R.h>

#include <iomanip>
#include <algorithm>
#include <limits>
#include <math.h>

#ifdef ENABLE_NLS
   #include <libintl.h>
   #define _(String) dgettext ("blmr", String)
#else
   #define _(String) (String)
#endif


using namespace JAMA;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::Rcout;
using Rcpp::stop;


enum MODEL{ M1, M2 };
enum METHOD { GEO, GEO2, AF, AF2, MC, INIT };
const double zero_eq = ldexp(1.,-40);
const double inf = numeric_limits<double>::infinity();     
const double und = numeric_limits<double>::quiet_NaN();


#endif

