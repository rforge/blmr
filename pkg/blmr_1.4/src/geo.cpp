//


#include "blmr.h"





double Cblmr::geo(const double th_2, double *const err)  const
// calculate probability that  gam*v > w  on the interval (th0, th2) or (th2, th0)
// given  gam0*v = z ,  using KSZ's geometric formula, 
{
	if(err!=0) *err= 0.;
	double  th2= th_2,  pr= 0.;
	if (th2 > x_d[n_d-2]) th2 = x_d[n_d-2]; else { if (th2 <x_d[k1]) th2 = x_d[k1]; }
	if (   variance_unknown  &&    cov_matrix_diagonal)  pr= geo_vu_D(th2,err);
	if (   variance_unknown  &&  ! cov_matrix_diagonal)  pr= geo_vu_ND(th2,err);
	if ( ! variance_unknown  &&    cov_matrix_diagonal)  pr= geo_vk_D(th2,err);
	if ( ! variance_unknown  &&  ! cov_matrix_diagonal)  pr= geo_vk_ND(th2,err);
	return min(pr,1.);
}





double Cblmr::geo_vu_D(const double th2, double *const err)  const
// for case variance unknown,  cov matrix diagonal
{
	if (fabs(th0-th2) < zero_eq || (th2<th0 && th0<x_d[k1]) || (th2>th0 && th0>x_d[n_d-2])) return 0.;

	const double r =sqrt((1-w*w)*(1-z*z)), rhoL =z*w-r, rhoZero =z/w, rhoU =z*w+r;

	int k2=0; 
	while (x_d[k2]<=th2 && k2<n_d) k2++; 
	const double ro2 = rho(th2,k2);
	if (ro2 > rhoU)  return 0.;

	double  pr=0,  arg;
	if(ro2<rhoZero) arg =(w*w-z*z)/(1-z*z);	else arg =(w-z*ro2)*(w-z*ro2)/(1-ro2*ro2)/(1-z*z);
	pr += F(m-2,-sqrt(arg));

// os = offset for kappa in function calls
	int kinc;
	if(th2<th0)  kinc= -1;  else  kinc= 1;
	const int  os= (1-kinc)/2;	

// x_d[k] = boundary of th0's data interval
	int  k= k0-os; 
	if(th0==x_d[k0-1]) k -=os; 
	if(th0<x_d[k1])  k= k1+1; 
	if(th0>x_d[n_d-2])  k= n_d-3;

	if ( (th2-x_d[k])*kinc < zero_eq  ) return pr;

	if ( mu_by_tau(th2,k2-1+os) > 1.- zero_eq )  return pr;		// assumes th2=x_d[j], a data point

// find th1  such that  mu_by_tau(th1) is > and near 1
	double th1;
	while ( mu_by_tau(x_d[k+kinc],k+1-os) > 1.-zero_eq ) k+=kinc;
	if ( mu_by_tau( x_d[k], k+1-os) < 1.+zero_eq)  th1 =x_d[k];
		else  th1 = bisect( x_d[k], x_d[k+kinc], &Cblmr::mu_by_tau, k+1-os, 1., inc_x);

// integrate piecewise from th1 to thZero
	double thZero;
	if ( ro2 > rhoZero-zero_eq ) thZero =th2;  else  if(rho(x_d[k],k)<rhoZero) thZero=x_d[k]; else {
		int kZ = k; 
		while(rho(x_d[kZ+kinc],kZ+kinc)>rhoZero && kZ+kinc>=0 && kZ+kinc<=n_d-1)  kZ += kinc; 
		thZero = rho_inv(rhoZero,kZ+1-os); 
	}

	double  error=0., er=0, *const per =&er,  x1 =th1,  x2 =x_d[k];
	while (  (x1-thZero)*kinc < -zero_eq )  {

		k += kinc;

		if(th2>th0) x2 =min(thZero,x_d[k]); else x2 =max(thZero,x_d[k]);

		pr += integrate( x1, x2, &Cblmr::Emupr, k+os, per);
		error += *per;

		x1 = x2;
	}
	if (x2!=x_d[k]) k -= kinc;

// integrate piecewise from thZero to thL
	while ( (x1-th2)*kinc < -zero_eq  &&  rho(x1,k+1-os) > rhoL + zero_eq )  {

		k += kinc;

		if ( mu_by_tau( x1, k+os) > -1. )  {

			if(th2>th0) x2 =min(th2,x_d[k]); else x2 =max(th2,x_d[k]);

			if ( mu_by_tau(x2, k+os) < -1.) x2 = bisect( x1, x2, &Cblmr::mu_by_tau, k+os, -1., -inc_x);

			pr += integrate( x1, x2, &Cblmr::Emupr, k+os, per);
			error += *per;
		}

		x1 = x_d[k];
	}

	if (err!=0) *err =error;
	return min(pr,1.);
}





double Cblmr::geo_vu_ND(const double th2, double *const err)  const
// case variance unknown,  cov matrix non-diagonal
{
	if (fabs(th0-th2) < zero_eq) return 0.;

	const double  r =sqrt((1-w*w)*(1-z*z)), rhoZero =z/w, rhoU =z*w+r;

// os = offset for kappa in function calls
	int kinc;  
	if(th2<th0) kinc =-1; else kinc =1;  
	const int  os = (1-kinc)/2;	

// x_d[k] = boundary of th0's data interval
	int  k = k0-os; 
	if(th0==x_d[k0-1])  k -=os; 
	if(th0<x_d[k1])  k= k1+1; 
	if(th0>x_d[n_d-2])  k= n_d-3;

// on th0's data interval, rho monotonic decreasing from rho=1
	double  pr=0.,  ro =rho(x_d[k],k);  if((th0-th2)*(x_d[k]-th2)<0.) ro =rho(th2);
	if ( ro < rhoU ) {
		if (ro < rhoZero) ro = rhoZero;
		const double arg = (w - z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
		pr +=  F(m-2,-arg);
	}

//determine next rho-monotonic interval, add intergral value
	double  error=0., er=0, *const per = &er;

	while ( (th2 - x_d[k])*kinc > zero_eq) {

		double a=x_d[k],b=x_d[k+kinc]; if ( (b-th2)*kinc >0) b=th2;

		const double rZp = a0[k+1-os]/b0[k+1-os];

// fabs(a-rZp)>zero_eq && fabs(b-rZp)>zero_eq && 

		if ( (a-rZp)*(b-rZp) < 0.) {
			pr += geo_vu_covNDab( k+1-os, a,  rZp,  -kinc,  per);
			error += *per;
			pr += geo_vu_covNDab( k+1-os, rZp,  b,  kinc,  per);
			error += *per;
		} else {
			pr += geo_vu_covNDab( k+1-os, a,  b,  1,  per);
			error += *per;
		}

		k += kinc;
	}


	if (err!=0) *err = error;
	return min(pr,1.);
}



double Cblmr::geo_vu_covNDab(const int k, const double a, const double b, const int lo, double *const err)  const
//get integral on rho-monotonic interval (a,b)
//assume 'a' is nearer th0
{
	if (err!=0) *err = 0.;
	if( fabs(b-a) < zero_eq ) return 0.;

	const double r =sqrt((1-w*w)*(1-z*z)), rhoU =z*w+r, rhoZero =z/w, rhoL =z*w-r, 
		rhoa =rho(a,k), rhob =rho(b,k);

	if(err!=0) *err= 0.;

	double Fa= 0., Fb= 0., ro, arg;

	if (rhoa > rhob) {

		if ( rhoa < rhoL || rhob > rhoU ) return 0.;

		if (rhoa > rhoZero) {
			if (rhoa > rhoU) Fa = 1.; else {
				ro= rhoa;
				arg = (w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
				Fa = F(m-2,arg);
			}
			if (rhob < rhoZero) ro = rhoZero; else ro= rhob;
			arg = (w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
			Fb = F(m-2,arg);
		}

	}  else  {

		if ( rhob < rhoL || rhoa > rhoU ) return 0.;

		if (rhoa < rhoZero) {
			if (rhoa < rhoL) Fa = 1.; else {
				ro= rhoa;
				arg = (w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
				Fa = F(m-2,arg);
			}
			if (rhob > rhoZero) ro = rhoZero; else ro= rhob;
			arg = (w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
			Fb = F(m-2,arg);
		}
	}

	double  pr =  Fa - Fb;

	if( fabs(pr) < zero_eq)  pr= 0.;

	const double  mbta= mu_by_tau(a,k), mbtb= mu_by_tau(b,k);	// mu_by_tau is always decreasing

	if ( mbta <-1+zero_eq || mbtb >1-zero_eq)  return pr;

	bool zero;  
	if ( (rhoa - rhoZero)*(rhob - rhoZero) < 0.) zero =true; else  zero =false;
	double thZero= 0.;
	if (zero) thZero = rho_inv(rhoZero,k,lo);

	double th1, thm1;
	if (mbta<1+zero_eq) th1 = a; else {
		if (zero) th1 = bisect(a,thZero,&Cblmr::mu_by_tau,k,1.,inc_x);
			else  th1 = bisect(a,b,&Cblmr::mu_by_tau,k,1.,inc_x);
	}
	if (mbtb>-1-zero_eq) thm1 = b; else {
		if (zero) thm1 = bisect(thZero,b,&Cblmr::mu_by_tau,k,-1.,-inc_x);
			else  thm1 = bisect(a,b,&Cblmr::mu_by_tau,k,-1.,-inc_x);
	}

	double  error=0., er=0, *const per = &er;
	if (zero) {
		pr += integrate( th1, thZero, &Cblmr::Emupr, k, per);
		error += *per;
		pr += integrate( thZero, thm1, &Cblmr::Emupr, k, per);
		error += *per;
	} else {
		pr += integrate( th1, thm1, &Cblmr::Emupr, k, per);
		error += *per;
	}

	if( fabs(pr) < zero_eq)  pr= 0.;


	if (err!=0) *err = error;
	return min(pr,1.);
}




double Cblmr::geo_vk_D(const double th2, double *const err)  const
// for case variance known,  cov matrix diagonal
{
	if (fabs(th0-th2) < zero_eq) return 0.;

	const double  rhoZero= z/w, ro2= rho(th2);

	const double  ro= max(rhoZero,ro2),  arg= (w - z*ro)/sqrt( 1-ro*ro );
	double  pr = Rf_pnorm5(-arg ,0,1,1,0);

// os = offset for kappa in function calls
	int kinc;  
	if(th2<th0) kinc= -1; else kinc= 1;
	const int  os= (1-kinc)/2;	

// x_d[k] = boundary of th0's data interval
	int k= k0-os; 
	if(th0==x_d[k0-1]) k -= os; 
	if(th0<x_d[k1]) k= k1+1; 
	if(th0>x_d[n_d-2]) k= n_d-3;

	if ( (th2-x_d[k])*kinc < zero_eq  )  return pr;	

// integrate piecewise from the th0 interval to thZero.  
// integrand negligible for abs(mu)/Omega > 7, which has 1 or no maxima.
	double thZero;
	if ( ro2 > rhoZero-zero_eq ) thZero= th2;  else  { 
		if( rho(x_d[k],k) < rhoZero ) thZero= x_d[k];  else  {
			int kZ= k; 
			while( rho(x_d[kZ+kinc],kZ+kinc) > rhoZero  &&  kZ+kinc >= 0  &&  kZ+kinc <= n_d-1 )  kZ += kinc; 
			thZero= rho_inv(rhoZero,kZ+1-os); 
		}
	}

	double  error=0., er=0, *const per =&er,  x1 =x_d[k],  x2 =x_d[k];
	while (  (x1-thZero)*kinc < zero_eq  &&  k+kinc >= 0  &&  k+kinc < n_d )  {

		if(th2>th0)  x2= min(thZero,x_d[k+kinc]);  else  x2= max(thZero,x_d[k+kinc]);

		const double a1= amu_by_Omega(x1, k+1-os), a2= amu_by_Omega(x2, k+1-os);

		if ( a1 <7.5 ||  a2 <7.5) {

			if(a1>7.5 && a2<6.5) x1= bisect( x1, x2, &Cblmr::amu_by_Omega, k+1-os, 7, inc_x);
			if(a1<6.5 && a2>7.5) x2= bisect( x1, x2, &Cblmr::amu_by_Omega, k+1-os, 7, inc_x);

			pr += integrate( x1, x2, &Cblmr::Emupr_vk, k+1-os, per);
			error += *per;
		}

		x1 = x_d[k+kinc];
		k += kinc;
	}
	if ( (x1-thZero)*kinc >= zero_eq )  { k -= kinc;  x1= thZero; }


// integrate piecewise from thZero to th2
// integrand negligible for abs(mu)/Omega > 7, which has 1 or no maxima.

	while ( (x1-th2)*kinc < zero_eq &&  k+kinc >= 0  &&  k+kinc < n_d )  {

		if(th2>th0) x2= min(th2,x_d[k+kinc]); else x2= max(th2,x_d[k+kinc]);

		const double a1= amu_by_Omega(x1, k+1-os), a2= amu_by_Omega(x2, k+1-os);

		if ( a1 <7.5 ||  a2 <7.5) {

			if(a1<6.5 && a2>7.5) x2= bisect( x1, x2, &Cblmr::amu_by_Omega, k+1-os, 7, inc_x);
			if(a2<6.5 && a1>7.5) x1= bisect( x1, x2, &Cblmr::amu_by_Omega, k+1-os, 7, inc_x);

			pr += integrate( x1, x2, &Cblmr::Emupr_vk, k+1-os, per);
			error += *per;
		}

		x1 = x_d[k+kinc];
		k += kinc;
	}

	if (err!=0) *err = error;
	return min(pr,1.);
}





double Cblmr::geo_vk_ND(const double th2, double *const err)  const
// case variance unknown,  cov matrix non-diagonal
{
	if (fabs(th0-th2) < zero_eq) return 0.;

// os = offset for kappa in function calls
	int kinc;  
	if(th2<th0)  kinc = -1;  else  kinc =1;
	const int  os = (1-kinc)/2;	

// x_d[k] = boundary of th0's data interval
	int k = k0-os; 
	if( th0 == x_d[k0-1] )  k -= os; 
	if( th0 < x_d[k1] )  k= k1+1; 
	if( th0 > x_d[n_d-2] )  k= n_d-3;

// on th0's data interval, rho is monotonic decreasing from rho(th0)=1
	const double  rhoZero= z/w;  
	double ro = rho(x_d[k]);
	if ( (th0-th2)*(x_d[k]-th2) < 0. ) ro = rho(th2);
	if (ro < rhoZero) ro = rhoZero;
	const double arg = (w - z*ro)/sqrt(1-ro*ro);
	double  pr =  Rf_pnorm5(-arg ,0,1,1,0) ;

//determine next rho-monotonic interval, add intergral value
	double  error=0., er=0, *const per = &er;

	while ( (th2 - x_d[k])*kinc > zero_eq) {

		double a= x_d[k], b= x_d[k+kinc]; if ( (b-th2)*kinc >0) b= th2;

		const double rZp = a0[k+1-os]/b0[k+1-os];

// fabs(a-rZp)>zero_eq && fabs(b-rZp)>zero_eq &&
		if ( (a-rZp)*(b-rZp) < 0.) {
			pr += geo_vk_covNDab( k+1-os, a,  rZp,  -kinc,  per);
			error += *per;
			pr += geo_vk_covNDab( k+1-os, rZp,  b,  kinc,  per);
			error += *per;
		} else {
			pr += geo_vk_covNDab( k+1-os, a,  b,  1,  per);
			error += *per;
		}

		k += kinc;
	}

	if (err!=0) *err = error;
	return min(pr,1.);
}



double Cblmr::geo_vk_covNDab(const int k, const double a, const double b, const int lo, double *const err)  const
//get integral on rho-monotonic interval (a,b)
// integrand negligible for abs(mu)/Omega > 7, which has 1 or no maxima.on rho>rhoZero and 
// likewise on rho<rhoZero
{
	if (err!=0) *err = 0.;
	if( fabs(b-a) < zero_eq ) return 0.;

	const double  rhoa = rho(a,k), rhob = rho(b,k), rhoZero= z/w;

//	if(rhoa < rhob)  { double tmp= a; a= b; b= tmp; tmp= rhoa; rhoa= rhob; rhob= tmp; }

	double ro, arga= 0., argb= 0., pr= 0.;

	if (rhoa > rhob)  {

		if (rhoa > rhoZero) {
			arga =(w-z*rhoa)/sqrt(1-rhoa*rhoa);
			if (rhob < rhoZero) ro = rhoZero; else ro= rhob;
			argb =(w-z*ro)/sqrt(1-ro*ro);
			pr =  Rf_pnorm5(arga ,0,1,1,0)   -   Rf_pnorm5(argb ,0,1,1,0) ;
		}

	}  else  {

		if (rhoa < rhoZero) {
			arga = (w-z*rhoa)/sqrt(1-rhoa*rhoa);
			if (rhob > rhoZero) ro = rhoZero; else ro= rhob;
			argb =(w-z*ro)/sqrt(1-ro*ro);
			pr =  Rf_pnorm5(arga ,0,1,1,0)   -   Rf_pnorm5(argb ,0,1,1,0) ;
		}

	}

	if( fabs(pr) < zero_eq)  pr= 0.;


	const double aa= amu_by_Omega(a,k), ab= amu_by_Omega(b,k);

	bool zero;  
	if ( (rhoa - rhoZero)*(rhob - rhoZero) < 0.)  zero = true;  else  zero = false;
	double thZero= 0.;
	if(zero)  thZero = rho_inv(rhoZero,k,lo);

	if(!zero && aa > 7. &&  ab > 7.)  return pr;

	double th1= a, thm1= b;
	if(zero) { 
		if(aa > 7.5) th1= bisect(a,thZero,&Cblmr::amu_by_Omega,k,7.,inc_x);
		if(ab > 7.5) thm1= bisect(thZero,b,&Cblmr::amu_by_Omega,k,7.,inc_x);
	} else {
		if( rhoa > rhoZero && aa > 7.5 && ab < 6.5)
			th1= bisect(a,b,&Cblmr::amu_by_Omega,k,7.,inc_x);
		if( rhoa < rhoZero && aa < 6.5 && ab > 7.5)
			thm1= bisect(a,b,&Cblmr::amu_by_Omega,k,7.,inc_x);
	}

	double  error=0., er=0, *const per = &er;
	if (zero) {
		pr += integrate( th1, thZero, &Cblmr::Emupr_vk, k, per);
		error += *per;
		pr += integrate( thZero, thm1, &Cblmr::Emupr_vk, k, per);
		error += *per;
	} else {
		pr += integrate( th1, thm1, &Cblmr::Emupr_vk, k, per);
		error += *per;
	}

	if( fabs(pr) < zero_eq)  pr= 0.;


	if (err!=0) *err = error;
	return min(pr,1.);
}


