//


#include "blmr.h"





double Cblmr::geo(double th2, double *err)
// calculate probability that  gam*v > w  on the interval (th0, th2) or (th2, th0)
// given  gam0*v = z ,  using KSZ's geometric formula, 
{
	if(err!=0) *err=0.;
	if (th2 > xs[ns-2]) th2 = xs[ns-2]; else { if (th2 <xs[k1]) th2 = xs[k1]; }
	double pr =0;
	if (   variance_unknown  &&    cov_matrix_diagonal)  pr= geo_vu_D(th2,err);
	if (   variance_unknown  &&  ! cov_matrix_diagonal)  pr= geo_vu_ND(th2,err);
	if ( ! variance_unknown  &&    cov_matrix_diagonal)  pr= geo_vk_D(th2,err);
	if ( ! variance_unknown  &&  ! cov_matrix_diagonal)  pr= geo_vk_ND(th2,err);
	return min(pr,1.);
}





double Cblmr::geo_vu_D(const double th2, double *err)
// for case variance unknown,  cov matrix diagonal
{
	if (fabs(th0-th2) < zero_eq || (th2<th0 && th0<xs[k1]) || (th2>th0 && th0>xs[ns-2])) return 0.;

	double r =sqrt((1-w*w)*(1-z*z)), rhoL =z*w-r, rhoZero =z/w, rhoU =z*w+r;

	int k2=0; while (xs[k2+1]<th2+zero_eq) k2++; double ro2 = rho(th2,k2);
	if (ro2 > rhoU)  return 0.;

	double  pr=0,  arg;
	if(ro2<rhoZero) arg =(w*w-z*z)/(1-z*z);	else arg =(w-z*ro2)*(w-z*ro2)/(1-ro2*ro2)/(1-z*z);
	pr += F(m-2,-sqrt(arg));

// os = offset for kappa in function calls
	int kinc, os;  {if(th2<th0) kinc =-1; else kinc =1;}  os =(1-kinc)/2;	

// xs[k] = boundary of th0's data interval
	int k =gk0-os+1; if(th0==xs[gk0]) k -=os; if(th0<xs[k1]) k=k1+1; if(th0>xs[ns-2]) k=ns-3;

	if ( (th2-xs[k])*kinc < zero_eq  ) return pr;

	if (mu_by_tau(th2,k2-1+os) > 1.-zero_eq)  return pr;

// find th1  such that  mu_by_tau(th1) is > and near 1
	double th1;
	while ( mu_by_tau(xs[k+kinc],k-os) > 1.-zero_eq ) k+=kinc;
	if ( mu_by_tau( xs[k], k-os) < 1.+zero_eq)  th1 =xs[k];
		else  th1 = bisect( xs[k], xs[k+kinc], &Cblmr::mu_by_tau, k-os, 1., inc_x);

// integrate piecewise from th1 to thZero
	double thZero;
	if ( ro2 > rhoZero-zero_eq ) thZero =th2;  else  if(rho(xs[k],k)<rhoZero) thZero=xs[k]; else
		{int kZ=k; while(rho(xs[kZ+kinc],kZ+kinc)>rhoZero) kZ+=kinc; thZero =rho_inv(rhoZero,kZ-os); }

	double  error=0., er=0, *per =&er,  x1 =th1,  x2 =xs[k];
	while (  (x1-thZero)*kinc < -zero_eq )  {

		k += kinc;

		if(th2>th0) x2 =min(thZero,xs[k]); else x2 =max(thZero,xs[k]);

		pr += integrate( x1, x2, &Cblmr::Emupr, k-1+os, per);
		error += *per;

		x1 = x2;
	}
	if (x2!=xs[k]) k -= kinc;

// integrate piecewise from thZero to thL
	while ( (x1-th2)*kinc < -zero_eq  &&  rho(x1,k) > rhoL + zero_eq )  {

		k += kinc;

		if ( mu_by_tau( x1, k-1+os) > -1. )  {

			if(th2>th0) x2 =min(th2,xs[k]); else x2 =max(th2,xs[k]);

			if ( mu_by_tau(x2, k-1+os) < -1.) x2 = bisect( x1, x2, &Cblmr::mu_by_tau, k-1+os, -1., -inc_x);

			pr += integrate( x1, x2, &Cblmr::Emupr, k-1+os, per);
			error += *per;
		}

		x1 = xs[k];
	}

	if (err!=0) *err =error;

	return min(pr,1.);
}





double Cblmr::geo_vu_ND(const double th2, double *err)
// case variance unknown,  cov matrix non-diagonal
{
	if (fabs(th0-th2) < zero_eq) return 0.;

	double r =sqrt((1-w*w)*(1-z*z)), rhoZero =z/w, rhoU =z*w+r;

// os = offset for kappa in function calls
	int kinc, os;  {if(th2<th0) kinc =-1; else kinc =1;}  os =(1-kinc)/2;	

// xs[k] = boundary of th0's data interval
	int k =gk0-os+1; if(th0==xs[gk0]) k -=os; if(th0<xs[k1]) k=k1+1; if(th0>xs[ns-2]) k=ns-3;

// on th0's data interval, rho monotonic decreasing from rho=1
	double  pr=0.,  ro =rho(xs[k],k);  if((th0-th2)*(xs[k]-th2)<0.) ro =rho(th2);
	if ( ro < rhoU ) {
		if (ro < rhoZero) ro = rhoZero;
		double arg = (w - z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
		pr +=  1. -  F(m-2,arg);
	}

//determine next rho-monotonic interval, add intergral value
	double  error=0., er=0, *per = &er;

	while ( (th2 - xs[k])*kinc > zero_eq) {

		double a=xs[k],b=xs[k+kinc]; if ( (b-th2)*kinc >0) b=th2;

		double rZp = a0[k-os]/b0[k-os];

		if ( fabs(a-rZp)>zero_eq && fabs(b-rZp)>zero_eq && (a-rZp)*(b-rZp) < 0.) {
			pr += geo_vu_covNDab( k-os, a,  rZp,  -kinc,  per);
			error += *per;
			pr += geo_vu_covNDab( k-os, rZp,  b,  kinc,  per);
			error += *per;
		} else {
			pr += geo_vu_covNDab( k-os, a,  b,  1,  per);
			error += *per;
		}

		k += kinc;
	}


	if (err!=0) *err = error;

	return min(pr,1.);
}



double Cblmr::geo_vu_covNDab(int k, double a, double b, int lo, double *err)
//get integral on rho-monotonic interval (a,b)
//assume a is nearer th0
{
	double r =sqrt((1-w*w)*(1-z*z)), rhoU =z*w+r, rhoZero =z/w, rhoL =z*w-r, 
		rhoa =rho(a,k), rhob =rho(b,k), thZero =0;

	if(err!=0) *err=0.;

	double Fmax, Fmin, ro, arg, pr=0, error=0., er=0, *per = &er;

	if (rhoa > rhob) {

		if ( rhoa < rhoL || rhob > rhoU ) return 0.;

		if (rhoa > rhoZero) {
			if (rhoa > rhoU) Fmax = 1.; else {
				ro=rhoa;
				arg =(w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
				Fmax = F(m-2,arg);
			}
			if (rhob < rhoZero) ro = rhoZero; else ro=rhob;
			arg =(w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
			Fmin = F(m-2,arg);
			pr += Fmax - Fmin;
		}

	}  else  {

		if ( rhob < rhoL || rhoa > rhoU ) return 0.;

		if (rhob > rhoZero) {
			if (rhob > rhoU) Fmax = 1.; else {
				ro=rhob;
				arg =(w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
				Fmax = F(m-2,arg);
			}
			if (rhoa < rhoZero) ro = rhoZero; else ro=rhoa;
			arg =(w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
			Fmin = F(m-2,arg);
			pr += Fmax - Fmin;
		}
	}


	if (mu_by_tau(a,k)<-1+zero_eq || mu_by_tau(b,k)>1-zero_eq) return pr;

	bool zero;  if ( (rhoa - rhoZero)*(rhob - rhoZero) < 0.) zero =true; else  zero =false;
	if (zero) thZero = rho_inv(rhoZero,k,lo);

	double th1, thm1;
	if (mu_by_tau(a,k)<1+zero_eq) th1 = a; else {
		if (zero) th1 = bisect(a,thZero,&Cblmr::mu_by_tau,k,1.,inc_x);
			else  th1 = bisect(a,b,&Cblmr::mu_by_tau,k,1.,inc_x);
	}
	if (mu_by_tau(b,k)>-1-zero_eq) thm1 = b; else {
		if (zero) thm1 = bisect(thZero,b,&Cblmr::mu_by_tau,k,-1.,-inc_x);
			else  thm1 = bisect(a,b,&Cblmr::mu_by_tau,k,-1.,-inc_x);
	}

	if (zero) {
		pr += integrate( th1, thZero, &Cblmr::Emupr, k, per);
		error += *per;
		pr += integrate( thZero, thm1, &Cblmr::Emupr, k, per);
		error += *per;
	} else {
		pr += integrate( th1, thm1, &Cblmr::Emupr, k, per);
		error += *per;
	}


	if (err!=0) *err = error;

	return min(pr,1.);
}




double Cblmr::geo_vk_D(const double th2, double *err)
// for case variance known,  cov matrix diagonal
{
	if (fabs(th0-th2) < zero_eq) return 0.;

	double  rhoZero = z/w, ro2 = rho(th2);

	double pr=0,  ro =max(rhoZero,ro2),  arg = (w - z*ro)/sqrt( 1-ro*ro );
	pr += 1.-  Rf_pnorm5(arg ,0,1,1,0);

// os = offset for kappa in function calls
	int kinc, os;  {if(th2<th0) kinc =-1; else kinc =1;}  os =(1-kinc)/2;	

// xs[k] = boundary of th0's data interval
	int k =gk0-os+1; if(th0==xs[gk0]) k -=os; if(th0<xs[k1]) k=k1+1; if(th0>xs[ns-2]) k=ns-3;

	if ( (th2-xs[k])*kinc < zero_eq  ) return pr;	

// integrate piecewise from the th0 interval to thZero
	double thZero;
	if ( ro2 > rhoZero-zero_eq ) thZero =th2;  else  if(rho(xs[k],k)<rhoZero) thZero=xs[k]; else
		{int kZ=k; while(rho(xs[kZ+kinc],kZ+kinc)>rhoZero) kZ+=kinc; thZero =rho_inv(rhoZero,kZ-os); }

	double  error=0., er=0, *per =&er,  x1 =xs[k],  x2 =xs[k];
	while (  (x1-thZero)*kinc < -zero_eq )  {

		k += kinc;

		if(th2>th0) x2 =min(thZero,xs[k]); else x2 =max(thZero,xs[k]);

		pr += integrate( x1, x2, &Cblmr::Emupr_vk, k-1+os, per);
		error += *per;

		x1 = x2;
	}
	if (x2!=xs[k]) k -= kinc;

// integrate piecewise from thZero to th2
	while ( (x1-th2)*kinc < -zero_eq )  {

		k += kinc;

		if(th2>th0) x2 =min(th2,xs[k]); else x2 =max(th2,xs[k]);

		pr += integrate( x1, x2, &Cblmr::Emupr_vk, k-1+os, per);
		error += *per;

		x1 = xs[k];
	}

	if (err!=0) *err = error;
	return min(pr,1.);
}





double Cblmr::geo_vk_ND(const double th2, double *err)
// case variance unknown,  cov matrix non-diagonal
{
	if (fabs(th0-th2) < zero_eq) return 0.;

// os = offset for kappa in function calls
	int kinc, os;  {if(th2<th0) kinc =-1; else kinc =1;}  os =(1-kinc)/2;	

// xs[k] = boundary of th0's data interval
	int k =gk0-os+1; if(th0==xs[gk0]) k -=os; if(th0<xs[k1]) k=k1+1; if(th0>xs[ns-2]) k=ns-3;

// on th0's data interval, rho is monotonic decreasing from rho=1
	double  rhoZero=z/w, pr=0.,  ro = rho(xs[k]);
	if ( (th0-th2)*(xs[k]-th2) < 0. ) ro = rho(th2);
	if (ro < rhoZero) ro = rhoZero;
	double arg = (w - z*ro)/sqrt(1-ro*ro);
	pr +=  1. -  Rf_pnorm5(arg ,0,1,1,0) ;

//determine next rho-monotonic interval, add intergral value
	double  error=0., er=0, *per = &er;

	while ( (th2 - xs[k])*kinc > zero_eq) {

		double a=xs[k],b=xs[k+kinc]; if ( (b-th2)*kinc >0) b=th2;

		double rZp = a0[k-os]/b0[k-os];

		if ( fabs(a-rZp)>zero_eq && fabs(b-rZp)>zero_eq && (a-rZp)*(b-rZp) < 0.) {
			pr += geo_vk_covNDab( k-os, a,  rZp,  -kinc,  per);
			error += *per;
			pr += geo_vk_covNDab( k-os, rZp,  b,  kinc,  per);
			error += *per;
		} else {
			pr += geo_vk_covNDab( k-os, a,  b,  1,  per);
			error += *per;
		}

		k += kinc;
	}

	if (err!=0) *err = error;
	return min(pr,1.);
}



double Cblmr::geo_vk_covNDab(int k, double a, double b, int lo, double *err)
//get integral on rho-monotonic interval (a,b)
{

	double  rhoa =rho(a,k), rhob =rho(b,k), rhoZero=z/w, thZero;

	bool zero;  if ( (rhoa - rhoZero)*(rhob - rhoZero) < 0.) zero =true; else  zero =false;
	if (zero) thZero = rho_inv(rhoZero,k,lo);


	double ro, arga, argb, pr=0, error=0., er=0, *per = &er;

	if (rhoa > rhob) {

		if (rhoa > rhoZero) {
			arga =(w-z*rhoa)/sqrt(1-rhoa*rhoa);
			if (rhob < rhoZero) ro = rhoZero; else ro=rhob;
			argb =(w-z*ro)/sqrt(1-ro*ro);
			pr +=  Rf_pnorm5(arga ,0,1,1,0)   -   Rf_pnorm5(argb ,0,1,1,0) ;
		}

	}  else  {

		if (rhob > rhoZero) {
			argb =(w-z*rhob)/sqrt(1-rhob*rhob);
			if (rhoa < rhoZero) ro = rhoZero; else ro=rhoa;
			arga =(w-z*ro)/sqrt(1-ro*ro);
			pr += Rf_pnorm5(argb ,0,1,1,0)   -   Rf_pnorm5(arga ,0,1,1,0) ;
		}
	}

	if (zero) {
		pr += integrate( a, thZero, &Cblmr::Emupr_vk, k, per);
		error += *per;
		pr += integrate( thZero, b, &Cblmr::Emupr_vk, k, per);
		error += *per;
	} else {
		pr += integrate( a, b, &Cblmr::Emupr_vk, k, per);
		error += *per;
	}


	if (err!=0) *err = error;

	return min(pr,1.);
}


