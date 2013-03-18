//

#include "blmr.h"



double Cblmr::bisect(double x1, double x2, double (Cblmr::*fn)(double,int), int k, double value, double crit)
// use bisection to find x such that value < fn(x) < value + crit
// if  crit < 0  find x such that value - crit < fn(x) < value
{
	double  f1 = (this->*fn)(x1,k) - value,  f2 = (this->*fn)(x2,k) - value;
	if (f1*f2>0 || f1==f2 || isnan(f1*f2)) {
		Rcout << _("input to 'bisect' routine are the boundary points  x1, x2, f(x1), f(x2) = ") << endl; 
		Rcout << x1 << " " << x2 << " " << f1 << " " << f2 << endl;
		stop( _("'bisect' cannot find interim point from starting values") );
	}
	int iteration=0;
	while (  fabs(x1-x2) > fabs(crit)  && (iteration<50) ) {
		double  xmean = (x1+x2)/2,  fx = (this->*fn)(xmean,k)-value;
		if(f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; } 
		iteration++;
	}
	if(iteration==50)  Rf_warning( _("'bisect' failed to reach tolerance after 50 iterations") );
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}




double Cblmr::bisect_sl(double x1, double x2, METHOD mode, double crit)
// use bisection to find x such that   SL < sl(x) < SL + crit
// if  crit < 0  find x such that SL - crit < sl(x) < SL
{
	double f1= sl(x1,mode,false)-SL, f2= sl(x2,mode,false)-SL;
	if( fabs(f1)<zero_eq && fabs(f1-f2)<zero_eq) return (x1+x2)/2.;
	double p12 = f1*f2;
	if (p12>0 || f1==f2 || fabs(p12)>1 || isnan(p12)) {
		Rcout << _("input to 'bisect_sl' routine are the boundary points  x1, x2, sl(x1), sl(x2) =") << endl;
		Rcout << x1 << " " << x2 << " " << f1+SL << " " << f2+SL << endl;
		stop( _("'bisect_sl' cannot find interim point from starting values") );
	}
	int iteration=0;
	while (  fabs(x1-x2) > fabs(crit)  && (iteration<50) ) {
		double  xmean= (x1+x2)/2,  fx= sl(xmean,mode,false) -SL;
		if (f1*fx<=0 && f1!=fx) { x2= xmean; f2= fx; }  else  { x1= xmean; f1= fx; }
		iteration++;
	}
	if(iteration==50)  Rf_warning( _("'bisect_sl' failed to reach tolerance after 50 iterations") );
	if (crit<0) { if (f1 <= 0) return x1; else return x2; }
		else  { if (f1 >= 0) return x1; else return x2; }
}




double Cblmr::integrate(double x1, double x2, double (Cblmr::*fn)(double, int), int k, double *er)
// estimate integral of fn(x)dx for  x = x1 to x2  by composite Simpson's rule
// based on sects 4.3, 4.4 and error estimate based on sect. 4.6  in
// Burden,R.L. and Faires,J.D. (2011) "Numerical Analysis, 9th ed," Boston: Brooks/Cole.
{

	double a= min(x1,x2), b= max(x1,x2), abs_error= acc_sl_abs/2/n, rel_error= acc_sl_rel/2;

	if (er != 0) *er=0.;
	if (b-a < zero_eq) return 0.;

// initial estimate
	double h = (b-a)/2;
	double S = 2*(this->*fn)(a+h,k);
	double sum =  (this->*fn)(a,k)  +  (this->*fn)(b,k)  +  S;
	double estimate = (sum + S)*h/3;
	double error = fabs(estimate) + fabs(abs_error) + 1.;
	int intervals = 2;

	while (intervals < 513 && error > abs_error  &&  fabs(error/estimate) > rel_error) {

		h = h/2;
		intervals *= 2;
		S = 0.;
		for(int i=1; i < intervals; i+=2 ) S += (this->*fn)(a+i*h,k);
		S *= 2;
		sum += S;
		double prev_estimate = estimate;
		estimate = (sum + S)*h/3;
		error = fabs(estimate - prev_estimate)/8.;
	}

	if(intervals==512) Rf_warning( _("'integrate' failed to reach tolerance after maximum iterations") );

	if (er != 0) *er= error;

	return estimate;
}







#include "qng.h"


double rescale_error (double err, const double result_abs, const double result_asc) ;


double Cblmr::integrate2(double par_a, double par_b, double (Cblmr::*fn)(double, double*), double *err2)
//  estimate integral of fn(x)dx for  x = x1 to x2  by Gauss-Konrod 21, 43 or 87 point rule 
//  or, if that fails to reach tolerance, by adaptive Gauss-Konrod 21 point rule.
//  The following routine modifies the GNU Scientific Library "gsl_integration_qng" 2007 Brian Gough
//  which is a C-language translation of the QUADPACK  "QNG" and "QAG"  Fortran routines by Piessens et al.
//  Here in the R package "blmr", this routine integrates another numerical integral, so it adds
//  the average error estimate of the integrand to the error estimate of the numeric integral.
{
  double a =min(par_a,par_b), b =max(par_a,par_b), epsabs =acc_sl_abs, epsrel = acc_sl_rel, ferr =0, er=0., *per=&er;
  if (!variance_unknown) epsabs /= lambda;
  if (model==M2)  epsabs /= 2; 

  if (err2 != 0) *err2=0.;
  if (b-a < zero_eq) return 0.;


  double fv1[5], fv2[5], fv3[5], fv4[5];
  double savfun[21];  // array of function values which have been computed 
  double res10, res21, res43, res87;    // 10, 21, 43 and 87 point results 
  double result_kronrod, err ; 
  double resabs; // approximation to the integral of abs(f) 
  double resasc; // approximation to the integral of abs(f-i/(b-a)) 

  const double half_length =  0.5 * (b - a);
  const double abs_half_length = fabs(half_length);
  const double center = 0.5 * (b + a);
  const double f_center = (this->*fn)(center,per);
			ferr += *per;


// Compute the integral using the 10- and 21-point formula. 

  res10 = 0;
  res21 = w21b[5] * f_center;
  resabs = w21b[5] * fabs(f_center);

  for (int k = 0; k < 5; k++)
    {
      const double abscissa = half_length * x1[k];
      const double fval1 = (this->*fn)(center + abscissa,per);  
			ferr += *per;
      const double fval2 = (this->*fn)(center - abscissa,per);  
			ferr += *per;
      const double fval = fval1 + fval2;
      res10 += w10[k] * fval;
      res21 += w21a[k] * fval;
      resabs += w21a[k] * (fabs (fval1) + fabs (fval2));
      savfun[k] = fval;
      fv1[k] = fval1;
      fv2[k] = fval2;
    }

  for (int k = 0; k < 5; k++)
    {
      const double abscissa = half_length * x2[k];
      const double fval1 = (this->*fn)(center + abscissa,per);  
			ferr += *per;
      const double fval2 = (this->*fn)(center - abscissa,per);  
			ferr += *per;
      const double fval = fval1 + fval2;
      res21 += w21b[k] * fval;
      resabs += w21b[k] * (fabs (fval1) + fabs (fval2));
      savfun[k + 5] = fval;
      fv3[k] = fval1;
      fv4[k] = fval2;
    }

  resabs *= abs_half_length ;

  { 
    const double mean = 0.5 * res21;
  
    resasc = w21b[5] * fabs (f_center - mean);
    
    for (int k = 0; k < 5; k++) 
		resasc += (w21a[k] * (fabs (fv1[k] - mean) + fabs (fv2[k] - mean))
          + w21b[k] * (fabs (fv3[k] - mean) + fabs (fv4[k] - mean)));

    resasc *= abs_half_length ;
  }

  result_kronrod = res21 * half_length;
  
  err = rescale_error ((res21 - res10) * half_length, resabs, resasc) + (b-a)*ferr/21;

//   test for convergence. 
  if (err < epsabs || err < epsrel * fabs (result_kronrod))
    {
 	  if (err2 != 0) *err2= err;
      return result_kronrod ;
    }


// compute the integral using the 43-point formula. 

  res43 = w43b[11] * f_center;

  for (int k = 0; k < 10; k++)  res43 += savfun[k] * w43a[k];

  for (int k = 0; k < 11; k++)
    {
      const double abscissa = half_length * x3[k];
      const double fval1 = (this->*fn)(center + abscissa,per);  
			ferr += *per;
      const double fval2 = (this->*fn)(center - abscissa,per);
			ferr += *per;
      const double fval = fval1 + fval2;
      res43 += fval * w43b[k];
      savfun[k + 10] = fval;
    }

  result_kronrod = res43 * half_length;
  err = rescale_error ((res43 - res21) * half_length, resabs, resasc) + (b-a)*ferr/43;

//   test for convergence. 
  if (err < epsabs || err < epsrel * fabs (result_kronrod))
    {
 	  if (err2 != 0) *err2= err;
      return result_kronrod ;
    }


// compute the integral using the 87-point formula. 

  res87 = w87b[22] * f_center;

  for (int k = 0; k < 21; k++)  res87 += savfun[k] * w87a[k];

  for (int k = 0; k < 22; k++)
    {
      const double abscissa = half_length * x4[k];
      const double fval1 = (this->*fn)(center + abscissa,per);
			ferr += *per;
      const double fval2 = (this->*fn)(center - abscissa,per);
			ferr += *per;
      const double fval = fval1 + fval2;
	  res87 += w87b[k] * fval;
    }

  result_kronrod = res87 * half_length ;
  err = rescale_error ((res87 - res43) * half_length, resabs, resasc) + (b-a)*ferr/87;

//   test for convergence. 
  if (err < epsabs || err < epsrel * fabs (result_kronrod))
    {
 	  if (err2 != 0) *err2= err;
      return result_kronrod ;
    }



// failed to converge by  21-43-87  composite Gauss-Kronrod-Patterson
// switch to adaptive quadrature, integrate subintervals by 21-point Gauss-Kronrod

//limit= maximum number of subintervals to be integrated
//ablist, intg, errs = lists of interval endpoints, integral values, and max. error estimates
  const int  limit=20;			
  double ablist[limit+1], intg[limit], errs[limit];	
  double intgsum=result_kronrod;

  ablist[0]=a;
  ablist[1]=b;
  intg[0]=result_kronrod;
  errs[0]=err;


// while loop:  bisect interval with largest error and integrate subintervals
// repeat until error fits condition or until reached maximum subintervals

  int numi=1;
  while( numi < limit  &&  err > epsabs  &&  err > epsrel*fabs(intgsum) ) {

	  int kmax=0;
	  for(int k=1;k<numi;k++) if(errs[k]>errs[kmax]) kmax=k;
	  for(int k=numi;k>kmax;k--) {
		  ablist[k+1]=ablist[k];
		  intg[k]=intg[k-1];
		  errs[k]=errs[k-1];
	  }
	  ablist[kmax+1] = (ablist[kmax] + ablist[kmax+1])/2;

	  double ai=ablist[kmax], bi=ablist[kmax+1];
	  intg[kmax] = qgk21(ai,bi, &Cblmr::ipr, per);
	  errs[kmax] =*per;

	  ai=bi; 
	  bi=ablist[kmax+2];
	  intg[kmax+1] = qgk21(ai,bi, &Cblmr::ipr, per);
	  errs[kmax+1] =*per;

	  numi++;
	  err =intgsum =0;
	  for(int k=0;k<numi;k++) { err += errs[k]; intgsum += intg[k]; }
  }

  if (err < epsabs || err < epsrel * fabs(intgsum)) {
 	  if (err2 != 0) *err2= err;
	  return intgsum ;
	}

  Rf_warning( _("'integrate2' failed to reach tolerance with highest-order rule nor with adaptive rule") );

  if (err2 != 0) *err2= err;
  return result_kronrod ;
}





double Cblmr::qgk21(double a, double b, double (Cblmr::*fn)(double, double*), double *err2)
//  estimate integral of fn(x)dx for  x = x1 to x2  by Gauss-Konrod 21 point rule
//  The following routine modifies the GNU Scientific Library "gsl_integration_qng"  2007 Brian Gough
//  which is a C-language translation of the QUADPACK  "QNG"  Fortran routine.
//  This routine adds the average error estimate of the integrand to the error estimate of this integral.
{
  double  ferr =0, er=0., *per=&er;

  if (err2 != 0)  *err2=0.;
  if (b-a < zero_eq)  return 0.;


  double fv1[5], fv2[5], fv3[5], fv4[5];
  double res10, res21;    // 10, 21 point results 
  double result_kronrod, err; 
  double resabs; // approximation to the integral of abs(f) 
  double resasc; // approximation to the integral of abs(f-i/(b-a)) 

  const double half_length =  0.5 * (b - a);
  const double abs_half_length = fabs(half_length);
  const double center = 0.5 * (b + a);
  const double f_center = (this->*fn)(center,per);
			ferr += *per;


// Compute the integral using the 10- and 21-point formula. 

  res10 = 0;
  res21 = w21b[5] * f_center;
  resabs = w21b[5] * fabs(f_center);

  for (int k = 0; k < 5; k++)
    {
      const double abscissa = half_length * x1[k];
      const double fval1 = (this->*fn)(center + abscissa,per);  
			ferr += *per;
      const double fval2 = (this->*fn)(center - abscissa,per);  
			ferr += *per;
      const double fval = fval1 + fval2;
      res10 += w10[k] * fval;
      res21 += w21a[k] * fval;
      resabs += w21a[k] * (fabs (fval1) + fabs (fval2));
      fv1[k] = fval1;
      fv2[k] = fval2;
    }

  for (int k = 0; k < 5; k++)
    {
      const double abscissa = half_length * x2[k];
      const double fval1 = (this->*fn)(center + abscissa,per);  
			ferr += *per;
      const double fval2 = (this->*fn)(center - abscissa,per);  
			ferr += *per;
      const double fval = fval1 + fval2;
      res21 += w21b[k] * fval;
      resabs += w21b[k] * (fabs (fval1) + fabs (fval2));
      fv3[k] = fval1;
      fv4[k] = fval2;
    }

  resabs *= abs_half_length ;

  { 
    const double mean = 0.5 * res21;
  
    resasc = w21b[5] * fabs (f_center - mean);
    
    for (int k = 0; k < 5; k++) 
		resasc += (w21a[k] * (fabs (fv1[k] - mean) + fabs (fv2[k] - mean))
          + w21b[k] * (fabs (fv3[k] - mean) + fabs (fv4[k] - mean)));

    resasc *= abs_half_length ;
  }

  result_kronrod = res21 * half_length;
  
  err = rescale_error ((res21 - res10) * half_length, resabs, resasc) + (b-a)*ferr/21;

  if (err2 != 0)  *err2=err;

  return  result_kronrod;

}





double
rescale_error (double err, const double result_abs, const double result_asc)
{
  err= fabs(err);

  if (result_asc != 0 && err != 0) {
        double scale= pow((200*err/result_asc), 1.5);
        if(scale < 1)  err= result_asc*scale; else  err= result_asc;
      }

  if (result_abs > DBL_MIN / (50*DBL_EPSILON))  {
      double min_err = 50 * DBL_EPSILON * result_abs ;
      if(min_err > err)  err = min_err;
    }
  
  return err;
}

