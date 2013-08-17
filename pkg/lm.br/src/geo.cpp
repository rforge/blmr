//


#include "lmbr.h"

void  igeo( double *const x,  const int n,  void *const ex );


//  'Rdqag' sometimes returns error code '5' for "doubtful convergence,"  
//  but the integral is convergent, O(theta^-2), and results have been reliable,  
//  so routines below ignore  ier = 5.




double Clmbr::geo(const double th2, double *const perr)  const
// calculate probability that  gam*v > w  on the interval (th0, th2) or (th2, th0)
// given  gam0*v = z ,  using KSZ's geometric formula, 
{
	if( perr != 0 )  *perr = 0.;

	double  pr= 0.;

	if (   variance_unknown  &&    cov_matrix_diagonal )   pr= geo_vu_D( th2, perr);
	if (   variance_unknown  &&  ! cov_matrix_diagonal )   pr= geo_vu_ND( th2, perr);
	if ( ! variance_unknown  &&    cov_matrix_diagonal )   pr= geo_vk_D( th2, perr);
	if ( ! variance_unknown  &&  ! cov_matrix_diagonal )   pr= geo_vk_ND( th2, perr);

	return  min( pr, 1. );
}




double Clmbr::geo_vu_D(const double th2, double *const err)  const
// for case variance unknown,  weights matrix diagonal
{
	if ( fabs(th0-th2) < zero_eq )  return 0.;


	const double  rad= sqrt((1-w*w)*(1-z*z)),  rU= z*w+rad,  rZ= z/w,  rL= z*w-rad,  r2= rho(th2);
//Rcout << "th0 th2 rU rZ rL r2  " << th0 << " " << th2 << " " << rU << " " << rZ << " " << rL << " " << r2 << endl;

	if(r2 > rU)  return 0.;

	double  arg;
	if( r2 < rZ )  arg= sqrt( (w*w-z*z)/(1-z*z) );  else  arg= (w-z*r2)/sqrt( (1-r2*r2)*(1-z*z) );
	double  pr =  F( m-2, -arg );
//Rcout << "pr " << pr << endl;


//   ki  is the first data-interval where  Omega != 0  after  th0 ,   kj  is the interval with outside boundary  th2
	int  ki, kj, kinc;
	if( th2 > th0 )   { ki= k0+1;  kj= ns-2;  kinc= 1;}   else   {ki= k0-1;  kj= k1+1;  kinc= -1;}
	if( th2 < th0 && k0 > 0 )  if( th0==xs[k0-1] )  ki= k0-2;
	if( k1 >= 0 )  if( th2 > th0 && th0 < xs[k1] )  ki= k1+2;
	if( th2 < th0 && th0 > xs[ns-2] )  ki= ns-3;


// integrate piecewise from ki to kj

	double  error= 0.;

	for( int k= ki; (k-kj)*kinc <= 0; k += kinc)  {

		double ra, rb;
		if( kinc > 0 )  {
			ra= rho(xs[k-1],k);
			rb= rho(xs[k],k);
		}  else  {
			ra= rho(xs[k],k);
			if(k>0)  rb= rho(xs[k-1],k);  else  rb= rho(-Inf,0);
		}


		bool  trunc= false;

//  check where  mu_by_tau  defined
		if( ra < rL )  break;
		if( rb > rU )  continue;
		if( ra > rU )  { ra= rU;  trunc= true; }
		if( rb < rL )  { rb= rL;  trunc= true; }


//  check where  -1 < mu_by_tau < 1
		const double  rd= sqrt((1-B[k])*(1-w*w))/w,  ru= rZ + rd,  rl= rZ - rd;

		if( ra < rl || rb > ru )  continue;

		if( ra > ru )  { ra= ru;  trunc= true; }
		if( rb < rl )  { rb= rl;  trunc= true; }

		double  th1,  th2;
		if( trunc )  { th1= rho_inv(ra,k);  th2= rho_inv(rb,k); }  else  
			{ if(k>0) {th1= xs[k-1]; th2= xs[k];}  else  {th1= -Inf; th2= xs[k];} }
		double  tha= max(th1,th2),  thb= min(th1,th2);

		bool  zero = false;
		if( (rZ-ra)*(rZ-rb) < 0 )  zero = true;
		double  thZ = NaN;
		if( zero )  thZ = rho_inv( rZ, k );
		if( fabs(thZ-tha) < zero_eq || fabs(thZ-thb) < zero_eq || isinf(thZ) )  zero= false;

//  'Rdqag' parameters
		int  inf_flag = -1,  neval =0,  ier =0,  limit =100,  lenw = 4*limit,  last =0,  iwork[limit];
		double  epsabs = acc_sl_abs/2/n,  epsrel = acc_sl_rel/2,  result =0,  abserr =0,  work[lenw];

		const void *  exc[2] = {  this,  &k  };
		void **  ex = const_cast< void** >( exc );


		if( zero ) {

			Rdqags( igeo, ex, &tha, &thZ, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );

			pr += fabs(result);
			error += abserr;
			if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 


//Rcout << "tha thZ pr error  " << tha << " " << thZ << " " << pr << " " << error << endl;

			if( isinf(thb) )  
				Rdqagi( igeo, ex, &thZ, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );
			else
				Rdqags( igeo, ex, &thZ, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );

			pr += fabs(result);
			error += abserr;
			if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 

//Rcout << "thZ thb pr error  " << thZ << " " << thb << " " << pr << " " << error << endl;

		} else {

			if( isinf(thb) )
				Rdqagi( igeo, ex, &tha, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );
			else
				Rdqags( igeo, ex, &tha, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );

			pr += fabs(result);
			error += abserr;
			if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 
//Rcout << "tha thb pr error  " << tha << " " << thb << " " << pr << " " << error << endl;

		}
//Rcout << "k ra rb pr error  " << k << " " << ra << " " << rb << " " << pr << " " << error << endl;
	}

	if (err!=0)  *err = error;
	return min(pr,1.);
}




double Clmbr::geo_vu_ND(const double th2, double *const err)  const
// case variance unknown,  weights matrix non-diagonal
{
	if ( fabs(th0-th2) < zero_eq )  return 0.;

	const double  rad= sqrt((1-w*w)*(1-z*z)),  rU= z*w + rad,  rZ= z/w;


//   ki  is the first data-interval where  Omega != 0  after  th0 ,   kj  is the interval with outside boundary  th2
	int  ki, kj, kinc;
	if( th2 > th0 )   ki= k0+1,  kj= ns-2,  kinc= 1;   else   ki= k0-1,  kj= k1+1,  kinc= -1;
	if( th2 < th0 && k0 > 0 )  if( th0==xs[k0-1] )  ki= k0-2;
	if( k1 >= 0 )  if( th2 > th0 && th0 < xs[k1] )  ki= k1+2;
	if( th2 < th0 && th0 > xs[ns-2] )  ki= ns-3;

//Rcout << endl << endl << "vkND:  th0 th2 ki kj rZ  " << th0 << " " << th2 << " " << ki << " " << kj << " " << rZ << endl;


// on th0's data interval,  rho  is monotonic decreasing from rho=1
	double  pr = 0.,  ro;
	if( th2 > th0 )  ro= rho( xs[ki-1], ki );  else  {
		if(ki<0)  ro= rho( -Inf, 0 );  else  ro= rho( xs[ki], ki );
	}
	if ( ro < rU ) {
		double  arg;
		if( ro < rZ )  arg= sqrt( (w*w-z*z)/(1-z*z) );  else  arg= (w-z*ro)/sqrt( (1-ro*ro)*(1-z*z) );
		pr +=  F( m-2, -arg );
	}
//Rcout << "ro ro2 arg  " << ro << " " << rho(th2) << " " << arg << endl;
//Rcout << "k0 pr  " << k0 << " " << pr << endl;


// integrate piecewise from ki to kj  on rho-monotonic sub intervals

	double  error= 0.;

	for( int k= ki; (k-kj)*kinc <= 0; k += kinc)  {

		double  tha,  thb;
		if( kinc > 0 )  tha= xs[k-1],  thb= xs[k];   else  {
			tha= xs[k];
			if(k>0)  thb= xs[k-1];  else  thb= -Inf;
		}

		const double  thZp = a0[k]/b0[k];	// theta  for  zero  of  derivative of 'rho'
		bool  pzero= false;
		if(k>0)  {  if(  xs[k-1] < thZp  &&  thZp < xs[k]  )  pzero= true;  }
		   else  {  if( thZp < xs[k] )  pzero= true;  }

		double  er=0,  *const per = &er;
		if(  pzero  )  {
			pr += geo_vu_NDab( k,  tha,  thZp,  -kinc,  per);
			error += *per;
			pr += geo_vu_NDab( k,  thZp,  thb,  kinc,  per);
			error += *per;
		}  else  {
			pr += geo_vu_NDab( k,  tha,   thb,  1,  per);
			error += *per;
		}
	}


	if ( err != 0 )  *err = error;
	return  min( pr, 1. );
}




double Clmbr::geo_vu_NDab(const int k, const double th_a, const double th_b, const int hilo, double *const err)  const
//get integral on rho-monotonic interval (a,b)  where  a  is nearer  'theta0'
{
	if (err!=0)  *err = 0.;

	if( fabs( th_a - th_b ) < zero_eq )  return 0.;

	const double  rad = sqrt((1-w*w)*(1-z*z)),  rU = z*w + rad,  rZ = z/w,  rL = z*w - rad;
	double  r1= rho(th_a,k),  r2= rho(th_b,k),  ra= max(r1,r2),  rb= min(r1,r2);
//if(k==2) Rcout << endl << "vuNDab:  k tha thb ra rb  " << k << " " << th_a << " " << th_b << " " << r1 << " " << r2 << endl;
	if ( ra < rL || rb > rU )  return 0.;


	double  F1= 0.,  F2= 0.,  arg;

	if (r1 > r2) {

		if (r1 > rZ) {
			if (r1 > rU) F1 = 1.;  else  {
				arg = (w-z*r1)/sqrt( (1-r1*r1)*(1-z*z) );
				F1 = F(m-2,arg);
			}
			if (r2 < rZ)  arg= sqrt( (w*w-z*z)/(1-z*z) );  else  arg = (w-z*r2)/sqrt( (1-r2*r2)*(1-z*z) );
			F2 = F(m-2,arg);
		}

	}  else  {

		if (r1 < rZ) {
			if (r1 < rL) F1 = 1.;  else  {
				arg = (w-z*r1)/sqrt( (1-r1*r1)*(1-z*z) );
				F1 = F(m-2,arg);
			}
			if (r2 > rZ)  arg= sqrt( (w*w-z*z)/(1-z*z) );  else  arg = (w-z*r2)/sqrt( (1-r2*r2)*(1-z*z) );
			F2 = F(m-2,arg);
		}
	}

	double  pr =  F1 - F2;


	bool  trunc= false;

//  check where  mu_by_tau  defined
	if( ra > rU )  ra= rU,  trunc= true;
	if( rb < rL )  rb= rL,  trunc= true;
//if(k==2) Rcout << "ra rb  " << ra << " " << rb << endl;


//  check where  -1 < mu_by_tau < 1
	const double  rd= sqrt((1-B[k])*(1-w*w))/w,  ru= rZ + rd,  rl= rZ - rd;

	if( ra < rl || rb > ru )  return pr;

	if( ra > ru )  ra= ru,  trunc= true;
	if( rb < rl )  rb= rl,  trunc= true;
//if(k==2) Rcout << "ra rb  " << ra << " " << rb << endl;


	double  th1,  th2;
	if( trunc )  th1= rho_inv(ra,k,hilo),  th2= rho_inv(rb,k,hilo);  else  th1= th_a,  th2= th_b;
	double  tha= max(th1,th2),  thb= min(th1,th2);

	bool  zero = false;
	if( (rZ-ra)*(rZ-rb) < 0 )  zero = true;
	double  thZ = NaN;
	if( zero )  thZ = rho_inv( rZ, k );
	if( fabs(thZ-tha) < zero_eq || fabs(thZ-thb) < zero_eq || isinf(thZ) )  zero= false;

//  'Rdqag' parameters
	int  inf_flag = -1,  neval =0,  ier =0,  limit =100,  lenw = 4*limit,  last =0,  iwork[limit];
	double  epsabs = acc_sl_abs/2/n,  epsrel = acc_sl_rel/2,  result =0,  abserr =0,  work[lenw];

	const void *  exc[2] = {  this,  &k  };
	void **  ex = const_cast< void** >( exc );

//Rcout << "k tha thb ra rb rZ pr  " << k << " " << tha << " " << thb << " " << ra << " " << rb << " " << rZ << " " << pr << endl;

	double  error= 0;
	if( zero )  {

		Rdqags( igeo, ex, &tha, &thZ, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
			&limit, &lenw, &last, iwork, work );

		pr += fabs(result);
		error += abserr;
		if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 


//Rcout << "tha thZ pr  " << tha << " " << thZ << " " << pr << endl;

		if( isinf(thb) )  
			Rdqagi( igeo, ex, &thZ, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );
		else
			Rdqags( igeo, ex, &thZ, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );

		pr += fabs(result);
		error += abserr;
		if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 

//Rcout << "thZ thb pr  " << thZ << " " << thb << " " << pr << endl;

	}  else  {

		if( isinf(thb) )
			Rdqagi( igeo, ex, &tha, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );
		else
			Rdqags( igeo, ex, &tha, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );

		pr += fabs(result);
		error += abserr;
		if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 
//Rcout << "tha thb pr  " << tha << " " << thb << " " << pr << endl;

	}


	if (err!=0)  *err = error;
	return  min( pr, 1. );
}




double  Clmbr::geo_vk_D( const double th2,  double *const err )  const
// for case variance known,  cov matrix diagonal
{
	if ( fabs(th0-th2) < zero_eq )  return 0.;

	const double  rZ= z/w,  r2= rho(th2);

	double  arg;
	if( r2 < rZ )  arg= sqrt(w*w-z*z);  else  arg= (w - z*r2)/sqrt( 1-r2*r2 );
	double  pr = Rf_pnorm5(-arg ,0,1,1,0);
//Rcout << "th0 th2 rZ r2 pr  " << th0 << " " << th2 << " " << rZ <<  " " << r2 << " " << pr << endl;


//   ki  is the first data-interval where  Omega != 0  after  th0 ,   kj  is the interval with outside boundary  th2
	int  ki, kj, kinc;
	if( th2 > th0 )   ki= k0+1,  kj= ns-2,  kinc= 1;   else   ki= k0-1,  kj= k1+1,  kinc= -1;
	if( th2 < th0 && k0 > 0 )  if( th0==xs[k0-1] )  ki= k0-2;
	if( k1 >= 0 )  if( th2 > th0 && th0 < xs[k1] )  ki= k1+2;
	if( th2 < th0 && th0 > xs[ns-2] )  ki= ns-3;


//  integrate piecewise from ki to kj

	double  error= 0.;

	for( int k= ki; (k-kj)*kinc <= 0; k += kinc)  {

		double  tha,  thb;
		if( kinc > 0 )  {
			tha= xs[k-1];
			thb= xs[k];
		}  else  {
			tha= xs[k];
			if(k>0)  thb= xs[k-1];  else  thb= -Inf;
		}

		const double  aa= amu_by_Omega( tha, k ),  ab= amu_by_Omega( thb, k );

// check for zero of 'mu'
		bool  zero = false;
		double  ra= rho(tha,k), rb= rho(thb,k);
		if( (rZ-ra)*(rZ-rb) < 0 )  zero = true;
		double  thZ = NaN;
		if( zero )  thZ = rho_inv( rZ, k );
		if( fabs(thZ-tha) < zero_eq || fabs(thZ-thb) < zero_eq || isinf(thZ) )  zero= false;
//Rcout << "k tha thb thZ aa ab aZ  " << k << " " << tha << " " << thb << " " << thZ << " " << aa << " " << ab << " " << amu_by_Omega( thZ, k ) << endl;  

		if( !zero  &&  aa > 6.5  &&  ab > 6.5 )  continue;

//  'Rdqag' parameters
		int  inf_flag = -1,  neval =0,  ier =0,  limit =100,  lenw = 4*limit,  last =0,  iwork[limit];
		double  epsabs = acc_sl_abs/2/n,  epsrel = acc_sl_rel/2,  result =0,  abserr =0,  work[lenw];

		const void *  exc[2] = {  this,  &k  };
		void **  ex = const_cast< void** >( exc );

		if( zero )  {	// mu = zero  on interval

			if( aa > 7.5 )  tha= bisect( tha, thZ, &Clmbr::amu_by_Omega, k, 7, inc_x);			
			if( ab > 7.5  &&  !isinf(thb) )  thb= bisect( thZ, thb, &Clmbr::amu_by_Omega, k, 7, inc_x);			

			Rdqags( igeo, ex, &tha, &thZ, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );

			pr += fabs(result);
			error += abserr;
			if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 


//Rcout << "tha thZ pr  " << tha << " " << thZ << " " << pr << endl;

			if( isinf(thb) )  
				Rdqagi( igeo, ex, &thZ, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );
			else
				Rdqags( igeo, ex, &thZ, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );

			pr += fabs(result);
			error += abserr;
			if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 

		}  else  {		// no zero 

			if(aa>7.5 && ab<6.5)  {
				if( !isinf(thb) )  
					tha= bisect( tha, thb, &Clmbr::amu_by_Omega, k, 7, inc_x);
				else  {
					double th= min( -1., tha);
					while( amu_by_Omega(th,k) > 6.8 )  th *= 2;
					tha= bisect( tha, th, &Clmbr::amu_by_Omega, k, 7, inc_x);
				}
			}
			if(aa<6.5 && ab>7.5 && !isinf(thb) )  thb= bisect( tha, thb, &Clmbr::amu_by_Omega, k, 7, inc_x);

//Rcout << "tha thb pr1  " << tha << " " << thb << " " << pr << endl;

			if( isinf(thb) )
				Rdqagi( igeo, ex, &tha, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );
			else
				Rdqags( igeo, ex, &tha, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
					&limit, &lenw, &last, iwork, work );

			pr += fabs(result);
			error += abserr;
			if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 
		}
//Rcout << "ki ra rb pr  " << ki << " " << ra << " " << rb << " " << pr << endl;
	}


	if (err!=0)  *err =error;
	return  min(pr,1.);
}





double Clmbr::geo_vk_ND(const double th2, double *const err)  const
// case variance unknown,  cov matrix non-diagonal
{
	if ( fabs(th0-th2) < zero_eq )  return 0.;

//   ki  is the first data-interval where  Omega != 0  after  th0 ,   kj  is the interval with outside boundary  th2
	int  ki, kj, kinc;
	if( th2 > th0 )   ki= k0+1,  kj= ns-2,  kinc= 1;   else   ki= k0-1,  kj= k1+1,  kinc= -1;
	if( th2 < th0 && k0 > 0 )  if( th0==xs[k0-1] )  ki= k0-2;
	if( k1 >= 0 )  if( th2 > th0 && th0 < xs[k1] )  ki= k1+2;
	if( th2 < th0 && th0 > xs[ns-2] )  ki= ns-3;

// on th0's data interval,  rho  is monotonic decreasing from rho=1
	const double  rZ= z/w;  
//Rcout << endl << endl << "vkND:  th0 th2 ki kj rZ w z  " << th0 << " " << th2 << " " << ki << " " << kj << " " << rZ <<
//		" " << w << " " << z << endl;
	double  ro,  arg;
	if( th2 > th0 )  ro= rho( xs[ki-1], ki );  else  {
		if(ki<0)  ro= rho( -Inf, 0 );  else  ro= rho( xs[ki], ki );
	}
	if( ro < rZ )  arg= sqrt( w*w - z*z );  else  arg= (w-z*ro)/sqrt( 1 - ro*ro );
//Rcout << "ro ro2 arg  " << ro << " " << rho(th2) << " " << arg << endl;
	double  pr =  Rf_pnorm5( -arg ,0,1,1,0) ;
//Rcout << "k0 pr  " << k0 << " " << pr << endl;



// integrate piecewise from ki to kj  on rho-monotonic sub intervals

	double  error= 0.;
	for( int k = ki; (k-kj)*kinc <= 0; k += kinc )  {

		double  tha,  thb;
		if( kinc > 0 )  tha= xs[k-1],  thb= xs[k];   else  {
			tha= xs[k];
			if(k>0)  thb= xs[k-1];  else  thb= -Inf;

		}

		const double  thZp= a0[k]/b0[k];	//  thZp = unique theta  for  zero of derivative of 'rho'
		bool  pzero= false;
		if(k>0)  {  if(  xs[k-1] < thZp  &&  thZp < xs[k]  )  pzero= true;  }
		   else  {  if( thZp < xs[k] )  pzero= true;  }


		double  er=0,  *const per = &er;
		if(  pzero  )  {
			pr += geo_vk_NDab( k,  tha,  thZp,  -kinc,  per);
			error += *per;
			pr += geo_vk_NDab( k,  thZp,  thb,  kinc,  per);
			error += *per;
		}  else  {
			pr += geo_vk_NDab( k,  tha,   thb,  1,  per);
			error += *per;
		}

//Rcout << endl << "k tha thb thZp pzero pr  " << k << " " << tha << " " << thb << " " << thZp << " " << pzero << " " << pr << endl;
	}


//Rcout << endl << "vkND:  pr  " << pr << endl;
	if (err!=0)  *err = error;
	return  min( pr, 1. );
}




double Clmbr::geo_vk_NDab( const int k, const double th_a, const double th_b, const int hilo, double *const err )  const 
//      Get integral on rho-monotonic interval  ( tha, thb )  where  tha  is nearer  theta0.
//  Integrand negligible for  abs(mu)/Omega > 7,  which has one or no maxima on  rho > z/w  and 
//  one or no maxima on  rho < z/w.
{
	if (err!=0)  *err = 0.;
	if( fabs( th_a - th_b ) < zero_eq )  return  0.;		

	double  tha =th_a,  thb =th_b,  ra= rho(tha,k),  rb= rho(thb,k),  rZ= z/w,  arga= 0., argb= 0.;

	bool  zero= false;
	if( (ra-rZ)*(rb-rZ) < 0 )  zero= true;		// check for zero of 'mu'


//Rcout << endl << "NDab:  k ra rb dif rZ  " << k << " " << ra << " " << rb << " " << ra-rb << " " << rZ << endl;

	if ( ra > rb )  {	// case of 'rho' decreasing as theta moves away from theta0

		if (ra > rZ) {
			if( ra >= 1 )  arga = Inf;  else  arga = (w-z*ra)/sqrt(1-ra*ra);
			if (rb < rZ)  argb = sqrt( w*w-z*z );  else  argb = (w-z*rb)/sqrt(1-rb*rb);
		}

	}  else  {		//  rho increasing

		if (ra < rZ) {
			arga = (w-z*ra)/sqrt(1-ra*ra);
			if (rb > rZ)  argb = sqrt( w*w-z*z );  else  argb = (w-z*rb)/sqrt(1-rb*rb);
		}
	}

	double  pr =  Rf_pnorm5(arga ,0,1,1,0)  -  Rf_pnorm5(argb ,0,1,1,0) ;
//Rcout << "arga argb 1-B pr  " << arga << " " << argb << " " << 1-B[k] << " " << pr << endl;



	const double  aa= amu_by_Omega(tha,k),  ab= amu_by_Omega(thb,k);
//Rcout << "thb tha zero aa ab pr  " << thb << " " << tha << " " << zero << " " << aa << " " << ab << " " << pr << endl;

	if( !zero  &&  aa > 6.5  &&  ab > 6.5)  return pr;

	double thZ= NaN;
	if( zero )  thZ= rho_inv( rZ, k, hilo );
	if( fabs(thZ-tha) < zero_eq || fabs(thZ-thb) < zero_eq || isinf(thZ) )  zero= false;

//Rcout << "thZ rho(thZ) z-w*ro  " << thZ << " " << rho(thZ,k) << " " << z-w*rho(thZ,k) << endl;


//  'Rdqag' parameters
	int  inf_flag = -1,  neval =0,  ier =0,  limit =100,  lenw = 4*limit,  last =0,  iwork[limit];
	double  epsabs = acc_sl_abs/2/n,  epsrel = acc_sl_rel/2,  result =0,  abserr =0,  work[lenw];

	const void *  exc[2] = {  this,  &k  };
	void **  ex = const_cast< void** >( exc );

	double  error =0;
	if( zero )  {

		if( aa > 7.5 )  tha= bisect( tha, thZ, &Clmbr::amu_by_Omega, k, 7, inc_x);			
		if( ab > 7.5  &&  !isinf(thb) )  thb= bisect( thZ, thb, &Clmbr::amu_by_Omega, k, 7, inc_x);			
//Rcout << "tha thZ thb  " << tha << " " << thZ << " " << thb << endl;

		Rdqags( igeo, ex, &tha, &thZ, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
			&limit, &lenw, &last, iwork, work );

		pr += fabs(result);
		error += abserr;
		if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 


//Rcout << "tha thZ pr  " << tha << " " << thZ << " " << pr << endl;

		if( isinf(thb) )  
			Rdqagi( igeo, ex, &thZ, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );
		else
			Rdqags( igeo, ex, &thZ, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );

		pr += fabs(result);
		error += abserr;
		if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 

//Rcout << "thZ thb pr  " << thZ << " " << thb << " " << pr << endl;

	}  else  {

		if(aa>7.5 && ab<6.5)  {
			if( !isinf(thb) )
				tha= bisect( tha, thb, &Clmbr::amu_by_Omega, k, 7, inc_x);
			else  {
				double th= min( tha, -1. );
				while( amu_by_Omega(th,k) > 6.8 )  th *= 2;
				tha= bisect( tha, th, &Clmbr::amu_by_Omega, k, 7, inc_x);
			}
		}
		if(aa<6.5 && ab>7.5 && !isinf(thb) )  thb= bisect( tha, thb, &Clmbr::amu_by_Omega, k, 7, inc_x);

//Rcout << "tha thb pr1  " << tha << " " << thb << " " << pr << endl;

		if( isinf(thb) )
			Rdqagi( igeo, ex, &tha, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );
		else
			Rdqags( igeo, ex, &tha, &thb, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
				&limit, &lenw, &last, iwork, work );

		pr += fabs(result);
		error += abserr;
		if( ier > 0  &&  ier!=5 )  Rf_warning( _("integration flag") ); 

//Rcout << "tha thb pr2  " << tha << " " << thb << " " << pr << endl;
	}

//Rcout << "NDab: k pr   " << k << " " << pr << endl << endl;

	if (err!=0)  *err = error;
	return  min( pr, 1. );
}




void  igeo(double *const x, const int n, void *const ex)
// integrand for Rdqags and Rdqagi
{
	const Clmbr  **const  ppObj  = static_cast< const Clmbr **const > ( ex );
	const Clmbr  *const  pObj = ppObj[0];

	const int  **const  ppk  = static_cast< const int **const > ( ex );
	const int  *const  pk = ppk[1];

	if( (*pObj).variance_unknown ) 
		for (int i =0; i < n; i++ )  x[i] = (*pObj).Emupr( x[i], *pk );
	else  
		for (int i =0; i < n; i++ )  x[i] = (*pObj).Emupr_vk( x[i], *pk );

	return;
}



