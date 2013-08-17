//

#include "lmbr.h"



int Clmbr::ci(const METHOD met, const double incr, const bool output, double *const bounds)
// check whether {theta  such that  sig. level  > SL}  is contiguous
// return number of contiguous segments and boundaries of each segment
// 'incr' specifies increments to cover in grid search for method 'GEO'
{
	int numr = 0;

	double *const bds= new (nothrow) double[2*ns];
	if(bds==NULL) {
		Rcout << _("message: ") << 11 << endl;
		stop( _("memory allocation failed") );
	}


	if(trivial) {

		numr= 1.;
		const double thmle= mle(false);
		if( isnan(thmle)  &&  !isinf(thmle) )  { bds[0]= -Inf; bds[1]= Inf; }  
			else  { if( Model==M2 && thmle==xs[0] )  { bds[0]= -Inf; bds[1]= thmle; }		// single line in Model 2
				else  bds[0] = bds[1] = thmle; }

	}  else  {

		if (met==GEO || met==GEO2)  numr = ci_geo(met,incr,bds);
		if (met==AF || met==AF2)  numr = ci_af(met,bds);

	}


	if (output)  {
		Rcout << _( "Confidence interval for changepoint,  confidence level,  and  method:" ) << endl;
		Rcout << "    ";
		if( model_in > 0 )  {
			for (int i=0;i<2*numr;i+=2) {
				Rcout << "[ ";
				if(bds[i]==-Inf) Rcout << "-Inf"; else Rcout << bds[i];
				Rcout  << ", ";
				if(bds[i+1]==Inf) Rcout << "Inf"; else Rcout << bds[i+1];
				Rcout << " ]";
				if (i+2<2*numr) Rcout << ",  ";
			}
		}  else  {
			for (int i=2*numr-2;i>=0;i-=2) {
				Rcout << "[ ";
				if(bds[i+1]==Inf) Rcout << "-Inf"; else Rcout << -bds[i+1];
				Rcout  << ", ";
				if(bds[i]==-Inf) Rcout << "Inf"; else Rcout << -bds[i];
				Rcout << " ]";
				if (i-2>=0) Rcout << ",  ";
			}

		}
		Rcout << "      " << 1-SL << "      ";
		if (met==GEO) Rcout << "CLR";
		if (met==AF) Rcout << "AF";

		Rcout << endl << endl;
	}


	if(bounds != 0)  for (int i=0;i<2*numr;i+=2)  { bounds[i] = bds[i];  bounds[i+1] = bds[i+1]; }

	delete[] bds;

	return numr;
}




int Clmbr::ci_geo( const METHOD met, const double incr, double *const bds )
// Using Knowles, Siegmund and Zhang's geometric formula to calculate significance level.
// In Model=M1, treat regions before x[0] and after x[n-1] as two seperate regions.
// Conditional SL(th,mle-alpha) is not constant on end-intervals.
// 'incr' specified so as to cover same theta values as in 'cr' routine.
{
	int  numi=0, ind=0;		// ind = indicator = {0 if sl_geo was below SL, 1 if above}
	double th, sl_th =2, thold;
	Rcpp::Function Rflush("flush.console");

	if( Model==M3 )  {

		const double  sl_inf= sl(-Inf,met,false);
		th= min( -1., xs[0] );
		int it=0;
		while( fabs( sl_th - sl_inf ) > acc_sl_abs  &&  it < 18 )  
			th *= 2,  sl_th= sl(th,met,false),  it++;
		if (sl_th > SL)  bds[numi++] = -Inf,  ind= 1;
		thold= th;
		const double  thi = th,  inc = (xs[0]-thi)/(subints+0.5);
//Rcout << "thi inc  " << thi << " " << inc << endl;

		for (th=thi;th<xs[0];th+=inc) {
			sl_th = sl(th,met,false);
//Rcout << "th sl ind thold  " << th << " " << sl_th << " " << ind << " " << thold << endl;
			if (sl_th > SL && ind==0) {
				bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
				ind=0;
			}
			thold = th;
		}


	}  else  {
	
// start with point  xs[0]-1.
		th= xs[0] - 1.;
		sl_th= sl(th,met,false);
		if (sl_th > SL) {
			bds[numi++] = -Inf;
			ind= 1;
		}
		thold= th;

		if(Model==M1) {		// in M1, boundary is discontinuous at x[0]

			if(met==GEO) {
				th= xs[1]; 
				sl_th= sl(th,met,false);
				if (sl_th > SL && ind==0) {
					bds[numi++] = xs[0];
					ind=1;
				}
				if (sl_th < SL && ind==1) {
					bds[numi++] = xs[0];
					ind=0;
				}
				thold= th;
			}

			if(met==GEO2) {
				th = xs[0] + acc_xb/2;
				sl_th = sl(th,met,false);
				if (sl_th > SL && ind==0) {
					bds[numi++] = xs[0];
					ind=1;
				}
				if (sl_th < SL && ind==1) {
					bds[numi++] = xs[0];
					ind=0;
				}
				thold= th;
			}
		}

	}


// get critical points

	double *const cpts= new (nothrow) double[ns+2];
	if( cpts==NULL )  {
		Rcout << _("message: ") << 12 << endl;
		stop( _("memory allocation failed") );
	}

	const double thmle = mle(false);

	double  thm = thmle;
	if(Model==M3)  thm = max( thmle, thold);

	int k= 0, ncp= 0;
	if(Model==M1) k= 1;
	if(met==GEO2 && Model==M1) cpts[ncp++]= xs[0]+acc_xb/2;
	if(met==GEO2 && Model==M2) cpts[ncp++]= xs[0]-1;
	while(xs[k]<thm) cpts[ncp++]= xs[k++];
	cpts[ncp++]= thm;
	if(thmle==xs[k]) k++;
	while( k < ns-1 ) cpts[ncp++]= xs[k++];
	if(met==GEO2) cpts[ncp++]= xs[ns-1]-acc_xb/2;

//Rcout << "ncp " << ncp << endl;
//for(k = 0; k < ncp; k++) Rcout << "k cpt  " << k << " " << cpts[k] << endl;
	bool msg= false;
	double  tstart= time(NULL);

// grid search
	for (k = 0; k < ncp - 1; k++) {
		th= cpts[k];
			sl_th = sl(th,met,false);
//Rcout << "th sl ind thold  " << th << " " << sl_th << " " << ind << " " << thold << endl;
			if (sl_th > SL && ind==0) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind= 1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind= 0;
			}
			thold = th;
		double inc= incr;
		if(Model==M3 && k==0) inc= (cpts[k+1]-cpts[k])/(subints+0.5);
		while( (cpts[k+1]-cpts[k])/inc < subints + 1 )  inc /= 2.;
		double  fth= floor(th);
		while( fth < cpts[k] + acc_xb )  fth += inc;
//Rcout << "fth inc  " << fth << " " << inc << endl;
		for (th=fth;th<cpts[k+1];th+=inc) {		// covers the same points as in 'cr' routine
			sl_th = sl(th,met,false);
//Rcout << "th sl ind thold  " << th << " " << sl_th << " " << ind << " " << thold << endl;
			if (sl_th > SL && ind==0) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind= 1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind= 0;
			}
			thold = th;

			double  tfinish = time( NULL ),  elapsed = tfinish - tstart;
			if( elapsed > 10 ) {
				const double  progress = floor( 100*(th-cpts[0])/(cpts[ncp-1]-cpts[0]) );
				if(!msg) { msg=true; if(met==GEO) Rcout << "   progress:   "; }
				Rcout << progress << "%...   ";  Rflush();
				tstart= tfinish;
			}
		}

	}
	if(msg && met==GEO) Rcout << endl << endl;

// check boundary of final end-interval
	th = cpts[k];
	sl_th = sl( th, met, false );
//Rcout << "th sl ind thold  " << th << " " << sl_th << " " << ind << " " << thold << endl;
	if (sl_th > SL && ind==0) {
		bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
		ind= 1;
	}
	if (sl_th < SL && ind==1) {
		bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
		ind= 0;
	}


// check region after xs[ns-1]
	th = xs[ns-1]+1.;
	sl_th = sl(th,met,false);
//Rcout << "th sl ind thold  " << th << " " << sl_th << " " << ind << " " << thold << endl;
	if (sl_th < SL && ind==1) bds[numi++] = xs[ns-1];
	if (sl_th > SL && ind==0) {
		bds[numi++] = xs[ns-1];
		bds[numi++] = Inf;
	}
	if (sl_th > SL && ind==1) bds[numi++] = Inf;


	delete[] cpts;
	return numi/2;
}





int Clmbr::ci_af( const METHOD met, double *const bds )
// using AF to calculate significance level
{
	double  th,  sl_th,  thold;
	int  k,  numi=0,  ind=0;	// ind = indicator = 0 if sl was below SL, 1 if above


	if( Model==M3 )  {

		const double  sl_inf= sl(-Inf,met,false);

		const double  yx = *pqy*pqx[0],  y1 = *pqy*pq1[0],  thzero =yx/y1; 
		const double  ya =yx*qx1[0]-y1*qxx[0],  yb =yx*q11[0]-y1*qx1[0],  thk =ya/yb;
		const double  th1 = min(thk,thzero),  th2 = max(thk,thzero);


		th= min( -1., xs[0] );
		th= min( th, th1 );
		th *= 2,  sl_th= sl( th, met, false );
		int it =0;
		while( fabs( sl_th - sl_inf ) > acc_sl_abs  &&  it < 18 )  
			th *= 2,  sl_th= sl(th,met,false),  it++;
		if (sl_th > SL)  bds[numi++] = -Inf,  ind=1;
		thold= th;

		if ( th1 < xs[0]) {
			th = th1;
			sl_th = sl( th, met, false );
			if (sl_th > SL && ind==0) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind=0;
			}
			thold = th;
		}

		if ( th2 < xs[0]) {
			th = th2;
			sl_th = sl( th, met, false );
			if (sl_th > SL && ind==0) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind=0;
			}
			thold = th;
		}

		th = xs[0];
		sl_th = sl( th, met, false );
		if (sl_th > SL && ind==0) {
			bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
			ind=0;
		}
		thold = th;


	}  else  {
	

// if Model=M1, treat regions before xs[0] and after xs[ns-1] as two seperate regions

// check region before xs[0]
		th = xs[0] - 1.;
		sl_th = sl( th, met, false );
		if (sl_th > SL) {
			bds[numi++] = -Inf;
			ind=1;
		}

// check first end interval
		if (Model==M1) {
			th = xs[1]; 
			sl_th = sl( th, met, false );
			if (sl_th > SL && ind==0) {
				bds[numi++] = xs[0];
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = xs[0];
				ind=0;
			}
		}

	}


//check critical points
// (gamma*sy)^2 is monotonic between each of  xs[k],  'thzero' where gamma*sy=0,  'thk' where d(gamma*sy)=0,  and  xs[k+1]

	int  ki = max(k1,0);
	double  yf = *pqy*q_f(xs[ki],ki);
	thold = xs[ki];

	for (k=ki;k<ns-2;k++) 
	{
		const double  y1 = *pqy*pq1[k+1],  yx = yf + y1*xs[k];	
		yf = yx - y1*xs[k+1];

		const double ya =yx*qx1[k+1]-y1*qxx[k+1], yb =yx*q11[k+1]-y1*qx1[k+1], thk =ya/yb, thzero =yx/y1;

		const double th1 = min(thk,thzero), th2 = max(thk,thzero);

		if (xs[k] < th1 && th1 < xs[k+1]) {
			th = th1;
			sl_th = sl( th, met, false );
			if (sl_th > SL && ind==0) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = bisect_sl( thold, th, met, -acc_xb );
				ind=0;
			}
			thold = th;
		}

		if (xs[k] < th2 && th2 < xs[k+1]) {
			th = th2;
			sl_th = sl(th,met,false);
			if (sl_th > SL && ind==0) {
				bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
				ind=0;
			}
			thold = th;
		}

		th = xs[k+1];
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
			ind=0;
		}
		thold = th;

	}



// check final end interval
	th = (xs[ns-1]+xs[ns-2])/2.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = bisect_sl(thold,xs[ns-2],met,-acc_xb);
		ind=1;
	}
	if (sl_th < SL && ind==1) {
		bds[numi++] = bisect_sl(thold,xs[ns-2],met,-acc_xb);
		ind=0;
	}

// check region after xs[ns-1]
	th = xs[ns-1]+1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = xs[ns-1];
		bds[numi++] = Inf;
	}
	if (sl_th < SL && ind==1) bds[numi++] = xs[ns-1];
	if (sl_th > SL && ind==1) bds[numi++] = Inf;


	return numi/2;
}

