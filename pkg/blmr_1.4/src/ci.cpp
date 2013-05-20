//

#include "blmr.h"



int Cblmr::ci(const METHOD met, const double incr, const bool output, double *const bounds)
// check whether {theta  such that  sig. level  > SL}  is contiguous
// return number of contiguous segments and boundaries of each segment
// 'incr' specifies increments to cover in grid search for method GEO
{
	int numr = 0;

	double *const bds= new (nothrow) double[2*n_d];
	if(bds==NULL) {
		Rcout << _("message: ") << 12 << endl;
		stop( _("memory allocation failed") );
	}


	if(trivial) {

		numr= 1.;
		const double thmle= mle(false);
		if( isnan(thmle) )  { bds[0]= -inf; bds[1]= inf; }  
			else  { if( thmle==x_d[0] )  { bds[0]= -inf; bds[1]= thmle; }		// single line in Model 2
				else  bds[0] = bds[1] = thmle; }

	}  else  {

		if (met==GEO || met==GEO2)  numr = ci_geo(met,incr,bds);
		if (met==AF || met==AF2)  numr = ci_af(met,bds);

	}


	if (output)  {
		Rcout << _( "confidence level,  confidence interval for theta,  and  method:" ) << endl;
		Rcout << "  " << 1-SL << "    ";
		if( model_in != -2 )  {
			for (int i=0;i<2*numr;i+=2) {
				Rcout << "[ ";
				if(bds[i]==-inf) Rcout << "-inf"; else Rcout << bds[i];
				Rcout  << ", ";
				if(bds[i+1]==inf) Rcout << "inf"; else Rcout << bds[i+1];
				Rcout << " ]";
				if (i+2<2*numr) Rcout << ",  ";
			}
		}  else  {
			for (int i=2*numr-2;i>=0;i-=2) {
				Rcout << "[ ";
				if(bds[i+1]==inf) Rcout << "-inf"; else Rcout << -bds[i+1];
				Rcout  << ", ";
				if(bds[i]==-inf) Rcout << "inf"; else Rcout << -bds[i];
				Rcout << " ]";
				if (i-2>=0) Rcout << ",  ";
			}

		}
		Rcout << "    ";
		if (met==GEO) Rcout << "CLR";
		if (met==AF) Rcout << "AF";

		Rcout << endl << endl;
	}


	if(bounds != 0)  for (int i=0;i<2*numr;i+=2)  { bounds[i] = bds[i];  bounds[i+1] = bds[i+1]; }

	delete[] bds;

	return numr;
}




int Cblmr::ci_geo( const METHOD met, const double incr, double *const bds )
// Using Knowles, Siegmund and Zhang's geometric formula to calculate significance level.
// In Model=M1, treat regions before x[0] and after x[n-1] as two seperate regions.
// Conditional SL(th,mle-alpha) is not constant on end-intervals, except on [x(n-1),x(n)] in M2.
// 'incr' specified so as to cover same theta values as in 'cr' routine.
{
	int  numi=0, ind=0;		// ind = indicator = {0 if sl_geo was below SL, 1 if above}
	double th, sl_th, thold;


// start with point  x_d[0]-1.
	th= x_d[0] - 1.;
	sl_th= sl(th,met,false);
	if (sl_th > SL) {
		bds[numi++] = -inf;
		ind= 1;
	}
	thold= th;

	if(Model==M1) {		// in M1, boundary is discontinuous at x(1)

		if(met==GEO) {
			th= x_d[1]; 
			sl_th= sl(th,met,false);
			if (sl_th > SL && ind==0) {
				bds[numi++] = x_d[0];
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = x_d[0];
				ind=0;
			}
			thold= th;
		}

		if(met==GEO2) {
			th = x_d[0] + acc_xb/2;
			sl_th = sl(th,met,false);
			if (sl_th > SL && ind==0) {
				bds[numi++] = x_d[0];
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = x_d[0];
				ind=0;
			}
			thold= th;
		}
	}



// get critical points

	double *const cpts= new (nothrow) double[n_d+1];
	if( cpts==NULL )  {
		Rcout << _("message: ") << 13 << endl;
		stop( _("memory allocation failed") );
	}

	const double thmle = mle(false);

	int k= 0, ncp= 0;
	if(Model==M1) k= 1;
	if(met==GEO2 && Model==M1) cpts[ncp++]= x_d[0]+acc_xb/2;
	if(met==GEO2 && Model==M2) cpts[ncp++]= x_d[0]-1;
	while(x_d[k]<thmle) cpts[ncp++]= x_d[k++];
	cpts[ncp++]= thmle;
	if(thmle==x_d[k]) k++;
	while( k < n_d-1 ) cpts[ncp++]= x_d[k++];
	if(met==GEO2 && Model==M1) cpts[ncp++]= x_d[n_d-1]-acc_xb/2;



// grid search
	for (k = 0; k < ncp - 1; k++) {
		th= cpts[k];
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
		double inc= incr;
		while( (cpts[k+1]-cpts[k])/inc < subints + 1 )  inc /= 2.;
		double  fth= floor(th);
		while( fth < cpts[k] + acc_xb )  fth += inc;
		for (th=fth;th<cpts[k+1];th+=inc) {		// covers the same points as in 'cr' routine
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
	}


// check boundary of final end-interval
	th = cpts[k];
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
		ind=1;
	}
	if (sl_th < SL && ind==1) {
		bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
		ind=0;
	}


// check region after x_d[n_d-1]
	th = x_d[n_d-1]+1.;
	sl_th = sl(th,met,false);
	if (sl_th < SL && ind==1) bds[numi++] = x_d[n_d-1];
	if (sl_th > SL && ind==0) {
		bds[numi++] = x_d[n_d-1];
		bds[numi++] = inf;
	}
	if (sl_th > SL && ind==1) bds[numi++] = inf;


	delete[] cpts;
	return numi/2;
}





int Cblmr::ci_af( const METHOD met, double *const bds )
// using AF to calculate significance level
{
	double th, sl_th;
	int k, numi=0, ind=0;	// ind = indicator = 0 if sl was below SL, 1 above

// if Model=M1, treat regions before x_d[0] and after x_d[n_d-1] as two seperate regions

// check region before x_d[0]
	th = x_d[0] - 1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL) {
		bds[numi++] = -inf;
		ind=1;
	}

// check first end interval
	if (Model==M1) {
		th = x_d[1]; 
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = x_d[0];
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = x_d[0];
			ind=0;
		}
	}

//check critical points
// (gamma*sy)^2 is monotonic between each of  x_d[k], thzero, thk, x_d[k+1]

	double  thold = x_d[k1], yf = *pqy*q_f(x_d[k1],k1);

	for (k=k1;k<n_d-2;k++) 
	{
		const double  y1 = *pqy*pq1[k+1],  yx = yf + y1*x_d[k];	
		yf = yx - y1*x_d[k+1];

		const double ya =yx*qx1[k+1]-y1*qxx[k+1], yb =yx*q11[k+1]-y1*qx1[k+1], thk =ya/yb, thzero =yx/y1;

		const double th1 = min(thk,thzero), th2 = max(thk,thzero);

		if (x_d[k] < th1 && th1 < x_d[k+1]) {
			th = th1;
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

		if (x_d[k] < th2 && th2 < x_d[k+1]) {
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

		th = x_d[k+1];
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
	th = (x_d[n_d-1]+x_d[n_d-2])/2.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = bisect_sl(thold,x_d[n_d-2],met,-acc_xb);
		ind=1;
	}
	if (sl_th < SL && ind==1) {
		bds[numi++] = bisect_sl(thold,x_d[n_d-2],met,-acc_xb);
		ind=0;
	}

// check region after x_d[n_d-1]
	th = x_d[n_d-1]+1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = x_d[n_d-1];
		bds[numi++] = inf;
	}
	if (sl_th < SL && ind==1) bds[numi++] = x_d[n_d-1];
	if (sl_th > SL && ind==1) bds[numi++] = inf;


	return numi/2;
}

