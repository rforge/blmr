//

#include "blmr.h"



int Cblmr::ci(METHOD met, bool output, double *bounds)
// check whether {theta  such that  sig. level  > SL}  is contiguous
// return number of contiguous segments, and start and finish boundary points for each segment
{
	int numr =0;
	double *bds;
	bds = new double[2*ns];

	if (met==GEO || met==GEO2) numr = ci_geo(met,bds);
	if (met==AF || met==AF2) numr = ci_af(met,bds);

	if (output) {
		Rcout << _( "confidence level,  confidence interval for theta,  and  method:" ) << endl;
		Rcout << "  " << 1-SL << "    ";
		for (int i=0;i<2*numr;i+=2) {
			Rcout << "[ ";
			if(bds[i]==-inf) Rcout << "-inf"; else Rcout << bds[i];
			Rcout  << ", ";
			if(bds[i+1]==inf) Rcout << "inf"; else Rcout << bds[i+1];
			Rcout << " ]";
			if (i+2<2*numr) Rcout << ",  ";
		}
		Rcout << "    ";
		if (met==GEO) Rcout << "CLR";
		if (met==AF) Rcout << "AF";

		Rcout << endl << endl;
	}

	if(bounds != 0)  for (int i=0;i<2*numr;i+=2)  {bounds[i] = bds[i]; bounds[i+1] = bds[i+1];}
	delete[] bds;
	return numr;
}




int Cblmr::ci_geo( METHOD met, double *bds )
// using Knowles, Siegmund and Zhang's geometric formula to calculate significance level
{
	double th, sl_th, thold;
	int k, numi=0, ind=0;	// ind = indicator = {0 if sl_geo was below SL, 1 if above}

// if Model=M1, treat regions before x[0] and after x[n-1] as two seperate regions


// start with point  xs[0]-1.
	th = xs[0] - 1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL) {
		bds[numi++] = -inf;
		ind=1;
	}
	thold=xs[0] - 1.;


// when met==GEO2, scan first interval in M1 or scan below xs[0] in M2, 
// because 2-parameter CLR SL's along ridge are not constant there

	if(met==GEO2) {

		if(Model==M2) {
// scan region before xs[0]
			double inc = 1./(subints+0.5);
			for (th=xs[0] - 1.+inc;th<xs[0];th+=inc) {
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

		if(Model==M1) {		// SL discontinuous at xs[0]

// check boundary of the region below xs[0]
			th = xs[0]*(1+rel_print_eps);
			sl_th = sl(th,met,false);
			if (sl_th > SL && ind==0) {
				bds[numi++] = xs[0];
				ind=1;
			}
			if (sl_th < SL && ind==1) {
				bds[numi++] = xs[0];
				ind=0;
			}
			thold = th;

// scan first end interval
			double inc = (xs[1]-xs[0])/(subints+0.5);
			for (th=xs[0]+inc;th<xs[1];th+=inc) {
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
	}



// check first end interval in M1
	if (met==GEO && Model==M1) {
		th = (xs[1]+xs[0])/2.; 
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = xs[0];
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = xs[0];
			ind=0;
		}
	}


	thold = xs[k1];		// k1 = { 1 in M1, 0 in M2 }


// two grid searches, before and after 'thmle'

// get critical points
	double thmle = mle(false);
	int kmle = 0;
	while ( xs[ kmle + 1 ] < thmle + zero_eq )  kmle++;
	bool thmle_eq_datapt = false;
	if ( fabs(xs[kmle] - thmle) < zero_eq) {
		thmle = xs[kmle];
		thmle_eq_datapt = true;
	}

	int num_cpts;
	double *cpts;
	cpts = new double[ns+1];
	if (thmle_eq_datapt) {
		for (k=k1;k<ns-1;k++) cpts[k-k1] = xs[k]; 
		num_cpts = ns-1-k1;
	} else {
		for (k=k1;k<kmle+1;k++) cpts[k-k1] = xs[k];
		cpts[k-k1] = thmle;
		for (k=kmle+1;k<ns-1;k++) cpts[k-k1+1] = xs[k];
		num_cpts = ns-k1;
	}

// grid search
	for (k = 0; k < num_cpts - 1; k++) {
		double inc = (cpts[k+1]-cpts[k])/(subints+0.5);
		for (th=cpts[k];th<cpts[k+1];th+=inc) {
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


	if(met==GEO || (met==GEO2 && Model==M2)) {

// check boundary of final end interval
		th = xs[ns-2];
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = bisect_sl(thold,th,met,-acc_xb);
			ind=0;
		}

// check final end interval
		th = (xs[ns-1]+xs[ns-2])/2.;
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = xs[ns-2];
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = xs[ns-2];
			ind=0;
		}

	} else {

// scan final end interval
		double inc = (xs[ns-1]-xs[ns-2])/(subints+0.5);
		for (th=xs[ns-2];th<xs[ns-1];th+=inc) {
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


// check region after xs[ns-1]
	th = xs[ns-1]+1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = xs[ns-1];
		bds[numi++] = inf;
	}
	if (sl_th < SL && ind==1) bds[numi++] = xs[ns-1];
	if (sl_th > SL && ind==1) bds[numi++] = inf;


	delete[] cpts;
	return numi/2;
}





int Cblmr::ci_af( METHOD met, double *bds )
// using AF to calculate significance level
{
	double th, sl_th;
	int k, numi=0, ind=0;	// ind = indicator = 0 if sl was below SL, 1 above

// if Model=M1, treat regions before xs[0] and after xs[ns-1] as two seperate regions

// check region before xs[0]
	th = xs[0] - 1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL) {
		bds[numi++] = -inf;
		ind=1;
	}

// check first end interval
	if (Model==M1) {
		th = (xs[1]+xs[0])/2.; 
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = xs[0];
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = xs[0];
			ind=0;
		}
	}

//check critical points
// (gamma*sy)^2 is monotonic between each of  xs[k], thzero, thk, xs[k+1]

	double  thold = xs[k1], yf = *pqy*q_f(xs[k1],k1);

	for (k=k1;k<ns-2;k++) 
	{
		double  y1 = *pqy*pq1[k],  yx = yf + y1*xs[k];	yf = yx - y1*xs[k+1];

		double ya =yx*qx1[k]-y1*qxx[k], yb =yx*q11[k]-y1*qx1[k], thk =ya/yb, thzero =yx/y1;

		double th1 = min(thk,thzero), th2 = max(thk,thzero);

		if (xs[k] < th1 && th1 < xs[k+1]) {
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
		bds[numi++] = inf;
	}
	if (sl_th < SL && ind==1) bds[numi++] = xs[ns-1];
	if (sl_th > SL && ind==1) bds[numi++] = inf;


	return numi/2;
}
