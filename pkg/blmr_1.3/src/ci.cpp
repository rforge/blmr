//

#include "blmr.h"



int Cblmr::ci(METHOD met, bool output, double *bounds)
// check whether {theta  such that  sig. level  > SL}  is contiguous
// return number of contiguous segments, and start and finish boundary points for each segment
{
	int k =0;
	double *bds;
	bds = new double[2*ns];

	if (met==GEO || met==GEO2) k = ci_geo(met,bds);
	if (met==AF || met==AF2) k = ci_af(met,bds);

	if (output) {
		Rcout << _( "confidence level,  confidence interval for theta,  and  method:" ) << endl;
		Rcout << "  " << 1-SL << "    ";
		for (int i=0;i<2*k;i+=2) {
			Rcout << "[ ";
			if(bds[i]==-inf) Rcout << "-inf"; else Rcout << bds[i];
			Rcout  << ", ";
			if(bds[i+1]==inf) Rcout << "inf"; else Rcout << bds[i+1];
			Rcout << " ]";
			if (i+2<2*k) Rcout << ",  ";
		}
		Rcout << "    ";
		if (met==GEO) Rcout << "CLR";
		if (met==AF) Rcout << "AF";

		Rcout << endl << endl;
	}

	if (bounds != 0)  for (int i=0;i<2*k;i+=2)  {bounds[i] = bds[i]; bounds[i+1] = bds[i+1];}
	delete[] bds;
	return k;
}




int Cblmr::ci_geo( METHOD met, double *bds )
// using Knowles, Siegmund and Zhang's geometric formula to calculate significance level
{
	double th, sl_th;
	int k, numi=0, ind=0;	// ind = indicator = 0 if sl_geo was below SL, 1 above

// if model=M1, treat regions before x[0] and after x[n-1] as two seperate regions

// check region before Xs[0]
	th = Xs[0] - 1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL) {
		bds[numi++] =-inf;
		ind=1;
	}

// check first end interval
	if (model==M1) {
		th = (Xs[1]+Xs[0])/2.; 
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = Xs[0];
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = Xs[0];
			ind=0;
		}
	}


// two grid searches, before and after thmle

	double thmle = mle(false);
	double thold = Xs[k1];

	int kmle = 0;
	bool thmle_eq_datapt = false;
	while (Xs[kmle+1] < thmle + zero_eq) kmle++;
	if ( fabs(Xs[kmle] - thmle) < zero_eq) {
		thmle = Xs[kmle];
		thmle_eq_datapt = true;
	}

// get critical points
	int num_cpts;
	double *cpts;
	cpts = new double[ns+1];
	if (thmle_eq_datapt) {
		for (k=k1;k<ns-1;k++) cpts[k-k1] = Xs[k]; 
		num_cpts = ns-1-k1;
	} else {
		for (k=k1;k<=kmle;k++) cpts[k-k1] = Xs[k];
		cpts[k-k1] = thmle;
		for (k=kmle+1;k<ns-1;k++) cpts[k-k1+1] = Xs[k];
		num_cpts = ns-k1;
	}

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

// check boundary of final end interval
	th = Xs[ns-2];
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
	th = (Xs[ns-1]+Xs[ns-2])/2.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = Xs[ns-2];
		ind=1;
	}
	if (sl_th < SL && ind==1) {
		bds[numi++] = Xs[ns-2];
		ind=0;
	}


// check region after Xs[ns-1]
	th = Xs[ns-1]+1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = Xs[ns-1];
		bds[numi++] = inf;
	}
	if (sl_th < SL && ind==1) bds[numi++] = Xs[ns-1];
	if (sl_th > SL && ind==1) bds[numi++] = inf;


	delete[] cpts;
	return numi/2;
}





int Cblmr::ci_af( METHOD met, double *bds )
// using AF to calculate significance level
{
	double th, sl_th;
	int k, numi=0, ind=0;	// ind = indicator = 0 if sl was below SL, 1 above

// if model=M1, treat regions before Xs[0] and after Xs[ns-1] as two seperate regions

// check region before Xs[0]
	th = Xs[0] - 1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL) {
		bds[numi++] = -inf;
		ind=1;
	}

// check first end interval
	if (model==M1) {
		th = (Xs[1]+Xs[0])/2.; 
		sl_th = sl(th,met,false);
		if (sl_th > SL && ind==0) {
			bds[numi++] = Xs[0];
			ind=1;
		}
		if (sl_th < SL && ind==1) {
			bds[numi++] = Xs[0];
			ind=0;
		}
	}

//check critical points
// (gamma*sy)^2 is monotonic between each of  Xs[k], thzero, thk, Xs[k+1]

	double  thold = Xs[k1], yf = *pqy*q_f(Xs[k1],k1);

	for (k=k1;k<ns-2;k++) 
	{
		double  y1 = *pqy*pq1[k],  yx = yf + y1*Xs[k];	yf = yx - y1*Xs[k+1];

		double ya =yx*qx1[k]-y1*qxx[k], yb =yx*q11[k]-y1*qx1[k], thk =ya/yb, thzero =yx/y1;

		double th1 = min(thk,thzero), th2 = max(thk,thzero);

		if (Xs[k] < th1 && th1 < Xs[k+1]) {
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

		if (Xs[k] < th2 && th2 < Xs[k+1]) {
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

		th = Xs[k+1];
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
	th = (Xs[ns-1]+Xs[ns-2])/2.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = bisect_sl(thold,Xs[ns-2],met,-acc_xb);
		ind=1;
	}
	if (sl_th < SL && ind==1) {
		bds[numi++] = bisect_sl(thold,Xs[ns-2],met,-acc_xb);
		ind=0;
	}

// check region after Xs[ns-1]
	th = Xs[ns-1]+1.;
	sl_th = sl(th,met,false);
	if (sl_th > SL && ind==0) {
		bds[numi++] = Xs[ns-1];
		bds[numi++] = inf;
	}
	if (sl_th < SL && ind==1) bds[numi++] = Xs[ns-1];
	if (sl_th > SL && ind==1) bds[numi++] = inf;


	return numi/2;
}
