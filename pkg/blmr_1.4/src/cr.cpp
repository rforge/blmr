//

#include "blmr.h"



int Cblmr::cr(METHOD met, double inc, bool output, double *bounds)
//
//      Checks whether { (theta,alpha)  such that  sig. level  > SL}  is contiguous.
// Returns N = number of rows for "bounds[N][3]" array of (th, low_alpha, high_alpha) values.
// If (met==GEO) using Knowles, Siegmund and Zhang's geometric formula to calculate  sig. levels;
// if (met==AF) using Approximate-F method.
//
//		This 'cr' routine determines the theta-boundaries by using 
//	the 'ci' subroutine to check the 2-parameter significance levels along 
//	the  ( theta0, MLE-alpha(theta0) )  ridge,  then finds the alpha-boundaries
//	at increments of theta0.  For a given theta0, in the broken line model, SL vs alpha
//  has a single maxima.  
//
{
	int i, numr, k, N=0;
	double  th =0;
	Rcpp::Function Rflush("flush.console");


// get theta-boundaries of confidence region(s)

	if(output) if(met==GEO) { Rcout << "   " << _("getting theta-boundaries...") << endl; Rflush();}

	METHOD  met2;
	if(met==GEO) met2 =GEO2; else met2 =AF2;
	double *tmp, *th_bds;
	tmp = new double[2*ns];
	numr = ci(met2,false,tmp);
	th_bds = new double[2*numr];
	for (i=0;i<2*numr;i++) th_bds[i] = tmp[i];
	delete[] tmp; tmp=NULL;

	if(th_bds[0]== -inf) th_bds[0]= xs[0]-1;
	if(th_bds[2*numr-1]== +inf) th_bds[2*numr-1]= xs[ns-1]+1;

// store boundary values in an  N x 3  array	
	double width=0;						
	for (i=0;i<numr;i++) width += th_bds[2*i+1] - th_bds[2*i];
	int Nmax = width/inc + 1 + ns + 2*numr + 2;
	double *bds;
	bds = new double[Nmax*3];


	if(output) if(met==GEO) { Rcout << "   " << _("getting alpha-boundaries..."); Rflush();}
	double min_th = max(th_bds[0],xs[1]), max_th = min(th_bds[2*numr-1],xs[ns-2]), mid_th = (min_th + max_th)/2.;
	bool progress_msg1 = false, progress_msg2 = false, progress_msg3 = false;


	for (i=0;i<numr;i++) {

		double tha = th_bds[2*i], thb = th_bds[2*i+1];

		if ( fabs(thb-tha) < zero_eq) continue;

// in M1:
// confidence region starts at "tha"  with a vertical line if tha=x(1) or x(n), 
// an alpha-MLE point if  x(1) < tha < x(n),  or an open-end if tha < x(1)

// in M2:
// confidence region starts at "tha"  with a vertical line if tha=x(n), 
// an alpha-MLE point if x(1)-1< tha <x(n), or an open-end if tha = x(1)-1

		double high_a;
		if ( Model==M1 && fabs(tha-xs[0])<zero_eq ) {
			th = xs[0]*(1.+ rel_print_eps);
			high_a = ahigh(met,th);
			*(bds+N*3+0) =tha; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
			*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++; 
		}

		else  if (fabs(tha-xs[ns-1])<zero_eq) {
			tha = xs[ns-1];
			th = xs[ns-1]*(1.+rel_print_eps);
			high_a = ahigh(met,th);
			*(bds+N*3+0) =tha; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++; 
			*(bds+N*3+0) =th;  *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++; 
		}

		else  if ( (Model==M1 && tha > xs[0] ) || (Model==M2 && tha > xs[0]-1.) ) {
			high_a = ahigh(met,tha);
			*(bds+N*3+0) =tha; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
		}

		else {
			*(bds+N*3+0) =tha; *(bds+N*3+1) =a_sl(met,tha,-1); *(bds+N*3+2) =a_sl(met,tha,1); N++; 
		}

		double th1 = floor(tha);
		if ( (Model==M1 && tha < xs[0]-zero_eq) || (Model==M2 && tha==xs[0]-1.) ) 
			th1 =tha; 
		else while (th1<tha+zero_eq) th1 +=inc;
		if( fabs(th1)<zero_eq) th1=0.;

		int k_1=0; while (xs[k_1+1]<th1 && k_1<ns-1) k_1++;
		int kb=0; while (xs[kb+1]<thb && kb<ns-1) kb++;
			
		double th2;
		if (th1 < xs[0]-zero_eq) th2 = xs[0]; 
			else if(k_1==ns-1) th2=thb; 
				else th2 = min(xs[k_1+1],thb);

		for (k=k_1;k<=kb;k++) {

			for (th = th1;th < th2*(1.-rel_print_eps); th += inc) {

				if( fabs(th-tha)< acc_xb || fabs(th-thb) < acc_xb ) continue;
				if( fabs(th)<zero_eq) th=0.;

				*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++; 

				if(output) if(met==GEO) {
					if(!progress_msg1) if(th>(min_th + mid_th)/2) { Rcout << "  25%...";  Rflush();  progress_msg1 = true;  }
					if(!progress_msg2) if(th>mid_th) { Rcout << "  50%...";  Rflush();  progress_msg2 = true;  }
					if(!progress_msg3) if(th>(mid_th + max_th)/2) { Rcout << "  75%...";  Rflush();  progress_msg3 = true;  }
				}

			}

			th1 = th;
			if ( fabs(th1-th2) < rel_print_eps*(th1+th2)/2. )  th1 = th+inc;


// confidence region boundary is discontinuous at x(1) in M1, and at x(n) in M1 and M2

			if ( th2==xs[0] && th2!=thb) {
				*(bds+N*3+0) =th2; *(bds+N*3+1) =a_sl(met,th2,-1); *(bds+N*3+2) =a_sl(met,th2,1); N++; 
				if(Model==M1) {
				  th = th2*(1.+rel_print_eps);
				  *(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
				}

			} else  if ( th2==xs[ns-1] && th2!=thb) {
				th = th2*(1.- rel_print_eps);
				*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
				*(bds+N*3+0) =th2; *(bds+N*3+1) =a_sl(met,th2,-1); *(bds+N*3+2) =a_sl(met,th2,1); N++;

			} else  if ( th2==xs[k+1] && th2!=thb) {
				*(bds+N*3+0) =th2; *(bds+N*3+1) =a_sl(met,th2,-1); *(bds+N*3+2) =a_sl(met,th2,1); N++;
			}

			if (th2==xs[0] && th2!=thb) {k--; th2=min(xs[1],thb);} 
				else if (th2==xs[ns-1] && th2!=thb) {th2=thb;}
					else if (th2!=thb) th2 =min(xs[k+2],thb);
		}


// confidence region ends at "thb" analagously to its start at "tha", that is 
//      in M1 or M2  with a vertical line if thb =x(1) or x(n), at an 
// alpha-MLE point if  x(1) < thb < x(n),  or an open end if  thb > x(n)

		if (fabs(thb-xs[ns-1])<zero_eq) {
			th = xs[ns-1]*(1.-rel_print_eps);
			high_a = ahigh(met,th);
			*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
			*(bds+N*3+0) =thb; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
		}

		else if (fabs(thb-xs[0])<zero_eq) {
			th = xs[0]*(1.-rel_print_eps);
			high_a = ahigh(met,th);
			*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
			*(bds+N*3+0) =thb; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
		}

		else if (thb < xs[ns-2] + zero_eq) {
			high_a = ahigh(met,thb);
			*(bds+N*3+0) =thb; *(bds+N*3+1) =high_a; *(bds+N*3+2) =high_a; N++;
		}

		else {
			*(bds+N*3+0) =th; *(bds+N*3+1) =a_sl(met,th,-1); *(bds+N*3+2) =a_sl(met,th,1); N++;
		}
	}


	if(output)  { 
		if(met==GEO) Rcout << endl << endl;
		Rcout << " " << 100*(1-SL) << _("-percent joint confidence region for  (theta, alpha)  by method  ");
		if(met==GEO)  Rcout << "CLR" << endl; else  Rcout << "AF" << endl;
		Rcout << endl << setw(10) << "theta" << setw(15) << "alpha-min" << setw(15) 
				<< "alpha-max" << endl << endl;
		for(i=0;i<N;i++)  {
			Rcout << setw(10) << *(bds+i*3+0) << "," << setw(14) << *(bds+i*3+1)
				<< "," << setw(14) << *(bds+i*3+2) << endl;
			for(int j=0;j<numr;j++)  if( *(bds+i*3+0) == th_bds[2*j+1] )  Rcout << endl;
			if( *(bds+i*3+0) > th_bds[2*numr-1] )  Rcout << endl;
		}
		Rflush(); 
	}


	if(bounds!=0) for(i=0;i<N*3;i++) *(bounds+i) = *(bds+i);


	delete[] th_bds; delete[] bds;
	return N;
}




double Cblmr::a_sl(METHOD met, double th, int high_low)
// return alpha high or low value for a given SL value and a given theta0
// if high_low > 0 return high value, if high_low < 0 return low value
// if met=GEO use AF alpha as first guess, then grid search + bisection
{
	if (met==AF)  {
		return a_af(th, high_low);  
	} else {
		if(th!=old_th) { ah =ahigh(met,th);  old_th =th; }
		if( sl(th,ah,met,false) < SL )  stop( _("'a_sl' initial point below critical SL") );
		double incr = inc_y*high_low, guess = a_af(th, high_low);
		if ( isnan(guess) || (guess-ah)*high_low < zero_eq) guess = ah + incr;
		double guess2 = ah;
		while ( sl(th,guess,met,false) > SL) { guess2 = guess;  guess += incr; }
		double a_geo =bisect(guess,guess2,&Cblmr::sl_a,1,SL,-acc_yb);
		return a_geo;
	}
}



double Cblmr::a_af(double th, int high_low)
// return alpha high or low value for a given SL and given theta0 by AF
// alpha-boundaries by AF are roots of a quadratic
// if 'high_low' > 0 return high value, if 'high_low' < 0 return low value
{

	if (th != prev_th) {
		prev_th = th;

		bool th_ex;
		if ( (Model==M1 && th<=xs[0]) || xs[ns-1]<=th) th_ex=true; else th_ex=false;

		double  a, bp, c, xc;
		if (th_ex) {
			if (variance_unknown) xc=x_vu_ex; else xc=x_vk_ex;
			if(Model==M1) {
				double ff = sxx-2*sx1*th+s11*th*th, fy = yx-th*y1, f1 = sx1-th*s11;
				a =s11-f1*f1/ff; bp =y1-fy*f1/ff; c =sysq-fy*fy/ff -xc;
			} else {
				a= s11; bp= y1; c= sysq - xc;
			}
		} else {					
			if (variance_unknown) xc=x_vu; else xc=x_vk;
			int k=0; while (xs[k+1]<th) k++;
			double  fry = *psy*gfr(th,k), fr1 = *psig1*gfr(th,k);
			if(Model==M1) {
				double smy = *psy*gsm(th,k), sm1=*psig1*gsm(th,k); 
				a =s11-sm1*sm1-fr1*fr1; bp =y1-smy*sm1-fry*fr1; c =sysq-smy*smy-fry*fry -xc;
			} else {
				a= s11-fr1*fr1; bp= y1-fry*fr1; c= sysq - fry*fry - xc;
			}
		}
		double rad = bp*bp-a*c, rrad; 
		if(fabs(rad)<zero_eq) rad=0.;
		if(rad<0) rrad=und; else rrad= sqrt(rad);
		a_low = (bp-rrad)/a;
		a_high = (bp+rrad)/a;
	}

	if (high_low <0) return a_low; else return a_high;
}



double Cblmr::ahigh(METHOD met, double th)
// return 'alpha' value that gives the highest significance level for a given theta
{
	bool th_ex;
	if ( (Model==M1 && th<=xs[0]) || xs[ns-1]<=th) th_ex=true; else th_ex=false;

	double amle;
	if (th_ex) {
		if(Model==M1) {
			double ff = sxx-2*sx1*th+s11*th*th, fy = yx-th*y1, f1 = sx1-th*s11;
			amle = y1/s11 - f1*(fy-f1*y1/s11)/(ff*s11-f1*f1);
		} else   amle = y1/s11;
	} else {
		int k=0; while (xs[k+1]<th) k++;
		double fr1 = *psig1*gfr(th,k), fry = *psy*gfr(th,k);
		if(Model==M1) {
			double sm1 = *psig1*gsm(th,k), smy = *psy*gsm(th,k);
			amle = (y1-fry*fr1-smy*sm1)/(s11-fr1*fr1-sm1*sm1);
		} else   amle = (y1 - fr1*fry)/(s11-fr1*fr1);
	}
	double high_alpha = amle;


	if ( met==GEO  &&  !th_ex ) {

		if(th!=th0) set_theta0(th);

// find maxima of sl_geo2 by grid search, starting from amle
		double incr=inc_y, ta=amle; 
		set_alpha0(ta);
		double tsl1=sl_geo2();
		ta +=incr;
		set_alpha0(ta);
		double tsl2=sl_geo2();

		if(tsl1 > tsl2) {
			incr = -incr;
			ta = amle;
			tsl2 = tsl1;
			tsl1 = tsl2-1.;
		}


		if(tsl2 < SL) {

			while( fabs(incr) > acc_yb ) {

				while (tsl2 > tsl1) {
					tsl1 = tsl2;
					ta += incr;
					set_alpha0(ta);
					tsl2=sl_geo2();
				}
				incr = -incr/10;
				if(tsl1 < SL/64) break;
				tsl1 = tsl2 - 1.;
			}

			ta -= incr*10;
		}

		high_alpha = ta;
	}

	return high_alpha;
}



double Cblmr::sl_a(double alpha, int k)
// wrapper for use in the bisection routine
// return sl_geo2(th0,alpha)
// assume th0 preset and use default accuracy
{
	set_alpha0(alpha, GEO2);
	return sl_geo2();
}

