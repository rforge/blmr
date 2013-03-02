//

#include "blmr.h"



int Cblmr::cr(METHOD met, double inc, bool output, double *bounds)
// checks whether { (theta,alpha)  such that  sig. level  > SL}  is contiguous
// returns N = number of rows for bounds[N][3] array of (th, low_alpha, high_alpha) values
// if(met==GEO) using Knowles, Siegmund and Zhang's geometric formula to calculate  sig. levels
// if(met==AF) using Approximate F method
//
//		The cr routine determines the theta-boundaries by using 
//	the ci subroutine to check the 2-parametre significance levels along 
//	the  ( theta, MLE-alpha(theta) )  ridge,  then finds the alpha-boundaries
//	at increments of theta0.
//
//		By AF, confidence regions can be disjoint verses theta, but not verses alpha.
//	For a given theta0 in the broken line model LR vs alpha has a single maxima.  By AF 
//  the alpha-values for a desired SL are roots of a quadratic.  Hence the confidence
//  region is contiguous versus alpha.
//		Pending proof, this routine assumes  SL-by-CLR vs alpha  also has a single maxima.
//
{

	int i, numr, k, N=0;
	double  th =0;
	Rcpp::Function Rflush("flush.console");


// get theta-boundaries of confidence region(s)
	if(output) if(met==GEO) { Rcout << "               getting theta-boundaries..."<< endl; Rflush();}
	METHOD  met2;
	if(met==GEO) met2 =GEO2; else met2 =AF2;
	double *tmp, *bd1;
	tmp = new double[2*ns];
	numr = ci(met2,false,tmp);
	bd1 = new double[2*numr];
	for (i=0;i<2*numr;i++) bd1[i] = tmp[i];
	delete[] tmp; tmp=NULL;

	if(bd1[0]== -inf) bd1[0]= Xs[0]-1;
	if(bd1[2*numr-1]== +inf) bd1[2*numr-1]= Xs[ns-1]+1;


// store boundary values in an  N x 3  array	
	double width=0;						
	for (i=0;i<numr;i++) width += bd1[2*i+1] - bd1[2*i];
	int Nmax = width/inc + 1 + ns + 2*numr + 2;
	double *bd2;
	bd2 = new double[Nmax*3];

	if (output) {
		Rcout << endl << "  " << 100*(1-SL) << "% joint confidence region for  (theta, alpha)  by ";
		if (met==GEO)  Rcout << "CLR:" << endl; else  Rcout << "AF:" << endl;
		Rcout << endl << setw(10) << "theta" << setw(13) << "low-alpha" << setw(16) 
				<< "high-alpha" << endl;
		Rflush();
	}



	for (i=0;i<numr;i++) {

		if(output) { Rcout << endl; Rflush(); }

		double tha = bd1[2*i], thb = bd1[2*i+1];

		if ( fabs(thb-tha) < zero_eq) continue;

// confidence region starts at "tha"  with a vertical line if tha=x(1) or x(n), 
// an alpha-MLE point if x(2)< tha <x(n-1), or an open-end if tha < x(1)

		double high_a;
		if (fabs(tha-Xs[0])<zero_eq) {
			th = Xs[0]*(1.+rel_print_eps);
			high_a = ahigh(met,th);
			*(bd2+N*3+0) =tha; *(bd2+N*3+1) =high_a; *(bd2+N*3+2) =high_a; 
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
			*(bd2+N*3+0) =th; *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1); 
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
		}

		else if (fabs(tha-Xs[ns-1])<zero_eq) {
			th = Xs[ns-1]*(1.+rel_print_eps);
			high_a = ahigh(met,th);
			*(bd2+N*3+0) =tha; *(bd2+N*3+1) =high_a; *(bd2+N*3+2) =high_a; 
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
			*(bd2+N*3+0) =th;  *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1); 
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
		}

		else if (tha > Xs[k1] - zero_eq) {
			high_a = ahigh(met,tha);
			*(bd2+N*3+0) =tha; *(bd2+N*3+1) =high_a; *(bd2+N*3+2) =high_a;
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
		}

		if(output) Rflush();

		double th1 = floor(tha);
		if (tha < Xs[0]-zero_eq) th1 =tha; else while (th1<tha+zero_eq) th1 +=inc;
		if( fabs(th1)<zero_eq) th1=0.;

		int k_1=0; while (Xs[k_1+1]<th1 && k_1<ns-1) k_1++;
		int kb=0; while (Xs[kb+1]<thb && kb<ns-1) kb++;
			
		double th2;
		if (th1 < Xs[0]-zero_eq) th2 = Xs[0]; else th2 = min(Xs[k_1+1],thb);

		for (k=k_1;k<=kb;k++) {

			for (th = th1;th < th2*(1.-rel_print_eps); th += inc) {
				if( fabs(th)<zero_eq) th=0.;
				*(bd2+N*3+0) =th; *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1); 
				if(output) { Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
					<< "," << setw(15) << *(bd2+N*3+2) << endl; Rflush();}  N++;
			}

			th1 = th;
			if ( fabs(th1-th2) < rel_print_eps*(th1+th2)/2. )  th1 = th+inc;

			if ( th2==Xs[0] && th2!=thb) {
				*(bd2+N*3+0) =th2; *(bd2+N*3+1) =a_sl(met,th2,-1); *(bd2+N*3+2) =a_sl(met,th2,1); 
				if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
					<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
				th = th2*(1.+rel_print_eps);
				*(bd2+N*3+0) =th; *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1);
				if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
					<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;

			} else  if ( th2==Xs[ns-1] && th2!=thb) {
				th = th2*(1.-rel_print_eps);
				*(bd2+N*3+0) =th; *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1);
				if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
					<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
				*(bd2+N*3+0) =th2; *(bd2+N*3+1) =a_sl(met,th2,-1); *(bd2+N*3+2) =a_sl(met,th2,1);
				if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
					<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;

			} else  if ( th2==Xs[k+1] && th2!=thb) {
				*(bd2+N*3+0) =th2; *(bd2+N*3+1) =a_sl(met,th2,-1); *(bd2+N*3+2) =a_sl(met,th2,1);
				if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
					<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
			}

			if(output) Rflush();

			if (th2==Xs[0] && th2!=thb) {k--; th2=Xs[1];} 
				else if (th2==Xs[ns-1] && th2!=thb) {th2=thb;}
					else if (th2!=thb) th2 =min(Xs[k+2],thb);
		}


// confidence region ends at "thb" analagously to start at "tha", that is with a vertical line
// if thb =x(1) or x(n), at an alpha-MLE point if  x(1) < thb < x(n),  or an open end if  thb > x(n)

		if (fabs(thb-Xs[ns-1])<zero_eq) {
			th = Xs[ns-1]*(1.-rel_print_eps);
			*(bd2+N*3+0) =th; *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1);
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
			high_a = ahigh(met,thb);
			*(bd2+N*3+0) =thb; *(bd2+N*3+1) =high_a; *(bd2+N*3+2) =high_a;
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
		}

		else if (fabs(thb-Xs[0])<zero_eq) {
			th = Xs[0]*(1.-rel_print_eps);
			*(bd2+N*3+0) =th; *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1);
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
			high_a = ahigh(met,th);
			*(bd2+N*3+0) =thb; *(bd2+N*3+1) =high_a; *(bd2+N*3+2) =high_a;
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
		}

		else if (thb < Xs[ns-2] + zero_eq) {
			high_a = ahigh(met,thb);
			*(bd2+N*3+0) =thb; *(bd2+N*3+1) =high_a; *(bd2+N*3+2) =high_a;
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
		}

		else {
			*(bd2+N*3+0) =th; *(bd2+N*3+1) =a_sl(met,th,-1); *(bd2+N*3+2) =a_sl(met,th,1);
			if(output) Rcout << setw(10) << *(bd2+N*3+0) << "," << setw(12) << *(bd2+N*3+1)
				<< "," << setw(15) << *(bd2+N*3+2) << endl;  N++;
		}
	}

	if(output)  { Rcout << endl << endl; Rflush(); }


	if(bounds!=0) for(i=0;i<N*3;i++) *(bounds+i) = *(bd2+i);


	delete[] bd1; delete[] bd2;
	return N;
}




double Cblmr::a_sl(METHOD met, double th, int high_low)
// return alpha high or low value for a given SL value and a given theta0
// if SL > 0 return high value, if SL < 0 return low value
// if met=GEO use AF alpha as first guess, then grid search + bisection
{
	if (met==AF)  return a_af(th, high_low);  else {
		if(th!=old_th) { ah =ahigh(met,th);  old_th =th; }
		double incr = inc_y*high_low, guess1 = a_af(th, high_low);
		if ( isnan(guess1) ) guess1 = ah;
		double sl1 = sl(th,guess1,met,false), sl2=sl1, guess2=guess1;
		if (sl1 < SL) incr = -incr;
		while ( (sl1-SL)*(sl2-SL)>0 ) {
			sl1=sl2;
			guess1=guess2;
			guess2 +=incr;
			if ( (guess1-ah)*(guess2-ah)<0 ) guess2=ah;
			sl2 =sl(th,guess2,met,false);
		}
		double a_geo =bisect(guess1,guess2,&Cblmr::sl_a,1,SL,-acc_yb);
		return a_geo;
	}
}



double Cblmr::a_af(double th, int high_low)
// return alpha high or low value for a given SL and given theta0 by AF
// alpha-boundaries by AF are roots of a quadratic
// if SL > 0 return high value, if SL < 0 return low value
{

	if (th != prev_th) {
		prev_th = th;

		bool th_ex;
		if ( (model==M1 && th<=Xs[0]) || Xs[ns-1]<=th) th_ex=true; else th_ex=false;

		double  a, bp, c, xc;
		if (th_ex) {
			if (variance_unknown) xc=x_vu_ex; else xc=x_vk_ex;
			if(model==M1) {
				double ff = sxx-2*sx1*th+s11*th*th, fy = yx-th*y1, f1 = sx1-th*s11;
				a =s11-f1*f1/ff; bp =y1-fy*f1/ff; c =sysq-fy*fy/ff -xc;
			} else {
				a= s11; bp= y1; c= sysq - xc;
			}
		} else {					
			if (variance_unknown) xc=x_vu; else xc=x_vk;
			int k=0; while (Xs[k+1]<th) k++;
			double  fry = *psy*gfr(th,k), fr1 = *psig1*gfr(th,k);
			if(model==M1) {
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
// return alpha value that gives the highest significance level for a given theta
{
	bool th_ex;
	if ( (model==M1 && th<=Xs[0]) || Xs[ns-1]<=th) th_ex=true; else th_ex=false;

	double amle;
	if (th_ex) {
		if(model==M1) {
			double ff = sxx-2*sx1*th+s11*th*th, fy = yx-th*y1, f1 = sx1-th*s11;
			amle = y1/s11 - f1*(fy-f1*y1/s11)/(ff*s11-f1*f1);
		} else   amle = y1/s11;
	} else {
		int k=0; while (Xs[k+1]<th) k++;
		double fr1 = *psig1*gfr(th,k), fry = *psy*gfr(th,k);
		if(model==M1) {
			double sm1 = *psig1*gsm(th,k), smy = *psy*gsm(th,k);
			amle = (y1-fry*fr1-smy*sm1)/(s11-fr1*fr1-sm1*sm1);
		} else   amle = (y1 - fr1*fry)/(s11-fr1*fr1);
	}
	double high_alpha = amle;


	if ( met==GEO && !th_ex && variance_unknown ) {

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
// return sl_geo2(th0,alpha), wrapper for use in the bisection routine
// assume th0 preset and use default accuracy
{
	set_alpha0(alpha, GEO2);
	return sl_geo2();
}

