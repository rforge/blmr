//

#include "blmr.h"




double Cblmr::sl_mc(void)
// calculate significance level by CLR, Monte Carlo evaluation method
{
	double acc= acc_sl_abs;
	Rcpp::Function Rflush("flush.console");

	Rcout << endl << "MC evaluation of CLR SL for theta0= " << th0 
		<< ",  target accuracy =  " << acc << ":" << endl << endl;
	Rcout << setw(10) << "iteration" << setw(12) << "est. SL" << setw(15) << "est. acc." << endl;
	Rflush();

	int  it, N =10000000,  count =0;
	double tstart =time(NULL), wsq =w*w, testw, s0;

	if (variance_unknown) { s0 = z/sqrt( 1.- z*z ); testw = wsq/(1-z*z);} 
					else  { s0 = z; testw = wsq;}


	Vector<double> s(m), &rs = s;

	for (it=1; it<N+1; it++)
	{
		if (th0ex)  s[0] =rnorm(); else  s[0] =0.;
		for (int i=1;i<m;i++) s[i] =rnorm();
		if (variance_unknown) s = 1./sqrt(s*s) * s;
		if (!th0ex) s[0] =s0;

		if ( m_gt_w(testw, rs) ) count++;

		if ( it==N/10 || !(it%(N/5)) ) {
			double  p = 1.*count/it,  err_est = 2*sqrt( p*(1.-p)/it );
			Rcout << setw(10) << it << setw(12) << p << setw(15) << err_est << endl; Rflush();
			if (err_est < acc &&  it >= N/5)  {it++; break;}
		}
	}
	it--;

	double sl = count*1./it;

		double tfinish = time( NULL );
		double elapsed = tfinish - tstart;
		Rcout << endl << "elapsed  " << elapsed << " sec  =  " << elapsed/60. << " min" << endl;
		Rcout << endl;

	return sl; 
}





double Cblmr::sl_mc2(void)
// calculate significance level by clr, Monte Carlo evaluation method
// for changepoint = th0 and   alpha = alpha0 
{
	double acc= acc_sl_abs;
	Rcpp::Function Rflush("flush.console");

	Rcout << endl << "MC evaluation of CLR SL for (th0,a0)= (" << th0 
		<< "," << alpha0 << "),  target accuracy = " << acc << ":" << endl << endl;
	Rcout << setw(10) << "iteration" << setw(14) << "est. SL" << setw(15) << "est. acc." << endl;
	Rflush();

	double  Fc=0, sl=0; 
	if (variance_unknown) {if (th0ex) Fc =F(m,-c); else  Fc =F(m-1,-c);} 
			else  Fc = Rf_pnorm5(-lambda*c ,0,1,1,0) ;

// generate mock results

	int it, count=0, N=10000000;
	double  sum=0., sumsqs=0., tstart = time( NULL ), s0, testw, den;
	Vector<double> s(m), &rs = s;

	for (it=1; it<N+1; it++)
	{
		double xi = c*2*(unif_rand()-0.5), wsq, z_tilde, deltasq;

		if (th0ex) z_tilde=0.; else  if (model==M1) z_tilde = xi*c1 + c2; else z_tilde =z_tilde_M2; 

		deltasq = lambdasq*(1-xi*xi) + z_tilde*z_tilde;

		if (th0ex) z =0.; else  if (variance_unknown) z =z_tilde/sqrt(deltasq); else z =z_tilde;

		if (variance_unknown) wsq = 1 - omega/deltasq; else wsq = deltasq - omega;
		w = sqrt(max(wsq,0.));


		if (variance_unknown) { s0 = z/sqrt(1-z*z); testw = wsq/(1-z*z); }
						else  { s0 = z; testw = wsq; }


		if (th0ex)  s[0] =rnorm();  else  s[0] =0.;
		for (int i=1;i<m;i++) s[i] = rnorm();
		if (variance_unknown) s =  1./sqrt(s*s) * s;
		if (!th0ex) s[0] = s0;

		if ( m_gt_w(testw, rs) ) {
			count++;
			if (variance_unknown) {if (th0ex) den =fk(m,xi); else  den =fk(m-1,xi);} 
					else  den = Rf_dnorm4(lambda*xi, 0,1,0) ;
			sum += den;
			sumsqs += den*den;
		}


		if ( it==N/10 || !(it%(N/5)) ) {
			double  p = 1.*sum/it;
			double sd = sqrt( (sumsqs/it - p*p)/it );
			double err_est = 2*c*sd;
			if (!variance_unknown)  err_est *= lambda;
			if (variance_unknown) sl = 2*Fc + 2*c*sum/it; else sl = 2*Fc + lambda*2*c*sum/it;
			Rcout << setw(10) << it << setw(14) << sl << setw(15) << err_est << endl; Rflush();
			if (err_est < acc && it >=N/5)  {it++; break;}
		}
	}
	it--;

		double tfinish = time( NULL );
		double elapsed = tfinish - tstart;
		Rcout << endl << "elapsed  " << elapsed << " sec  =  " << elapsed/60. << " min" << endl;
		Rcout << endl;

	return sl;
}

