//

#include "blmr.h"




double Cblmr::sl_mc(void)  const
// calculate significance level by CLR, Monte Carlo evaluation method
{
	Rcpp::Function Rflush("flush.console");
	const double acc= acc_sl_abs;
	GetRNGstate();

	Rcout << endl << _("MC evaluation of conditional likelihood-ratio SL") << endl;
	if(model_in != -2)
	  Rcout << _("for theta0= ") << th0 << _(",  target accuracy =  ") << acc << ":" << endl << endl;
	else
	  Rcout << _("for theta0= ") << -th0 << _(",  target accuracy =  ") << acc << ":" << endl << endl;
	Rcout << setw(10) << "iteration" << setw(12) << "est. SL" << setw(15) << "est. acc." << endl;
	Rflush();

	const double  tstart= time(NULL),  wsq= w*w; 

	double  s0, testw;
	if (variance_unknown) { s0 = z/sqrt( 1.- z*z );  testw = wsq/(1-z*z); } 
					else  { s0 = z; testw = wsq; }


	const int  N = 10000000;
	int  it,  count= 0;

	for (it=1; it<N+1; it++)
	{
		Vector<double>  s(m);
		if (th0ex)  s[0] =norm_rand(); else  s[0] =0.;
		for (int i=1;i<m;i++) s[i] =norm_rand();
		if (variance_unknown) s = 1./sqrt(s*s) * s;
		if (!th0ex) s[0] =s0;

		const Vector<double>  &rs = s;
		if ( m_gt_w(testw, rs) ) count++;

		if ( it==N/10 || !(it%(N/5)) ) {
			const double  p = 1.*count/it,  err_est = 2*sqrt( p*(1.-p)/it );
			Rcout << setw(10) << it << setw(12) << p << setw(15) << err_est << endl; Rflush();
			if (err_est < acc &&  it >= N/5)  {it++; break;}
		}
	}
	it--;

	const double sL = count*1./it;

	const double  tfinish = time( NULL ),  elapsed = tfinish - tstart;
	Rcout << endl << "elapsed  " << elapsed << " sec" << endl << endl;

	PutRNGstate();
	return sL; 
}





double Cblmr::sl_mc2(void)  const
// calculate significance level by clr, Monte Carlo evaluation method
// for changepoint = th0 and   alpha = alpha0 
{
	Rcpp::Function Rflush("flush.console");
	const double  acc= acc_sl_abs;
	GetRNGstate();
 
	if( model_in != -2 )
	  Rcout << endl << _( "MC evaluation of conditional likelihood-ratio SL for (th0,a0)= (" ) << th0 
		<< "," << alpha0 << _( "),  target accuracy = " ) << acc << ":" << endl << endl;
	else
	  Rcout << endl << _( "MC evaluation of conditional likelihood-ratio SL for (th0,a0)= (" ) << -th0 
		<< "," << alpha0 << _( "),  target accuracy = " ) << acc << ":" << endl << endl;
	Rcout << setw(10) << "iteration" << setw(14) << "est. SL" << setw(15) << "est. acc." << endl;
	Rflush();

	double  Fc=0, sL=0; 
	if (variance_unknown) {if (th0ex) Fc =F(m,-c); else  Fc =F(m-1,-c);} 
			else  Fc = Rf_pnorm5(-lambda*c ,0,1,1,0) ;

// generate mock results

	const double   tstart = time( NULL );

	const int  N = 10000000;
	int  it,  count= 0;
	double  sum=0., sumsqs=0.;

	for (it=1; it<N+1; it++)
	{
		const double xi = 2*c*(unif_rand()-0.5); 

		double z_tilde;
		if (th0ex)  z_tilde= 0.;  else  z_tilde = xi*c1 + c2; 

		const double  deltasq = lambdasq*(1-xi*xi) + z_tilde*z_tilde;

		double  z_;
		if (th0ex)  z_ = 0.;  else  { if (variance_unknown) z_ = z_tilde/sqrt(deltasq);  else  z_ = z_tilde; }

		double  wsq;
		if (variance_unknown)  wsq = 1 - omega/deltasq;  else  wsq = deltasq - omega;
		if (wsq<0.)  wsq= 0.;


		double  s0, testw;
		if (variance_unknown) { s0 = z_/sqrt(1-z_*z_); testw = wsq/(1-z_*z_); }
						else  { s0 = z_; testw = wsq; }


		Vector<double> s(m);
		if (th0ex)  s[0] = norm_rand();  else  s[0] = 0.;
		for (int i=1;i<m;i++) s[i] = norm_rand();
		if (variance_unknown) s =  1./sqrt(s*s) * s;
		if (!th0ex) s[0] = s0;

		const Vector<double>  &rs = s;
		if ( m_gt_w(testw, rs) ) {
			count++;
			double den;
			if (variance_unknown) {if (th0ex) den =fk(m,xi); else  den =fk(m-1,xi);} 
					else  den = Rf_dnorm4(lambda*xi, 0,1,0) ;
			sum += den;
			sumsqs += den*den;
		}


		if ( it==N/10 || !(it%(N/5)) )  {
			const double  p = 1.*sum/it;
			const double sd = sqrt( (sumsqs/it - p*p)/it );
			double err_est = 2*c*sd;
			if (!variance_unknown)  err_est *= lambda;
			if (variance_unknown) sL = 2*Fc + 2*c*sum/it; else sL = 2*Fc + lambda*2*c*sum/it;
			Rcout << setw(10) << it << setw(14) << sL << setw(15) << err_est << endl; Rflush();
			if (err_est < acc && it >=N/5)  {it++; break;}
		}
	}
	it--;
	if (variance_unknown)  sL = 2*Fc + 2*c*sum/it;  else  sL = 2*Fc + lambda*2*c*sum/it;

	const double  tfinish = time( NULL ),  elapsed = tfinish - tstart;
	Rcout << endl << "elapsed  " << elapsed << " sec" << endl << endl;

	PutRNGstate();
	return sL;
}

