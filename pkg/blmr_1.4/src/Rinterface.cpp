//

#include "blmr.h"

#define CLmsg   _( "confidence level must be between 0 and 1" )
#define methods2msg   _( "'method' must be \"clr\" or \"af\",  example:  ci( 0.99, \"af\" )" )
#define methods3msg   _( "'method' must be \"clr\", \"mc\" or \"af\",  example:  sl( 5.8, \"af\" )" )



void Cblmr::slR(double theta0) {  sl(theta0); return; }


void Cblmr::slR(double theta0, double alpha0) {  sl(theta0,alpha0); return; }


void Cblmr::slR(double theta0, string met ) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	sl(theta0, MET );

	return;
}


void Cblmr::slR(double theta0, double alpha0, string met ) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	sl(theta0, alpha0, MET );

	return;
}



void Cblmr::slR(double theta0, string met, double acc) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	sl(theta0, MET );

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return;
}



void Cblmr::slR(double theta0, double alpha0, string met, double acc) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	sl(theta0, alpha0, MET );

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return;
}




double Cblmr::slR(double theta0, string met, double acc, bool output) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	double result= sl(theta0, MET, output);

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return result;
}



double Cblmr::slR(double theta0, double alpha0, string met, double acc, bool output) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	double result= sl(theta0, alpha0, MET, output);

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return result;
}




void Cblmr::ciR(double CL) { 

	if(CL <0. || CL >1.) stop( CLmsg );

	double tmp = SL;
	set_SL(1.-CL);
	ci(); 
	set_SL(tmp);
	return; 
}



void Cblmr::ciR(double CL, string met) { 

	if(CL <0. || CL >1.) stop( CLmsg );

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			stop( methods2msg );
		}
	}

	double tmp = SL;
	set_SL(1.-CL);
	ci(MET); 
	set_SL(tmp);

	return; 
}



void Cblmr::ciR(void) { 

	double tmp = SL;
	set_SL(0.05);
	ci(); 
	set_SL(tmp);

	return;
}




void Cblmr::crR(double CL) { 

	if(CL <0. || CL >1.) stop( CLmsg );

	double tmp = SL;
	set_SL(1.-CL);
	cr(); 
	set_SL(tmp);

	return; 
}



void Cblmr::crR(double CL, string met) { 

	if(CL <0. || CL >1.) stop( CLmsg );

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			stop( methods2msg );
		}
	}

	double tmp = SL;
	set_SL(1.-CL);
	cr(MET); 
	set_SL(tmp);

	return;
}



void Cblmr::crR(double CL, string met, double incr) { 

	if(CL <0. || CL >1.) stop( CLmsg );

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			stop( methods2msg );
		}
	}

	double tmp = SL;
	set_SL(1.-CL);
	cr(MET,incr); 
	set_SL(tmp);

	return;
}



void Cblmr::crR(void) { 

	double tmp = SL;
	set_SL(0.05);
	cr(); 
	set_SL(tmp);

	return; 
}



void Cblmr::MLE(void) { mle(); return; }




void Cblmr::SET_irSy(NumericVector irSy) {

	int yn =irSy.size();
	if(yn!=n) stop( _("y vector has wrong dimension") );

	double *Ytmp;
	Ytmp = new double[n]; 
	for (int i=0;i<n;i++) Ytmp[i] = irSy[i];

	set_sy( Ytmp, GEO );

	delete[] Ytmp;

	return;
}


