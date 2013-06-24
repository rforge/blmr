//

#include "blmr.h"

#define CLmsg			_( "confidence level must be between 0 and 1" )
#define methods2msg		_( "'method' must be \"clr\" or \"af\",  example:  ci( 0.99, \"af\" )" )
#define methods3msg		_( "'method' must be \"clr\", \"mc\" or \"af\",  example:  sl( 5.8, \"af\" )" )
#define model_msg		_( "not applicable for this model" )



void Cblmr::slR(const double theta0) { 
	if(model_in > 0 ) sl(theta0); else sl(-theta0);
	return; 
}


void Cblmr::slR(const double theta0, const double alpha0) {
	if(Model==M3)  stop( model_msg );  
	if(model_in > 0 ) sl(theta0,alpha0); else sl(-theta0,alpha0);
	return; 
}


void Cblmr::slR(const double theta0, const string met ) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	if(model_in > 0 ) sl(theta0, MET ); else sl(-theta0, MET );

	return;
}


void Cblmr::slR(const double theta0, const double alpha0, const string met ) { 

	if(Model==M3)  stop( model_msg );  

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	if(model_in > 0 ) sl(theta0, alpha0, MET ); else sl(-theta0, alpha0, MET );

	return;
}



void Cblmr::slR(const double theta0, const string met, const double acc) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	const double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	if(model_in > 0 ) sl(theta0, MET ); else sl(-theta0, MET );

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return;
}



void Cblmr::slR(const double theta0, const double alpha0, const string met, const double acc) { 

	if(Model==M3)  stop( model_msg );  

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	const double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	if(model_in > 0 ) sl(theta0, alpha0, MET ); else sl(-theta0, alpha0, MET ); 

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return;
}




double Cblmr::slR(const double theta0, const string met, const double acc, const bool output) { 

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	const double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	double result;
	if( model_in  > 0 ) 
		result= sl(theta0, MET, output);
	else
		result= sl(-theta0, MET, output);

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return result;
}



double Cblmr::slR(const double theta0, const double alpha0, const string met, const double acc, const bool output) { 

	if(Model==M3)  stop( model_msg );  

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			if(met=="MC" || met=="mc") MET=MC; else {
				stop( methods3msg );
			}
		}
	}

	const double tmp1 = acc_sl_abs, tmp2 = acc_sl_rel;
	acc_sl_abs= acc;
	acc_sl_rel= min(10*acc_sl_abs,0.01);

	double result;
	if( model_in  > 0 ) 
		result= sl(theta0, alpha0, MET, output);
	else
		result= sl(-theta0, alpha0, MET, output);

	acc_sl_abs= tmp1;
	acc_sl_rel= tmp2;

	return result;
}




void Cblmr::ciR(const double CL) { 

	if(CL <=0. || CL >=1.) stop( CLmsg );

	const double tmp = SL;
	set_SL(1.-CL);
	ci(); 
	set_SL(tmp);
	return; 
}



void Cblmr::ciR(const double CL, const string met) { 

	if(CL <=0. || CL >=1.) stop( CLmsg );

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			stop( methods2msg );
		}
	}

	const double tmp = SL;
	set_SL(1.-CL);
	ci(MET); 
	set_SL(tmp);

	return; 
}



void Cblmr::ciR(void) { 

	const double tmp = SL;
	set_SL(0.05);
	ci(); 
	set_SL(tmp);

	return;
}




void Cblmr::crR(const double CL) { 

	if(Model==M3)  stop( model_msg );  
	if(CL <=0. || CL >=1.) stop( CLmsg );

	const double tmp = SL;
	set_SL(1.-CL);
	cr(); 
	set_SL(tmp);

	return; 
}



void Cblmr::crR(const double CL, const string met) { 

	if(Model==M3)  stop( model_msg );  
	if(CL <=0. || CL >=1.) stop( CLmsg );

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			stop( methods2msg );
		}
	}

	const double tmp = SL;
	set_SL(1.-CL);
	cr(MET); 
	set_SL(tmp);

	return;
}



void Cblmr::crR(const double CL, const string met, const double incr) { 

	if(Model==M3)  stop( model_msg );  
	if(CL <=0. || CL >=1.) stop( CLmsg );

	METHOD MET;
	if(met=="CLR" || met=="clr") MET=GEO; else {
		if(met=="AF" || met=="af") MET=AF; else { 
			stop( methods2msg );
		}
	}

	const double tmp = SL;
	set_SL(1.-CL);
	cr(MET,incr); 
	set_SL(tmp);

	return;
}



void Cblmr::crR(void) { 

	if(Model==M3)  stop( model_msg );  

	const double tmp = SL;
	set_SL(0.05);
	cr(); 
	set_SL(tmp);

	return; 
}



void Cblmr::MLE(void)  const
{ 
	mle(); 
	return; 
}




void Cblmr::SET_irSy(const NumericVector irSy)  {

	const int yn =irSy.size();
	if(yn!=n) stop( _("'irSy' vector has wrong dimension") );

	double *const  Ytmp= new (nothrow) double[n];
	if(Ytmp==NULL) {
		Rcout << _("message: ") << 8 << endl;
		stop( _("memory allocation failed") );
	}
	
	for (int i=0;i<n;i++) Ytmp[i] = irSy[i];

	set_sy( Ytmp, GEO2 );

	delete[] Ytmp;

	return;
}


