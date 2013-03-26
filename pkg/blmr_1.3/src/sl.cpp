//

#include "blmr.h"




double Cblmr::sl(double th_0, METHOD met, bool output)
{
	set_theta0(th_0, met);

	double sl, er=0.,*per=&er;
	if (met==AF) sl = sl_af();
	if (met==AF2) sl = sl_af(2);
	if (met==GEO) sl = sl_geo1(per);
	if (met==GEO2) {
		set_alpha0( ahigh(GEO,th_0), met );
		sl = sl_geo2(per);
	}
	if (met==MC) sl = sl_mc();

	if (output) {
		Rcout << _("  for theta0= ") << th0 << "  SL= " << sl << _("  by method ");
		if (met==AF)  Rcout << "AF" << endl;  
		if (met==GEO)  Rcout << "CLR  int.er.< " << *per << endl;  
		if (met==MC)  Rcout << "CLR-MC" << endl;
		Rcout << endl;
	}

	return sl ;
}



double Cblmr::sl(double th_0, double a0, METHOD met, bool output)
{

	set_theta0(th_0, met);
	set_alpha0(a0, met);

	double sl, er=0.,*per=&er;
	if (met==AF) sl = sl_af2();
	if (met==GEO) sl = sl_geo2(per); 
	if (met==MC) sl = sl_mc2();


	if (output) {
		Rcout << _(" for (th0,a0)= ( ") << th0 << ", " << alpha0 << " )   SL= " << sl << _("  by method "); 
		if (met==AF) Rcout << "AF" << endl;
		if (met==GEO) Rcout << "CLR  int.er.< " << *per << endl;
		if (met==MC) Rcout << "CLR-MC" << endl;
		Rcout << endl;
	}


	return sl;
}

