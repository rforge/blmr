//

#include "blmr.h"




double Cblmr::sl(const double th_0, const METHOD met, const bool output)
{
	double sL;

	if(trivial) { 
		const double thmle= mle(false);
		if( isnan(thmle) || th_0==thmle || (thmle==x_d[0] && th_0 <= thmle) ) sL= 1.;  else  sL= 0.;
		if (output)  {
			if( model_in != -2) 
				Rcout << _("  for theta0= ") << th_0 << "  SL= " << sL << endl << endl;
			else
				Rcout << _("  for theta0= ") << -th_0 << "  SL= " << sL << endl << endl;
		}

	}  else  {

		double  er=0., *const per= &er;

		set_theta0(th_0, met);

		if (met==AF) sL = sl_af();
		if (met==AF2) sL = sl_af(2);
		if (met==GEO) sL = sl_geo1(per);
		if (met==GEO2) {
			const double ah = ahigh(GEO,th_0);
			set_alpha0( ah, met );
			sL = sl_geo2(per);
		}
		if (met==MC) sL = sl_mc();

		if (output) {
			if( model_in != -2) 
				Rcout << _("  for theta0= ") << th_0 << "  SL= " << sL << _("  by method ");
			else
				Rcout << _("  for theta0= ") << -th_0 << "  SL= " << sL << _("  by method ");
			if (met==AF)  Rcout << "AF" << endl;  
			if (met==GEO)  Rcout << "CLR  int.er.< " << *per << endl;  
			if (met==MC)  Rcout << "CLR-MC" << endl;
			Rcout << endl;
		}
	}

	return sL ;
}



double Cblmr::sl(const double th_0, const double a0, const METHOD met, const bool output)
{
	set_theta0(th_0, met);
	set_alpha0(a0, met);

	double sL;

	if(trivial) { 
		const double thmle= mle(false);
		if( isnan(thmle)  ||  (thmle==x_d[0] && th_0 <= thmle) ) {
			const double  slope= ( (*py)[ i_d[1] ] - (*py)[ i_d[0] ] )/( x_d[1]-x_d[0]),  intercept= (*py)[ i_d[0] ]  - slope*x_d[0],
				amle = slope*th_0 + intercept;
			if( fabs(a0 - amle) < zero_eq )  sL= 1.;  else  sL= 0.;
		}  else  {
			if( lambdasq < zero_eq )  sL= 1.;  else  sL= 0.;
		}	
		if (output)  {
			if( model_in != -2 )
			Rcout << _(" for (th0,a0)= ( ") << th_0 << ", " << a0 << " )   SL= " << sL << endl << endl; 
			else
			Rcout << _(" for (th0,a0)= ( ") << -th_0 << ", " << a0 << " )   SL= " << sL << endl << endl; 
		}

	}  else  {

		double  er=0., *const per= &er;

		if (met==AF) sL = sl_af2();
		if (met==GEO) sL = sl_geo2(per); 
		if (met==MC) sL = sl_mc2();


		if (output) {
			if( model_in != -2 )
			Rcout << _(" for (th0,a0)= ( ") << th_0 << ", " << a0 << " )   SL= " << sL << _("  by method "); 
			else
			Rcout << _(" for (th0,a0)= ( ") << -th_0 << ", " << a0 << " )   SL= " << sL << _("  by method "); 
			if (met==AF) Rcout << "AF" << endl;
			if (met==GEO) Rcout << "CLR  int.er.< " << *per << endl;
			if (met==MC) Rcout << "CLR-MC" << endl;
			Rcout << endl;
		}
	}

	return sL;
}

