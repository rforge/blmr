//

#include "lmbr.h"




double Clmbr::sl(const double th_0, const METHOD met, const bool output)
{
	double sL;

	if(trivial) { 
		const double thmle= mle(false);
		if( (isnan(thmle) && !isinf(thmle)) || th_0==thmle || 
				(Model==M2 && thmle==xs[0] && th_0 <= thmle) )  sL= 1.;  else  sL= 0.;
		if (output)  {
			if( model_in > 0 ) 
				Rcout << "  SL= " << sL << _("  for theta0 = ") << th_0 << endl << endl;
			else
				Rcout << "  SL= " << sL << _("  for theta0 = ") << -th_0 << endl << endl;
		}

	}  else  {

		set_theta0( th_0, met );

		double  err=0., *const per= &err;

		if ( fabs(z) >= w )  sL= 1.;  else  {

			if (met==AF) sL = sl_af();
			if (met==AF2) sL = sl_af(2);
			if (met==GEO) sL = sl_geo(per);
			if (met==GEO2) {
				const double ah = ahigh(GEO,th_0);
				set_alpha0( ah, met );
				sL = sl_geo2(per);
			}
			if (met==MC) sL = sl_mc();
		}

		if (output) {
			if( model_in > 0 ) 
				Rcout << "  SL= " << sL << _("  for theta0 = ") << th_0 << _("  by method ");
			else
				Rcout << "  SL= " << sL << _("  for theta0 = ") << -th_0 << _("  by method ");
			if (met==AF)  Rcout << "AF" << endl;  
			if (met==GEO)  Rcout << "CLR  int.er.< " << err << endl;  
			if (met==MC)  Rcout << "CLR-MC" << endl;
			Rcout << endl;
		}
	}

if( isinf(sL) || isnan(sL) ) {
  Rcout << "sL= " << sL << endl;
  Rcout << "th0  " << th0 << endl;
  Rcout << "y  " << *py << endl;
  stop( "debug" );
}

	return sL ;
}



double Clmbr::sl(const double th_0, const double a0, const METHOD met, const bool output)
{
	double sL;

	if(trivial) { 
		const double thmle= mle(false);
		if( (isnan(thmle) && !isinf(thmle))  ||  (thmle==xs[0] && th_0 <= thmle) ) {
			const double  slope= ( (*py)[ is[1] ] - (*py)[ is[0] ] )/( xs[1]-xs[0]),  intercept= (*py)[ is[0] ]  - slope*xs[0],
				amle = slope*th_0 + intercept;
			if( fabs(a0 - amle) < zero_eq )  sL= 1.;  else  sL= 0.;
		}  else  {
			if( lambdasq < zero_eq )  sL= 1.;  else  sL= 0.;
		}	
		if (output)  {
			if( model_in > 0 )
			Rcout << "  SL= " << sL << _(" for (th0,a0)= ( ") << th_0 << ", " << a0 << " ) " << endl << endl; 
			else
			Rcout << "  SL= " << sL << _(" for (th0,a0)= ( ") << -th_0 << ", " << a0 << " ) " << endl << endl; 
		}

	}  else  {

		set_theta0(th_0, met);
		set_alpha0(a0, met);

		double  er=0., *const per= &er;

		if (met==AF) sL = sl_af2();
		if (met==GEO) sL = sl_geo2(per); 
		if (met==MC) sL = sl_mc2();


		if (output) {
			if( model_in > 0 )
			Rcout << "  SL= " << sL << _(" for (th0,a0)= ( ") << th_0 << ", " << a0 << " ) " << _("  by method "); 
			else
			Rcout << "  SL= " << sL << _(" for (th0,a0)= ( ") << -th_0 << ", " << a0 << " ) " << _("  by method "); 
			if (met==AF) Rcout << "AF" << endl;
			if (met==GEO) Rcout << "CLR  int.er.< " << *per << endl;
			if (met==MC) Rcout << "CLR-MC" << endl;
			Rcout << endl;
		}
	}
if( isinf(sL) || isnan(sL) ) {
  Rcout << "sL= " << sL << endl;
  Rcout << "th0  " << th0 << endl;
  Rcout << "y  " << *py << endl;
  stop( "debug" );
}
	return sL;
}


