//

#include "blmr.h"



double Cblmr::sl_geo(double *err)
// calculate significance level for changepoint = theta0 by CLR 
// using Knowles, Siegmund, Zhang's geometric formula
// assume model-object's  th0, y, z, w  values already set
{
	double  sL = 0., error=0., er=0, *per = &er;
	if (err!=0) *err=0.;


	if (th0ex) {

		sL = geo_ex(per);  error += *per;

	} else {

		if ( fabs(z) > w - zero_eq ) return 1.;

		if (th0 < xs[ns-2]) {
			sL += geo(xs[ns-2], per); error += *per;
			z = -z;
			sL += geo(xs[ns-2], per); error += *per;
		}
		if (th0 > xs[k1]) {
			sL += geo(xs[k1], per); error += *per;
			z = -z;
			sL += geo(xs[k1], per); error += *per;
		}
	}


	if (err!=0) *err = error;
	return min(sL,1.);
}



double Cblmr::sl_geo1(double *err)
{
		double wsq;
		if (variance_unknown)  wsq =1-omega/qysq; else  wsq = qysq - omega;
		w = sqrt(max(0.,wsq));

		if(th0ex) z=0.; else 
			{ if (variance_unknown)  z = star_z/sqrt(qysq); else z = star_z; }


		double sL, er=0.,*per=&er;

		sL = sl_geo(per);

		if(err!=0) *err = *per;


		return  sL;
}



double Cblmr::sl_geo2(double *err)
// calculate significance level for changepoint = (theta0,alpha0) by CLR 
// using Knowles, Siegmund, Zhang's geometric formula
// assume model-object's  th0, y  values already set
{
	double sL;
	if (variance_unknown) {if (th0ex) sL =2*F(m,-c); else  sL =2*F(m-1,-c);} 
		else  sL =2*Rf_pnorm5(-lambda*c ,0,1,1,0);

	double er=0, *per =&er;
//if(th0 > 5.608) for(double ts=-c;ts<c;ts+=2*c/20) Rcout << "xi ipr  " << ts << " " << ipr(ts,per) << endl; 
	double integral, error;
	if(Model==M1) { integral= integrate2(-c, c, &Cblmr::ipr, per); error= *per; }
	else 
		{ integral = 2*integrate2(-c, 0, &Cblmr::ipr, per); error= (*per)*2; }

	if (!variance_unknown) { integral *=lambda;  error *=lambda; }

	sL += integral;

	if(err!=0) *err = error;

	return min(sL,1.);
}




double Cblmr::ipr(double xi, double *err)
// integrand for  sl(theta,alpha)  geometric formula
{
	if(err!=0) *err =0.;

	double den;
	if(variance_unknown) {if (th0ex) den =fk(m,xi); else  den =fk(m-1,xi);} 
		else  den = Rf_dnorm4(lambda*xi ,0,1,0) ;

	if(xi==-c || xi==c) return den;

	double z_tilde =0;
	if( ! th0ex ) { if(Model==M1) z_tilde = xi*c1 + c2; else z_tilde = z_tilde_M2; }

	double deltasq = lambdasq*(1-xi*xi) + z_tilde*z_tilde;

	if(th0ex) z =0.; else { if(variance_unknown) z =z_tilde/sqrt(deltasq); else z =z_tilde;}

	double wsq;
	if(variance_unknown) wsq = 1 - omega/deltasq; else wsq = deltasq - omega;
	w = sqrt(max(wsq,0.));


	double  er=0,  *per=&er,  pr = sl_geo(per);

	if(err!=0) *err = (*per)*den;

	return  pr*den;
}




