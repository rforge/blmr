//

#include "lmbr.h"



double  Clmbr::mle( const bool output,  double *const  max_gy,  double *const  par  )  const
// find theta_mle for maximum value of (gamma(theta)*qy)^2
// calculate  alphamle, betamle, betapmle, vmle
{
	bool line= false;
	if( qysq/m < zero_eq )  line= true;


// find theta_mle, which is the value of theta that gives
// the maximum value of   ( gam(theta).Qy )^2 
	double  vf= (xs[ns-1]-xs[ns-2])*(*pqy*pq1[ns-1]), max= vf*vf/qff[ns-1], thmle= xs[ns-2];
	int kmle= ns-1;

	for (int k= ns-2; k > k1; k--) {

		const double  v1 = (*pqy)*pq1[k],  vx = vf + v1*xs[k];  
		if(k>0) vf = vx - v1*xs[k-1];

		const double  va = vx*qx1[k] - v1*qxx[k],  vb = vx*q11[k] - v1*qx1[k],  thk = va/vb;

		double mk;
		if(k>0)  {
			if(xs[k-1]<thk && thk<xs[k]) {
				mk = (vx*vb - v1*va)/ck[k]; 
				if (mk > max) {max = mk; thmle = thk; kmle=k;}
			} else { 
				mk = vf*vf/qff[k];
				if (mk > max) {max = mk; thmle = xs[k-1]; kmle=k;}
			}
		}  else  {
			if(  thk < xs[k]  ) {
				mk = (vx*vb - v1*va)/ck[k]; 
				if (mk > max) {max = mk; thmle = thk; kmle=k;}
			} else {
				const double  num= (*pqy)*gam(-Inf,0); 
				mk = num*num;
				if (mk > max)  { max = mk; thmle = -Inf; kmle= k; }
			}
		}

	}

	if(kmle>0)  if ( fabs(xs[kmle-1] - thmle) < zero_eq )  thmle = xs[kmle-1];
	if ( fabs(xs[kmle] - thmle) < zero_eq )  {  kmle++; thmle = xs[kmle-1]; }

	if( max_gy != NULL )  *max_gy =  max;


// calculate  alphamle, betamle, betapmle, vmle
	double alphamle, betamle, betapmle, vmle;

	if (Model==M1) {
		const double  amle = *psy*gbar(thmle,kmle);
		const double  bmle = *psy*gsm(thmle,kmle);
		const double  bpmle = *psy*gfr(thmle,kmle);

		double psi;
		if (cov_matrix_diagonal) psi=0.; else psi = gfr(thmle,kmle)*sfc(thmle,kmle);
		const double sm1 = *psig1*gsm(thmle,kmle);
		const double fr1 = *psig1*gfr(thmle,kmle);
		const double r1sf = sqrt(s11-sm1*sm1-fr1*fr1);
		alphamle = amle/r1sf;
		betamle = (bmle - alphamle*sm1)/sqrt(sfc(thmle,kmle)*sfc(thmle,kmle)-psi*psi);
		betapmle = (bpmle - alphamle*fr1 - betamle*psi)/sqrt(sf(thmle,kmle)*sf(thmle,kmle));

	}

	if (Model==M2) {

		const double sf1 = *pv1h*sf(thmle,kmle);
		const double y1 = *pv1h*(*psy);

		betamle = 0.;
		betapmle = ( *psy*sf(thmle,kmle) - sf1*y1 )/( sf(thmle,kmle)*sf(thmle,kmle) - sf1*sf1 );
		alphamle = (y1 - betapmle*sf1)/n1;

		if(model_in < 0) {
			betamle= -betapmle;
			betapmle= 0.;
		}
	}

	if (Model==M3) {

		alphamle = 0.;
		betamle = 0.;
		if( !isinf(thmle) )  betapmle =  ( *psy*sf(thmle,kmle) )/( sf(thmle,kmle)*sf(thmle,kmle) );  else  betapmle= 0.;

		if(model_in < 0) {
			betamle= -betapmle;
			betapmle= 0.;
		}
	}

	if (variance_unknown)  vmle = omega/(m-2);  else  vmle = 1.;

	if( fabs(thmle) < zero_eq) thmle= 0.;
	if( fabs(alphamle) < zero_eq) alphamle= 0.;
	if( fabs(betamle) < zero_eq) betamle= 0.;
	if( fabs(betapmle) < zero_eq) betapmle= 0.;
	if( fabs(vmle) < zero_eq) vmle= 0.;

	if(line) {
		if(Model==M1) {
			thmle= NaN;
			alphamle= NaN;
			const double  slope= ( (*py)[ n-1 ] - (*py)[ 0 ] )/( (*px)[n-1]-(*px)[0] );
			betamle = betapmle = slope;
		}
		if(Model==M2) {
			thmle= NaN;
			alphamle= (*py)[0];
			betamle = betapmle = 0.;
		}
		if(Model==M3) {
			thmle= NaN;
			alphamle= 0.;
			betamle = betapmle = 0.;
		}
	}

	const int  reflect = copysign( 1, model_in );

	if (output) {
		Rcout << endl << _("maximum-likelihood estimates of parameters:") << endl << endl;
		Rcout << setw(20) << "         theta =" << setw(12) << thmle*reflect << "     ( " << _("x-coordinate of changepoint") << " )" << endl;
		Rcout << setw(20) << "         alpha =" << setw(12) << alphamle << "     ( " << _("y-coordinate of changepoint") << " )" << endl;
		Rcout << setw(20) << "          beta =" << setw(12) << betamle << "     ( " << _("slope of first line") << " )" << endl;
		Rcout << setw(20) << "    beta-prime =" << setw(12) << betapmle << "     ( " << _("slope of second line") << " )" << endl;
		Rcout << setw(20) << "      variance =" << setw(12) << vmle << endl << endl;
	}


	if( par != NULL )  *par= thmle*reflect,  *(par+1)= alphamle,  *(par+2)= betamle,  *(par+3)= betapmle;

	return thmle;
}

