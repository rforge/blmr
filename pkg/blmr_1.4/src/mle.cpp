//

#include "blmr.h"



double Cblmr::mle( const bool output ) const
// find theta_mle for maximum value of (gamma(theta)*qy)^2
// calculate  alphamle, betapmle, betamle, vmle
{
	bool line= false;
	if( qysq/m < zero_eq )  line= true;


// find theta_mle
	double  vf= (x_d[n_d-1]-x_d[n_d-2])*(*pqy*pq1[n_d-1]), max= vf*vf/qff[n_d-1], thmle= x_d[n_d-2];
	int kmle=n_d-2;

	for (int k=n_d-3;k>=k1;k--) {

		const double  v1 = (*pqy)*pq1[k+1],  vx = vf + v1*x_d[k+1];  vf = vx - v1*x_d[k];

		const double  va = vx*qx1[k+1] - v1*qxx[k+1], vb = vx*q11[k+1] - v1*qx1[k+1], thk = va/vb;

		double mk;
		if(x_d[k]<thk && thk<x_d[k+1]) {
			mk = (vx*vb - v1*va)/ck[k+1]; 
			if (mk > max) {max = mk; thmle = thk; kmle=k+1;}
		} else { 
			mk = vf*vf/qff[k+1];
			if (mk > max) {max = mk; thmle = x_d[k]; kmle=k+1;}
		}
	}
	if ( fabs(x_d[kmle-1] - thmle) < zero_eq )  thmle = x_d[kmle-1];
	if ( fabs(x_d[kmle] - thmle) < zero_eq )  {  kmle++; thmle = x_d[kmle-1]; }
	if ( line )  thmle= x_d[n_d-1];


// calculate  alphamle, betapmle, betamle, vmle
	double alphamle, betapmle, betamle, vmle;

	if (Model==M1) {
		const double  amle = *psy*gbar(thmle,kmle);
		const double bpmle = *psy*gsm(thmle,kmle);
		const double  bmle = *psy*gfr(thmle,kmle);

		double psi;
		if (cov_matrix_diagonal) psi=0.; else psi = gfr(thmle,kmle)*sfc(thmle,kmle);
		const double sm1 = *psig1*gsm(thmle,kmle);
		const double fr1 = *psig1*gfr(thmle,kmle);
		const double r1sf = sqrt(s11-sm1*sm1-fr1*fr1);
		alphamle = amle/r1sf;
		betapmle = (bpmle - alphamle*sm1)/sqrt(sfc(thmle,kmle)*sfc(thmle,kmle)-psi*psi);
		betamle = (bmle - alphamle*fr1 - betapmle*psi)/sqrt(sf(thmle,kmle)*sf(thmle,kmle));

	}  else  {

		const double sf1 = *pv1h*sf(thmle,kmle);
		const double y1 = *pv1h*(*psy);

		betapmle = 0.;
		betamle = ( *psy*sf(thmle,kmle) - sf1*y1 )/( sf(thmle,kmle)*sf(thmle,kmle) - sf1*sf1 );
		alphamle = (y1 - betamle*sf1)/n1;
	}

	if (variance_unknown)  vmle = omega/(m-2);  else  vmle = 1.;

	if( fabs(thmle) < zero_eq) thmle= 0.;
	if( fabs(alphamle) < zero_eq) alphamle= 0.;
	if( fabs(betapmle) < zero_eq) betapmle= 0.;
	if( fabs(betamle) < zero_eq) betamle= 0.;
	if( fabs(vmle) < zero_eq) vmle= 0.;

	if(line) {
		if(Model==M1) {
			thmle= und;
			alphamle= und;
			const double slope= ( (*py)[ i_d[1] ] - (*py)[ i_d[0] ] )/( x_d[1]-x_d[0]);
			betapmle = betamle = slope;
		}
		if(Model==M2) {
			thmle= und;
			alphamle= (*py)[0];
			betapmle = betamle = 0.;
		}
	}

	if (output) {
		Rcout << _("maximum likelihood estimates of parameters") << endl;
		Rcout << setw(10) << "theta" << setw(13) << "alpha" << setw(18) << "beta-prime" << 
			setw(12) << "beta" << setw(13) << "var" << endl;
		if( model_in != -2 )
			Rcout << setw(10) << thmle << setw(14) << alphamle << setw(16) << betapmle << 
				setw(14) << betamle << setw(14) << vmle << endl << endl;
		else
			Rcout << setw(10) << -thmle << setw(14) << alphamle << setw(16) << betamle << 
				setw(14) << betapmle << setw(14) << vmle << endl << endl;
	}


	return thmle;
}

