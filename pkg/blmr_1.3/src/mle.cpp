//

#include "blmr.h"



double Cblmr::mle( bool output ) const
// find theta_mle for maximum value of (gamma(theta)*qy)^2
// calculate  alphamle, betapmle, betamle, vmle
{

// find theta_mle
	double  vf= (Xs[ns-1]-Xs[ns-2])*(*pqy*pq1[ns-2]), max= vf*vf/qff[ns-2], thmle= Xs[ns-2];
	int kmle=ns-2;

	for (int k=ns-3;k>=k1;k--) {

		double  mk,  v1 = (*pqy)*pq1[k],  vx = vf + v1*Xs[k+1];  vf = vx - v1*Xs[k];

		double va = vx*qx1[k] - v1*qxx[k], vb = vx*q11[k] - v1*qx1[k], thk = va/vb;

		if(Xs[k]<thk && thk<Xs[k+1]) {
			mk =(vx*vb - v1*va)/ck[k]; 
			if (mk > max) {max = mk; thmle = thk; kmle=k;}
		} else { 
			mk =vf*vf/qff[k];
			if (mk > max) {max = mk; thmle = Xs[k]; kmle=k;}
		}
	}


// calculate  alphamle, betapmle, betamle, vmle
	double alphamle, betapmle, betamle, vmle;

	if (model==M1) {
		double  amle = *psy*gbar(thmle,kmle);
		double bpmle = *psy*gsm(thmle,kmle);
		double  bmle = *psy*gfr(thmle,kmle);

		double psi;
		if (cov_matrix_diagonal) psi=0.; else psi = gfr(thmle,kmle)*sfc(thmle,kmle);
		double sm1 = *psig1*gsm(thmle,kmle);
		double fr1 = *psig1*gfr(thmle,kmle);
		double r1sf = sqrt(s11-sm1*sm1-fr1*fr1);
		alphamle = amle/r1sf;
		betapmle = (bpmle - alphamle*sm1)/sqrt(sfc(thmle,kmle)*sfc(thmle,kmle)-psi*psi);
		betamle = (bmle - alphamle*fr1 - betapmle*psi)/sqrt(sf(thmle,kmle)*sf(thmle,kmle));

	}  else  {

		double sf1 = *pv1h*sf(thmle,kmle);
		double y1 = *pv1h*(*psy);

		betapmle = 0.;
		betamle = ( *psy*sf(thmle,kmle) - sf1*y1 )/( sf(thmle,kmle)*sf(thmle,kmle) - sf1*sf1 );
		alphamle = (y1 - betamle*sf1)/n1;
	}

	if (variance_unknown)  vmle = omega/(m-2);  else  vmle = 1.;


	if (output) {
		Rcout << _("maximum likelihood estimates of parameters") << endl;
		Rcout << setw(10) << "theta" << setw(11) << "alpha" << setw(15) << "beta-prime" << 
			setw(9) << "beta" << setw(11) << "var" << endl;
		Rcout << setw(10) << thmle << setw(12) << alphamle << setw(12) << betapmle << 
			setw(12) << betamle << setw(12) << vmle << endl << endl;
	}


	return thmle;
}
