//


#include "lmbr.h"




void Clmbr::set_Sigma( void ) 
// 'Sigma' = covariate or covariance vector or matrix,  errors ~ N( 0, var * Sigma )
// check the input 'weights' vector or matrix,  make symetric if almost symetric, 
// then get 'rS' = square root of Sigma, and 'irS' = inverse square root of Sigma 
{
	int i,j;

	if( vectorS )  {
		for (i=0;i<n;i++) {
			const double  wi = *(w_in + i);
			if ( isinf(wi) || isnan(wi) )  stop( _("'weights' has invalid entries") );
			if ( wi <= 0 )  stop( _("'weights' has invalid entries") );
		}

	}  else  {
		Array2D<double>  W2D(n,n);
		for (i=0;i<n;i++) for (j=0;j<n;j++) {
			double  wij = *(w_in+j*n+i),  wji = *(w_in+i*n+j);
			if ( isinf(wij) || isnan(wij) )  stop( _("'weights' has invalid entries") );
			if( fabs(wij) < zero_eq )  wij= 0.; 
			if( fabs(wij - wji) < zero_eq )  if( i < j)  wij = wji;
			W2D[i][j] = wij;
		}

		Cholesky<double> testw(W2D);
		if ( !testw.is_spd() )  stop( _("'weights' matrix is not symetric positive-definite") );
	}


	double  maxD =0.,  minD =Inf;

	if ( vectorS )  {

		for (i=0;i<n;i++) {
			double  Di;
			if( model_in > 0 )  Di= *( w_in + i );  else  Di= *( w_in + (n-1-i) );
			const double  rDi = sqrt( Di );
			if( inverse )
				{  *( rS + i ) = rDi;  *( irS + i ) = 1./rDi; }
			else
				{ *( irS + i ) = rDi;   *( rS + i ) = 1./rDi; }
			if( Di > maxD )  maxD= Di;
			if( Di < minD )  minD= Di;
		}


	}  else  {

// geteigenvalues and eigenvectors of Sigma using LAPACK routine DSYEVR
		double  D[n], Q_[n*n];
		
		{
			const char  job ='V',  range ='A',  uplo ='L';
			const int  id =0;	
			const double  tol = 0, dd = 0;
			int  ne, isuppZ[2*n], lwork= -1, itmp[1], liwork =-1, info =0;
			double  tmp[1];

//  use 'irS' as the working input matrix, so copy 'w_in' to 'irS'
			for(i=0;i<n;i++) for(j=0;j<=i;j++)  
				if( model_in > 0 )  
					*(irS + j*n + i) = *(w_in + j*n + i);  
				else
					*(irS + j*n + i) = *(w_in + (n-1-j)*n + n-1-i);

			F77_CALL(dsyevr)( &job, &range, &uplo, &n, irS, &n, &dd, &dd, &id, &id, &tol,
								&ne, D, Q_, &n, isuppZ, tmp, &lwork, itmp, &liwork, &info );

			if( info )  stop( "LAPACK routine 'dsyevr' failed" );  else  { lwork= *tmp; liwork= *itmp; }
			double  work[lwork];
			int  iwork[liwork];

			F77_CALL(dsyevr)( &job, &range, &uplo, &n, irS, &n, &dd, &dd, &id, &id, &tol,
								&ne, D, Q_, &n, isuppZ, work, &lwork, iwork, &liwork, &info );

			if( info || ne < n )  stop( "LAPACK routine 'dsyevr' failed" );
		}

		double  rD[n];
		for (i=0;i<n;i++) {
			if( D[i] > maxD )  maxD= D[i];
			if( D[i] < minD )  minD= D[i];
			rD[i] =  sqrt( D[i] );
		}

// compute  'rS' = Q*sqrt(D)*QT  and  'irS' = Q*1/sqrt(D)*QT
		for (i=0;i<n;i++)  for(j=0;j<n;j++)  {
			*(rS + j*n + i) = 0.;
			*(irS + j*n + i) = 0.;
			for(int k=0;k<n;k++)  
				if( inverse )  {
					 *(rS + j*n + i) +=  *(Q_+k*n+i) * rD[k] * ( *(Q_+k*n+j) );
					*(irS + j*n + i) +=  *(Q_+k*n+i) / rD[k] * ( *(Q_+k*n+j) );
				}  else  {
					 *(rS + j*n + i) +=  *(Q_+k*n+i) / rD[k] * ( *(Q_+k*n+j) );
					*(irS + j*n + i) +=  *(Q_+k*n+i) * rD[k] * ( *(Q_+k*n+j) );
				}
		}

	}


	if( minD/maxD < 1.e-7 )  Rf_warning( _("weights matrix ill-conditioned for 'clr' method") );


	if ( px != NULL )  set_x();

	return;
}




