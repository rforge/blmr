//
//  define class  Cblmr  for
//  data-and-model objects 
//  which are applications of the broken line model to observed data sets
//


#if !defined Cblmr_h			//prevents compiler from repeating this code in other files
#define Cblmr_h

#include "globals.h"


class Cblmr
{

private:

	MODEL Model;
	int model_in, n, m, n_d;
	bool variance_unknown;

	bool cov_matrix_I, cov_matrix_diagonal, th0ex, trivial;
	int k1, k0, subints;
	double s11, sx1, sxx, n1, se1sq, y1, yx, sysq, qysq, omega, c, th0, alpha0, z, w;
	double prime_z, c1, c2, Lgamma, lambdasq, lambda;

	double ah, old_th, prev_th, a_low, a_high, rel_print_eps;
	double SL, prev_SL, cFex, cCHIex, cF, cCHI, x_vu_ex, x_vk_ex, x_vu, x_vk;
	double acc_rho, acc_sl_abs, acc_sl_rel, acc_xb, acc_yb, inc_x, inc_y;

	double *X_in, *Y_in, *S_in;
	int *i_d;
	double *x_d;
	double *q11, *qx1, *qxx, *ck, *qff;
	double *q10, *qx0, *a0, *b0;
	double *C;

	Vector<double> *px, *psig1, *psigx, *pv1h, *pxh;
	Vector<double> *und_n, *pnse1, *pnuse1, *pusen;
	Vector<double> *und_m, *puqe1, *puqen, *puqx;
	Vector<double> *ps1, *psx, *ps1c, *psxc;
	Vector<double> *pq1, *pqx;
	Vector<double> *pmq1; 

	Vector<double> *py, *psy;
	Vector<double> *pqy;

	Matrix<double> *prS, *pirS, *pQ;  



// function prototypes:

//initializing functions
	void initialize( void );
	void get_Q(void) const;
	void get_rS_irS(void) const;
	void get_M(Matrix<double> *pM) const;
	void set_theta0(double th_0, METHOD met =INIT);
	void set_alpha0(double a_0, METHOD met =INIT);
	void set_Sigma( void );
	void set_x(void);
	void set_y(void);
	void set_SL( double cSL =0.05);
	void set_acc( double acc =0.001);
	double max_mk() const;
	bool m_gt_w(double wsq, const Vector<double> &s) const;

// gamma and f functions
	Vector<double>  gam(double th, int data_interval) const;
	Vector<double>  gfr(double th, int data_interval) const;
	Vector<double>  gsm(double th, int data_interval) const;
	Vector<double>  gbar(double th, int data_interval) const;
	Vector<double>  gbar_prime(double th, int data_interval) const;
	Vector<double>  q_f(double th, int data_interval) const;
	Vector<double>  sf(double th, int data_interval) const;
	Vector<double>  sfc(double th, int data_interval) const;
	double  ff(double th, int data_interval) const;

//in file rho_etc:
	double rho(double th) const;
	double rho(double theta, int data_interval) const;
	double rhosq(double th, int data_interval) const;
	double drho(double th, int data_interval) const;
	double drhosq(double th, int data_interval) const;
	double dgsq(double th, int data_interval) const;
	double rho_inv(double s, int data_interval, int hi_lo =1) const;

// sl computing functions
	double sl_af(int mode =0) const;
	double sl_af2(void) const;
	double sl_geo(double *err=0);
	double sl_geo1(double *err=0);
	double sl_geo2(double *err=0);
	double ipr(double xi, double *err);
	double sl_mc(void) const;
	double sl_mc2(void) const;

//in files geo_ex and geo
	double geo_ex(void) const;
	double geo_vu_ex(void) const;
	double geo_vk_ex(void) const;
	double geo(double th_2, double *err) const;
	double geo_vu_D(double th2, double *err) const;
	double geo_vu_ND(double th2, double *err) const;
	double geo_vu_covNDab(int data_interval, double a, double b, int hi_lo, double *err) const;
	double geo_vk_D(double th2, double *err) const;
	double geo_vk_ND(double th2, double *err) const;
	double geo_vk_covNDab(int data_interval, double a, double b, int hi_lo, double *err) const;

//in file Fm_fm
	double F(int k, double arg) const;
	double fk(int k, double arg) const;
	double get_C(int k) const;
	double sF(int k, double arg) const;

//in file Emu_etc
	double mu_by_tau(double th, int data_interval) const;
	double amu_by_Omega(double th, int data_interval) const;
	double Emupr(double th, int data_interval) const;
	double Emupr_vk(double th, int data_interval) const;

//in file 'bis_int' bisection and integration routines
	double bisect(double x1, double x2, double (Cblmr::*fn)(double,int), int k, double value, double crit);
	double bisect(double x1, double x2, double (Cblmr::*fn)(double,int) const, int k, double value, double crit) const;
	double bisect_sl(double x1, double x2, METHOD met, double crit);
	double integrate(double x1, double x2, double (Cblmr::*fn)(double, int) const, int k, double *err) const;
	double integrate2(double x1, double x2, double (Cblmr::*fn)(double, double*), double *err);
	double qgk21(double a, double b, double (Cblmr::*fn)(double, double*), double *err);

//private versions of interface functions
	double sl(double theta0, METHOD met = GEO, bool output =true);
	double sl(double theta0, double alpha0, METHOD met = GEO, bool output =true);
	int ci(METHOD met =GEO, double increments =0.2, bool output =true, double *bounds =0);
	int cr(METHOD met =GEO, double increments =0.2, bool output =true, double *bounds =0);
	double mle( bool output =true ) const;
	void set_sy(double *irSy, METHOD met =INIT);

// ci and cr sub-functions
	int ci_geo( METHOD met, double increments =0.2, double *bounds =0);
	int ci_af( METHOD met, double *bounds =0);
	double a_sl(METHOD met, double th, int high_low);
	double a_af(double th, int high_low);
	double ahigh(METHOD met, double th);
	double sl_a(double alpha, int k);


public:
// constructors
	Cblmr(NumericVector yR, NumericVector xR, int model =1, bool var_known =false, NumericMatrix SigmaR =NULL);
	Cblmr(const Cblmr &initM);	// copy constructor
	~Cblmr();		// destructor


// interface functions
	void slR(double theta0);
	void slR(double theta0, double alpha0);
	void slR(double theta0, string met );
	void slR(double theta0, double alpha0, string met );
	void slR(double theta0, string met, double acc);
	void slR(double theta0, double alpha0, string met, double acc);
	double slR(double theta0, string met, double acc, bool output);
	double slR(double theta0, double alpha0, string met, double acc, bool output ); 
	void ciR(double CL);
	void ciR(double CL, string met);
	void ciR(void);
	void crR(double CL);
	void crR(double CL, string met);
	void crR(double CL, string met, double incr);
	void crR(void);
	void MLE(void) const;
	void SET_irSy(NumericVector irSy);

};




#endif

