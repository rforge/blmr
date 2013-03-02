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

	MODEL model;
	int n, m, ns;
	bool variance_unknown;

	bool cov_matrix_I, cov_matrix_diagonal, th0ex;
	int k1, gk0, subints;
	double s11, sx1, sxx, n1, se1sq, y1, yx, sysq, qysq, omega, c, th0, alpha0, z, w;
	double star_z, c1, c2, Lgamma, z_tilde_M2, lambdasq, lambda;

	double ah, old_th, prev_th, a_low, a_high, rel_print_eps;
	double SL, prev_SL, cFex, cCHIex, cF, cCHI, x_vu_ex, x_vk_ex, x_vu, x_vk;
	double acc_sl_abs, acc_sl_rel, acc_xb, acc_yb, inc_x, inc_y, aex;

	double *C;
	double *Xd, *Yd, *Sig;
	int *is;
	double *Xs;
	double *q11, *qx1, *qxx, *ck, *qff;
	double *q10, *qx0, *a0, *b0;

	Vector<double> *px, *py, *psy, *pqy;
	Vector<double> *pmq1, *pq1, *pqx, *ps1, *psx, *ps1c, *psxc;
	Vector<double> *psig1, *psigx, *pv1h, *pxh;
	Vector<double> *vund, *vund2, *pnse1, *pnuse1, *puqe1, *puqen, *puqx, *pusen;

	Matrix<double> *prS, *pirS, *pQ;  


//function prototypes
	void initialize(const double *Xdata, const double *Ydata, const double *Sigma =0);
	void get_Q(void) const;
	void get_rS_irS(void) const;
	void get_M(Matrix<double> *pM) const;
	void set_theta0(double th_0, METHOD mode =INIT);
	void set_alpha0(double a_0, METHOD mode =INIT);
	void generate_S(double *rS, int numel =0);
	double max_mk() const;
	bool m_gt_w(double wsq, Vector<double> &s);
	double unif_rand(void) const;
	double rnorm(void) const;

// gamma and f functions
	Vector<double>  gam(double th, int k) const;
	Vector<double>  gfr(double th, int k) const;
	Vector<double>  gsm(double th, int k) const;
	Vector<double>  gbar(double th, int k) const;
	Vector<double>  q_f(double th, int k) const;
	Vector<double>  sf(double th, int k) const;
	Vector<double>  sfc(double th, int k) const;
	double  ff(double th, int k) const;

// sl computing functions
	double sl_af(int mode =0);
	double sl_af2(void);
	double sl_geo(double *err=0);
	double sl_geo1(double *err=0);
	double sl_geo2(double *err=0);
	double sl_mc(void);
	double sl_mc2(void);
	int ci_geo( METHOD met, double *bounds);
	int ci_af( METHOD met, double *bounds );
	double a_sl(METHOD met, double th, int high_low);
	double a_af(double th, int high_low);
	double ahigh(METHOD met, double th);
	double sl_a(double alpha, int k);

//in file geo
	double geo(double th2, double *err);
	double geo_ex(double *err);
	double geo_vu_ex(void);
	double geo_vk_ex(double *err);
	double ipr_ex(double s, int k);
	double geo_vu_D(const double th2, double *err);
	double geo_vu_ND(const double th2, double *err);
	double geo_vu_covNDab(int k, double a, double b, int lo, double *er);
	double geo_vk_D(const double th2, double *err);
	double geo_vk_ND(const double th2, double *err);
	double geo_vk_covNDab(int k, double a, double b, int lo, double *err);
	double ipr(double xi, double *err);

//in file rho_etc:
	double mu_by_tau(double th, int k);
	double mu_by_tau_sq(double th, int k);
	double rho_inv(double s, int k);
	double rho(double th);
	double rho(double theta, int k);
	double rhosq(double th, int k);
	double drho(double th, int k);
	double drhosq(double th, int k);
	double dgsq(double th, int k);
	double F(int k, double arg);
	double fk(int k, double arg);
	double get_C(const int k) const;


//in file Emu_etc
	double fk_arg_invsq(int k, double arg);
	double Emupr(double th, int k);
	double Emupr_vk(double th, int k);
	double bisect(double x1, double x2, double (Cblmr::*fn)(double,int), int k, double value, double crit);
	double bisect_sl(double x1, double x2, METHOD mode, double crit);
	double integrate(double x1, double x2, double (Cblmr::*fn)(double, int), int k, double *err);
	double integrate2(double x1, double x2, double (Cblmr::*fn)(double, double*), double *err2);
	double qgk21(double a, double b, double (Cblmr::*fn)(double, double*), double *err2);
	double pt_m_argsq(double x, double n);
	double F_arg_neginvsq(int k, double arg);

//private versions of interface functions
	void set_Sigma(const double *Sigma_matrix =0);
	void set_x(const double *Xdata);
	void set_xs(void);
	void set_y(const double *Ydata);
	void set_sy(double *irSy, METHOD met =INIT);
	void set_SL( double cSL =0.05);
	void set_acc( double acc =0.001);
	double mle( bool output =true ) const;
	double sl(double theta0, METHOD met = GEO, bool output =true);
	double sl(double theta0, double alpha0, METHOD met = GEO, bool output =true);
	int ci(METHOD met =GEO, bool output =true, double *bounds =0);
	int cr(METHOD met =GEO, double increments =0.2, bool output =true, double *bounds =0);


public:
// constructors
	Cblmr(NumericVector Xda, NumericVector Yda, int n_model =1, bool var_known =false);
	Cblmr(NumericVector Xda, NumericVector Yda, NumericMatrix S, int n_model =1, bool var_known =false);
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
	void MLE(void);
	void SET_irSy(NumericVector irSy);

};





#endif

