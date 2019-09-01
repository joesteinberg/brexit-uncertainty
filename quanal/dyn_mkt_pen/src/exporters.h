#ifndef __EXPORTERS_H__
#define __EXPORTERS_H__

#include "globals.h"

///////////////////////////////////////////////////////////////////////////////
// Key structures and simple utilities
///////////////////////////////////////////////////////////////////////////////

int calc_all_moments;

typedef struct
{
  double kappa; // tail parameter for productivity distribution
  double alpha; // returns to scale in marketing
  double beta; // diminishing returns in marketing
  double psi; // probability of converting new customers
  double delta; // probability of losing customers who don't see ads
  double xi; // non-death exit rate
  double omega; // death rate
  double theta; // elasticity

  double Zi;
  double phi_grid[NPHI]; // productivity grid
  double phi_probs[NPHI]; // probabilities
  double phi_cumprobs[NPHI]; // cumulative probabilities
  double phi_pow_grid[NPHI]; // phi^{sigma-1}
}exporter_params;

typedef struct
{
  // prices and market size
  double W;
  double MC;
  double B;
  double Q;
  double L;
  double La;
  double tau;
  double Bii;
  double Lii;
  double psi_mult;

  // cross-sectional moments
  double exports;
  double export_participation_rate;
  double theoretical_export_participation_rate;
  double top5_share;
  double avg_rel_entrant_size;
  double lf;
  double exit_rate;
  double entrant_exit_rate;
  double rel_entrant_growth_rate;
  double startup_cost;
  double continuation_cost;
  double startup_cost2;
  double continuation_cost2;

  // destination-specific equilibrium objects
  double cutoff;
  double n_grid[NPHI][NN]; // customer base grid
  double hne[NPHI]; // marketing policy function for entrants
  double hn[NPHI][NN]; // marketing policy function for incumbents
  double diste[NPHI]; // distribution for entrants
  double dist[NPHI][NN]; // distribution for incumbents
  double Z;
  double n;
  double n2;

  // full mkt pen version
  double vn[NPHI]; // distribution for incumbents (full mkt pen version)
  double vx[NPHI]; // distribution for incumbents (full mkt pen version)
  double disti[NPHI]; // distribution for incumbents (full mkt pen version)

}exporter_vars;

typedef struct
{
  int ip;
  const exporter_params * ep;
  exporter_vars * ev;
  exporter_vars * evp;
  gsl_interp_accel * acc;
}foc_params;

typedef struct
{
  int ip;
  const exporter_params * ep;
  exporter_vars * ev;
  exporter_vars * evp1;
  exporter_vars * evp2;
  double prob;
  gsl_interp_accel * acc1;
  gsl_interp_accel * acc2;
}foc_stoch_params;

typedef struct
{
  int ip;
  const exporter_params * ep;
  exporter_vars * ev;
  exporter_vars ** evplist;
  double * probs;
  gsl_interp_accel ** acclist;
}foc_stoch_idio_params;

typedef struct
{
  double v;
  double w;
  double cw;
}export_obs;

void copy_exporter_params(exporter_params * dest, const exporter_params * src);
void copy_exporter_vars(exporter_vars * dest, const exporter_vars * src);

void set_exporter_consts1(exporter_vars * ev, double W, double Q, double MC, double psi_mult);
void set_exporter_consts2(exporter_vars * ev, double L, double P2, double Y2, double Ybar2, double theta, double alpha, double tau);
void set_exporter_consts3(exporter_vars * ev, double Lii, double P2ii, double Y2ii, double Ybarii2, double theta);

void discretize_phi(exporter_params * ep);

double mp_f(const exporter_params * ep, double np, double n, double L);
double mp_f1(const exporter_params * ep, double np, double n, double L);
double mp_f2(const exporter_params * ep, double np, double n, double L);

int set_entry_cutoff(const exporter_params * ep, exporter_vars * ev);
int check_cutoff(const exporter_params * ep, exporter_vars  * ev);
void init_policy(exporter_vars * ev);

double interp(gsl_interp_accel * acc, const double *xa, const double *ya, int n, double x);

double entrant_foc(double np, void * params);
double entrant_foc_stoch(double np, void * params);
double entrant_foc_stoch_idio(double np, void * params);

double invert_foc(const exporter_params * ep,
		  exporter_vars * ev,
		  exporter_vars * evp,
		  int ip, double np,
		  double npp);

double invert_foc_stoch(const exporter_params * ep,
			exporter_vars * ev,
			exporter_vars * evp1,
			exporter_vars * evp2,
			int ip,
			double np,
			double npp1,
			double npp2,
			double prob);

double invert_foc_stoch_idio(const exporter_params * ep,
			     exporter_vars * ev,
			     exporter_vars * evplist[6],
			     int ip,
			     double np,
			     double npplist[6],
			     double probs[6]);

int find_root_1d(gsl_function * f, double xlo, double xhi, double * x);

int iterate_entrant_policy_fn(const exporter_params * ep,
			      exporter_vars * ev,
			      exporter_vars * evp,
			      double * supnorm);

int iterate_entrant_policy_fn_stoch(const exporter_params * ep,
				    exporter_vars * ev,
				    exporter_vars * evp1,
				    exporter_vars * evp2,
				    double prob);

int iterate_entrant_policy_fn_stoch_idio(const exporter_params * ep,
					 exporter_vars * ev,
					 exporter_vars * evplist[6],
					 double probs[6]);

int iterate_incumbent_policy_fn(const exporter_params * ep,
				exporter_vars * ev,
				exporter_vars * evp,
				double * supnorm);

int iterate_incumbent_policy_fn_stoch(const exporter_params * ep,
				      exporter_vars * ev,
				      exporter_vars * evp1,
				      exporter_vars * evp2,
				      double prob);

int iterate_incumbent_policy_fn_stoch_idio(const exporter_params * ep,
					   exporter_vars * ev,
					   exporter_vars * evplist[6],
					   double probs[6]);

int iterate_vf_policy_full_mkt_pen(const exporter_params * ep,
				   exporter_vars * ev,
				   exporter_vars * evp,
				   double * supnorm);

int iterate_vf_policy_full_mkt_pen_stoch(const exporter_params * ep,
					 exporter_vars * ev,
					 exporter_vars * evp1,
					 exporter_vars * evp2,
					 double prob);

int iterate_vf_policy_full_mkt_pen_stoch_idio(const exporter_params * ep,
					      exporter_vars * ev,
					      exporter_vars * evplist[6],
					      double probs[6]);

int solve_steady_state_policies(const exporter_params * ep,
				exporter_vars * ev);

void init_dist(const exporter_params * ep,
	       exporter_vars * ev);

int update_dist(const exporter_params * ep,
		exporter_vars * ev,
		exporter_vars * evp,
		double diste[NPHI],
		double dist[NPHI][NN],
		double * supnorm);

int update_dist_full_mkt_pen(const exporter_params * ep,
			     exporter_vars * ev,
			     exporter_vars * evp,
			     double diste[NPHI],
			     double disti[NPHI],
			     double * supnorm);

int solve_steady_state_dist(const exporter_params * ep, exporter_vars * ev);
int compute_moments(const exporter_params * ep, exporter_vars * ev);
int compute_moments_full_mkt_pen(const exporter_params * ep, exporter_vars * ev);
int compute_moments_static_version(const exporter_params * ep, exporter_vars * ev);
int convert_psi(exporter_params * ep, exporter_vars * ev);

#endif
