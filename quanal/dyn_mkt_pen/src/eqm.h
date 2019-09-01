#ifndef __EQM_H__
#define __EQM_H__

#include "globals.h"
#include "calibrate.h"
#include "exporters.h"
#include "solver.h"

extern const uint nbgp; // number of BGP variables... used to solve for BGP to create initial guess for equilibrium
uint neqm; // number of equilibrium vars/eqns... depends on situation
double bbgp[NC-1]; // BGP bond-holdings... state variable
uint scenario; // scenario... may be redundant in stochastic model

// eqm contains all the vars associated with a deterministic equilibrium, or a single history of a stochastic equilibrium
typedef struct
{
  double Q_t[NT+1][NC];

  double B_t[NT+1][NC];
  double C_t[NT+1][NC];
  double MUC_t[NT+1][NC];
  double Chabit_t[NT+1][NC];
  double X_t[NT+1][NC];
  double Lf_t[NT+1][NC];
  double L_t[NT+1][NC];
  double K_t[NT+1][NC];
  double Y_t[NT+1][NC];
  double M_t[NT+1][NC];
  double Pi_t[NT+1][NC];
  double VA_t[NT+1][NC];
  double T_t[NT+1][NC];
  double Inc_t[NT+1][NC];
  double Exp_t[NT+1][NC];

  double P_t[NT+1][NC];
  double W_t[NT+1][NC];
  double R_t[NT+1][NC];
  double Lambda_t[NT+1][NC];
  double MC_t[NT+1][NC];

  double NGDP_t[NT+1][NC];
  double RGDP_t[NT+1][NC];
  double XY_t[NT+1][NC];

  double EX_t[NT+1][NC][NC];
  double IM_t[NT+1][NC][NC];
  double EXR_t[NT+1][NC][NC];
  double IMR_t[NT+1][NC][NC];
  double NX_t[NT+1][NC][NC];
  double RER_t[NT+1][NC][NC];
  double TOT_t[NT+1][NC][NC];

  double P2_t[NT+1][NC][NC];
  double Y2_t[NT+1][NC][NC];
  double Y2s_t[NT+1][NC][NC];
  double Kd_t[NT+1][NC];

  double D_t[NT+1][NC][NC];
  double Db_t[NT+1][NC][NC];
  double Dh_t[NT+1][NC][NC];
  double Dt_t[NT+1][NC][NC];

  double exrate_t[NT+1][NC][NC-1];
  double texrate_t[NT+1][NC][NC-1];
  double mkt_pen_t[NT+1][NC][NC-1];
  double mkt_pen_w_t[NT+1][NC][NC-1];
  exporter_vars ev_t[NT+1][NC][NC-1];
  exporter_vars ev_good_t[NT+1][2];
  exporter_vars ev_bad_t[NT+1][2];

  double welfare_t[NT+1][NC];
  double welfare2_t[NT+1][NC];
  double welfare3_t[NT+1][NC];
  double welfare4_t[NT+1][NC];
  double welfare_cost_t[NT+1][NC];  
}eqm;

// array of NTH eqm structs for use in solving deterministic equilibrium (like no-Brexit counterfactual)
// we use the array to parallelize the process of evaluating the jacobian matrix in the solver
eqm eee0[NTH];
eqm eee1[NTH];
eqm eee2[NTH];

// stochastic equilibrium structure... contains NHIST or NHIST_CH eqm structs (one for each  possible history)
typedef struct
{
  eqm eee[NHIST];
}stoch_eqm;

// stochastic equilibrium structure... contains NHIST or NHIST_CH eqm structs (one for each  possible history)
typedef struct
{
  eqm eee[NHIST_RB];
}stoch_eqm_rb;

// array of NTH stoch_eqm structs... serves same parallelization purpose as eee0[NTH]
stoch_eqm sss[NTH];

stoch_eqm_rb sss_rb[NTH];
eqm eee3[NTH];
eqm eee4[NTH];

void set_neqm();
void init_vars(eqm * e);
void copy_vars(eqm * e1, const eqm * e0);

uint stack_bgp_vars(double * myx, const eqm * e);
uint unstack_bgp_vars(eqm * e, const double * myx);
uint stack_eqm_vars(double * myx, const eqm * e);
uint unstack_eqm_vars(eqm * e, const double * myx);
uint stack_stoch_eqm_vars(double * myx, const stoch_eqm * s);
uint unstack_stoch_eqm_vars(stoch_eqm * s, const double * myx);
uint stack_stoch_eqm_rb_vars(double * myx, const stoch_eqm_rb * s);
uint unstack_stoch_eqm_rb_vars(stoch_eqm_rb * s, const double * myx);

uint set_initial_eqm_guess();
uint set_initial_bpg_guess();

uint write_bgp_vars(const eqm * e, const char * fname);
uint write_eqm_vars(const params * p, const eqm * e, const char * fname, uint i);

uint set_vars1(eqm * e, const params * p, uint t, uint bgp);
uint set_vars2a(eqm * e, const params * p, uint t, uint bgp);
uint set_stoch_vars2a(eqm * e, eqm * e_good, eqm * e_bad, const params * p, uint t, double pi);
uint set_stoch_vars2a_idio(eqm * e, eqm * e_good, eqm * e_bad, const params * p, uint t, double pi);
uint set_vars2b(eqm * e, const params * p, uint t, uint bgp);
uint set_vars3(eqm * e, const params * p, uint t, uint bgp);

uint eval_bgp_conds(const double * myx, double * myf, uint tn);
uint solve_bgp(double bb[NC-1]);

uint eval_eqm_conds(const double * myx, double * myf, uint tn);
uint solve_eqm();

uint eval_stoch_eqm_conds(const double * myx, double * myf, uint tn);
uint solve_stoch_eqm();

uint eval_stoch_eqm_rb_conds(const double * myx, double * myf, uint tn);
uint solve_stoch_eqm_rb();

void calc_welfare(eqm * e, const params * p);
//void calc_stoch_welfare();
//void calc_stoch_welfare_rb();
//void calc_ce_welfare1();
//void calc_ce_welfare2();

static inline double muc_mul(const params * p, const eqm * e, uint t, uint i)
{
  if(fixl==1 || (fixl==2 && (t<TREF || scenario == 0)))
    {
      return (e->L_t[t][i] - p->Lbar[i]/3.0)/(p->Lbar[i]/3.0);
    }
  else if(fixl==2 && t>=TREF && scenario > 0)
    {
      return e->W_t[t][i] - e->W_t[TREF-1][i];
    }
  else
    {
      return 1.0 * (e->MUC_t[t][i]/e->P_t[t][i] - 
	mul(
	    e->C_t[t][i],
	    e->L_t[t][i],
	    p->Lbar[i],
	    p->phi[i],
	    p->psi)/e->W_t[t][i]);
    }
}

static inline double euler(const params * p, const eqm * e, uint t, uint i)
{
  double retval = 0.0;
  if(habit_formation)
    {
      retval=100.0 * (e->Q_t[t][i] * e->MUC_t[t][i]
		      - p->beta * e->P_t[t+1][0] * (e->P_t[t][i]/e->P_t[t+1][i])
		      * e->MUC_t[t+1][i]);
    }
  else
    {
      retval=10000.0 * (e->Q_t[t][i] * e->MUC_t[t][i]
		      - p->beta * e->P_t[t+1][0] * (e->P_t[t][i]/e->P_t[t+1][i])
		      * e->MUC_t[t+1][i]);
    }
  return retval;
}

static inline double bop(const params * p, const eqm * e, uint t, uint i)
{
  uint j;
  double net_iceberg = 0.0;
  for(j=0; j<NC; j++)
    {
      if(j!=i)
	{
	  net_iceberg += p->ntb_ts[t][j][i]*e->P2_t[t][j][i]*e->Y2_t[t][j][i];
	  net_iceberg -= p->ntb_ts[t][i][j]*e->P2_t[t][i][j]*e->Y2_t[t][i][j];
	}
    }

  if(t<NT)
    {
      return SUM(e->NX_t[t][i],NC) + net_iceberg + 
	e->B_t[t][i]*e->P_t[t][0] - e->Q_t[t][i]*e->B_t[t+1][i];
    }
  else
    {
      return SUM(e->NX_t[t][i],NC) + net_iceberg + 
	e->B_t[t][i]*e->P_t[t][0] - e->Q_t[t][i]*e->B_t[t][i];
    }
}

static inline double mkt_clear_Y2(const eqm * e, uint t, uint i, uint j)
{
  return (e->Y2_t[t][i][j] - e->Y2s_t[t][i][j])/e->Y2s_t[t][i][j];
}

static inline double mkt_clear_K(const eqm * e, uint t, uint i)
{
  return (e->K_t[t][i] - e->Kd_t[t][i])/e->Kd_t[t][i];
}

static inline double price_norm(const eqm * e, uint t)
{
  return e->P_t[t][0] - 1.0;
}

static inline double fin_mkt(const eqm * e, uint t, uint i)
{
  if(fin_aut==1)
    {
      return e->B_t[t+1][i-1];
    }
  else if(comp_mkts==1 && scenario==1)
    {
      return 10000. * ((e->MUC_t[t][i]/e->P_t[t][i]) /  (e->MUC_t[t][0]/e->P_t[t][0])
		       -  (e->MUC_t[0][i]/e->P_t[0][i]) /  (e->MUC_t[0][0]/e->P_t[0][0]));
    }
  else
    {
      return e->Q_t[t][i] - e->Q_t[t][0];
    }
}

static inline double noarb(const params * p, const eqm * e, uint t, uint i)
{
  double lam = p->beta * (e->MUC_t[t+1][i]/e->MUC_t[t][i]) * (e->P_t[t][i]/e->P_t[t+1][i]);

  if(p->cap_adj_cost==0)
    {
      //e->R_t[t][i] = e->P_t[t][0]*e->P_t[t-1][i]/e->Q_t[t-1] - e->P_t[t][i]*(1.0-p->delta);
      //return e->Q_t[t] - (e->P_t[t+1][0]*e->P_t[t][i]/(e->R_t[t+1][i]+e->P_t[t+1][i]*(1.0-p->delta)));
      return e->P_t[t][i] - lam * (e->R_t[t+1][i] + e->P_t[t+1][i]*(1.0-p->delta));
    }
  else
    {
      double x = e->X_t[t][i]/e->K_t[t][i];
      double xp = e->X_t[t+1][i]/e->K_t[t+1][i];
      double dH = dphiK(x,p->delta,p->etaK);
      double dHp = dphiK(xp,p->delta,p->etaK);
      double Hp = phiK(xp,p->delta,p->etaK);
      double RR = e->R_t[t+1][i] + (e->P_t[t+1][i]/dHp) * (1.0-p->delta + Hp - xp*dHp);

      return e->P_t[t][i]/dH - lam * RR;

      //e->R_t[t][i] = e->P_t[t][0]*e->P_t[t-1][i]/e->Q_t[t-1]/dHm
      //  + (e->P_t[t][i]/dH) * (dH * x - H - (1.0-p->delta));

      //return e->Q_t[t] - (e->P_t[t+1][0]*e->P_t[t][i]/dH/RR);
    }
}

int bgp_func_f(const gsl_vector * x, void * data, gsl_vector * f);
int bgp_func_df(const gsl_vector * x, void * data, gsl_matrix * J);
int bgp_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);
int eqm_func_f(const gsl_vector * x, void * data, gsl_vector * f);
int eqm_func_df(const gsl_vector * x, void * data, gsl_matrix * J);
int eqm_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);
int stoch_eqm_func_f(const gsl_vector * x, void * data, gsl_vector * f);
int stoch_eqm_func_df(const gsl_vector * x, void * data, gsl_matrix * J);
int stoch_eqm_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);
int stoch_eqm_rb_func_f(const gsl_vector * x, void * data, gsl_vector * f);
int stoch_eqm_rb_func_df(const gsl_vector * x, void * data, gsl_matrix * J);
int stoch_eqm_rb_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);

#endif
