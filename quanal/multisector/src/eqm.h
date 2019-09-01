#ifndef __EQM_H__
#define __EQM_H__

#include "globals.h"
#include "calibrate.h"
#include "solver.h"

extern const uint nbgp; // number of BGP variables... used to solve for BGP to create initial guess for equilibrium
uint neqm; // number of equilibrium vars/eqns... depends on situation
double bbgp[NC-1]; // BGP bond-holdings... state variable
uint scenario; // scenario... may be redundant in stochastic model

// eqm contains all the vars associated with a deterministic equilibrium, or a single history of a stochastic equilibrium
typedef struct
{
  double pb_t[NT+1];

  double b_t[NT+1][NC];
  double cc_t[NT+1][NC];
  double ii_t[NT+1][NC];
  double ll_t[NT+1][NC];
  double kk_t[NT+1][NC];
  double cpi_t[NT+1][NC];
  double pi_t[NT+1][NC];
  double w_t[NT+1][NC];
  double rk_t[NT+1][NC];
  double ngdp_t[NT+1][NC];
  double rgdp_t[NT+1][NC];
  double iy_t[NT+1][NC];
  double lp_agg_t[NT+1][NC];

  double ex_t[NT+1][NC][NC];
  double im_t[NT+1][NC][NC];
  double nx_t[NT+1][NC][NC];
  double exf_t[NT+1][NC][NC];
  double imf_t[NT+1][NC][NC];
  double nxf_t[NT+1][NC][NC];
  double exm_t[NT+1][NC][NC];
  double imm_t[NT+1][NC][NC];
  double nxm_t[NT+1][NC][NC];
  double rer_t[NT+1][NC][NC];

  double rex_t[NT+1][NC][NC];
  double rim_t[NT+1][NC][NC];
  double rexf_t[NT+1][NC][NC];
  double rimf_t[NT+1][NC][NC];
  double rexm_t[NT+1][NC][NC];
  double rimm_t[NT+1][NC][NC];

  double y_t[NT+1][NC][NS];
  double py_t[NT+1][NC][NS];
  double va_t[NT+1][NC][NS];
  double md_t[NT+1][NC][NS][NS];
  double k_t[NT+1][NC][NS];
  double l_t[NT+1][NC][NS];
  double is_t[NT+1][NC][NS];
  double lp_t[NT+1][NC][NS];

  double exs_t[NT+1][NC][NS][NC];
  double ims_t[NT+1][NC][NS][NC];
  double nxs_t[NT+1][NC][NS][NC];
  double rexs_t[NT+1][NC][NS][NC];
  double rims_t[NT+1][NC][NS][NC];

  double pm_t[NT+1][NC][NS];
  double m_t[NT+1][NC][NS];  
  double m2_t[NT+1][NC][NS][NC];

  double p_t[NT+1][NC][NS];
  double q_t[NT+1][NC][NS];  
  double q2_t[NT+1][NC][NS][NC];

  double c_t[NT+1][NC][NS];
  double i_t[NT+1][NC][NS];

  double welfare_t[NT+1][NC];
  double welfare2_t[NT+1][NC];
  double welfare3_t[NT+1][NC];
  double welfare4_t[NT+1][NC];
  double welfare_cost_t[NT+1][NC];

  double Q_t[NT+1][NC];
  
}eqm;

// array of NTH eqm structs for use in solving deterministic equilibrium (like no-Brexit counterfactual)
// we use the array to parallelize the process of evaluating the jacobian matrix in the solver
eqm eee0[NTH];
eqm eee1[NTH];
eqm eee2[NTH];

// stochastic equilibrium structure... contains NHIST eqm structs (one for each  possible history)
typedef struct
{
  //double pi;
  eqm eee[NHIST];
}stoch_eqm;

// array of NTH stoch_eqm structs... serves same parallelization purpose as eee0[NTH]
stoch_eqm sss[NTH];

void set_neqm(); // sets the dimension of the equilibrium solution space
void init_vars(eqm * e); // initializes all the variables of an eqm struct to zero (or 1 where appropriate)
void copy_vars(eqm * e1, const eqm * e0); // copies all the variables from one eqm struct to another
uint stack_bgp_vars(double * myx, const eqm * e); // stacks the BGP variables (last period of an eqm struct) into an array
uint unstack_bgp_vars(eqm * e, const double * myx); // unstacks BGP vars from array to last period of an eqm struct
uint stack_eqm_vars(double * myx, const eqm * e); // stacks deterministic equilibrium vars into an array
uint unstack_eqm_vars(eqm * e, const double * myx); // unstacks deterministic equillibrium vars
uint stack_stoch_eqm_vars(double * myx, const stoch_eqm * s); // stacks stochastic equilibrium vars into an array
uint unstack_stoch_eqm_vars(stoch_eqm * s, const double * myx); // unstacks stochastic equillibrium vars
uint set_initial_bpg_guess(); // constructs initial guess for a BGP... the "initial guess for the initial guess" function
uint write_bgp_vars(const eqm * e, const params * p); // write BGP vars to text file
uint write_eqm_vars(const eqm * e, char * fname, uint i); // write main deterministic equilibrium vars for country i to file
void set_vars(eqm * e, const params * p, uint t, uint bgp); // sets all the variables for a given period t
uint eval_bgp_conds(const double * myx, double * myf, uint tn); // evaluates the BGP equations
uint solve_bgp(double bb[NC-1]); // solves for the balanced growth path
uint eval_eqm_conds(const double * myx, double * myf, uint tn); // evaluates the deterministic equilibrium conditions
uint solve_eqm(); // solves for the deterministic equilibrium
uint eval_stoch_eqm_conds(const double * myx, double * myf, uint tn);
uint solve_stoch_eqm();
void calc_welfare(eqm * e, const params * p);
void calc_stoch_welfare();
void calc_ce_welfare1();
void calc_ce_welfare2();

// inlined equilibrium equations
static inline double mpk_rk(const params * p, const eqm * e, uint t, uint i, uint s)
{
  if(t>=(NT-1) || p->cap_adj_cost==0)
    {
      return 100.0
	* ( (e->py_t[t][i][s] - DOT_PROD(p->lam[i][s],e->pm_t[t][i],NS))
	    * (p->a_ts[t][i][s]) * (p->alpha[i][s]) * (p->A[i][s]/p->lam_va[i][s])
	    * (pow(e->k_t[t][i][s]/(p->gam_ts[t][i][s] * e->l_t[t][i][s]),p->alpha[i][s]-1.0))
	    - e->rk_t[t][i]/(1.0-p->tauk[i]) );
    }
  else
    {
      return 100.0
	*( (e->py_t[t][i][s] - DOT_PROD(p->lam[i][s],e->pm_t[t][i],NS))
	   * (1.0-p->tauk[i]) *  (p->a_ts[t][i][s]) * (p->alpha[i][s]) * (p->A[i][s]/p->lam_va[i][s])
	   * (pow(e->k_t[t][i][s]/(p->gam_ts[t][i][s] * e->l_t[t][i][s]),p->alpha[i][s]-1.0))
	   - e->pi_t[t-1][i] * (e->cpi_t[t][0]/e->pb_t[t-1])
	   / dphiK(e->is_t[t-1][i][s]/e->k_t[t-1][i][s],p->delta,p->gbgp,p->etaK)
	   - (e->pi_t[t][i]/dphiK(e->is_t[t][i][s]/e->k_t[t][i][s],p->delta,p->gbgp,p->etaK))
	   * ( dphiK(e->is_t[t][i][s]/e->k_t[t][i][s],p->delta,p->gbgp,p->etaK) * e->is_t[t][i][s]/e->k_t[t][i][s]
	       - phiK(e->is_t[t][i][s]/e->k_t[t][i][s],p->delta,p->gbgp,p->etaK) - (1.0-p->delta) ) );
    }
}

static inline double mpl_w(const params * p, const eqm * e, uint t, uint i, uint s)
{
  return (e->py_t[t][i][s] - DOT_PROD(p->lam[i][s],e->pm_t[t][i],NS))
    * (p->a_ts[t][i][s]) * (1.0-p->alpha[i][s]) * (p->A[i][s]/p->lam_va[i][s]) * p->gam_ts[t][i][s]
    * (pow(e->k_t[t][i][s]/(p->gam_ts[t][i][s] * e->l_t[t][i][s]),p->alpha[i][s]))
    - e->w_t[t][i];
}

static inline double mkt_clear_y(const params * p, const eqm * e, uint t, uint i, uint s)
{
  if(!iceberg)
    {
      return e->y_t[t][i][s] - SUM_3D_DIM1(e->m2_t[t],s,i) - SUM_3D_DIM1(e->q2_t[t],s,i);
    }
  else
    {
      return e->y_t[t][i][s] 
	- (1.0+p->ntb_m_ts[t][0][s][i])*e->m2_t[t][0][s][i] - (1.0+p->ntb_f_ts[t][0][s][i])*e->q2_t[t][0][s][i]
	- (1.0+p->ntb_m_ts[t][1][s][i])*e->m2_t[t][1][s][i] - (1.0+p->ntb_f_ts[t][1][s][i])*e->q2_t[t][1][s][i]
	- (1.0+p->ntb_m_ts[t][2][s][i])*e->m2_t[t][2][s][i] - (1.0+p->ntb_f_ts[t][2][s][i])*e->q2_t[t][2][s][i];
    }
}

static inline double mkt_clear_m(const params * p, const eqm * e, uint t, uint i, uint s)
{
  return e->m_t[t][i][s] - e->md_t[t][i][0][s] - e->md_t[t][i][1][s];
}

static inline double mkt_clear_q(const params * p, const eqm * e, uint t, uint i, uint s)
{
  return e->q_t[t][i][s] - e->c_t[t][i][s] - e->i_t[t][i][s];
}

static inline double mucg_mucs(const params * p, const eqm * e, uint t, uint i)
{
  return p->eps[i][0][0]*pow(e->c_t[t][i][0],p->rho-1.0)/e->p_t[t][i][0]
    - p->eps[i][0][1]*pow(e->c_t[t][i][1],p->rho-1.0)/e->p_t[t][i][1];
}

static inline double muc_mul(const params * p, const eqm * e, uint t, uint i)
{
  if(fixl==1 || (fixl==2 && (t<TREF || scenario == 0)))
    {
      return e->ll_t[t][i] - p->lbar_ts[t][i]/3.0;
    }
  else if(fixl==2 && t>=TREF && scenario > 0)
    {
      return e->w_t[t][i] - e->w_t[TREF-1][i];
    }
  else
    {
      return 1.0 * (
		    muc
		    (
		     e->c_t[t][i],
		     e->ll_t[t][i],
		     p->lbar_ts[t][i],
		     p->pope_ts[t][i],
		     p->popw_ts[t][i],
		     p->eps[i][0],
		     p->rho,
		     p->phi[i],
		     p->psi,
		     0)/e->p_t[t][i][0] - 
		    mul
		    (
		     e->c_t[t][i],
		     e->ll_t[t][i],
		     p->lbar_ts[t][i],
		     p->pope_ts[t][i],
		     p->popw_ts[t][i],
		     p->eps[i][0],
		     p->rho,
		     p->phi[i],
		     p->psi)/e->w_t[t][i] );
    }
}

static inline double euler(const params * p, const eqm * e, uint t, uint i)
{
  return 10000.0 * (e->pb_t[t] * muc(
			  e->c_t[t][i],
			  e->ll_t[t][i],
			  p->lbar_ts[t][i],
			  p->pope_ts[t][i],
			  p->popw_ts[t][i],
			  p->eps[i][0],
			  p->rho,
			  p->phi[i],
			  p->psi,
			  0)
    - p->beta[i] * e->cpi_t[t+1][0] * (e->p_t[t][i][0]/e->p_t[t+1][i][0])
    * muc(
	  e->c_t[t+1][i],
	  e->ll_t[t+1][i],
	  p->lbar_ts[t+1][i],
	  p->pope_ts[t+1][i],
	  p->popw_ts[t+1][i],
	  p->eps[i][0],
	  p->rho,
	  p->phi[i],
	  p->psi,
	  0));
}

static inline double bop(const params * p, const eqm * e, uint t, uint i)
{
  double net_iceberg = 0.0;
  if(iceberg==1)
    {
      uint j;
      for(j=0; j<NC; j++)
	{
	  if(j!=i)
	    {
	      uint s;
	      for(s=0; s<NS; s++)
		{
		  net_iceberg = net_iceberg + 
		    + e->py_t[t][i][s] * ( p->ntb_m_ts[t][j][s][i]*e->m2_t[t][j][s][i] + 
					p->ntb_f_ts[t][j][s][i]*e->q2_t[t][j][s][i] )
		    - e->py_t[t][j][s] * ( p->ntb_m_ts[t][i][s][j]*e->m2_t[t][i][s][j] + 
					p->ntb_f_ts[t][i][s][j]*e->q2_t[t][i][s][j] );
		}
	    }
	}
    }

  if(t<NT)
    {
      return SUM(e->nx_t[t][i],NC) + net_iceberg + e->b_t[t][i]*e->cpi_t[t][0] - e->pb_t[t]*e->b_t[t+1][i];
    }
  else
    {
      return SUM(e->nx_t[t][i],NC) + net_iceberg + e->b_t[t][i]*e->cpi_t[t][0] - e->pb_t[t]*e->b_t[t][i]*p->gbgp;
    }
}

static inline double price_norm(const eqm * e, uint t)
{
  return e->cpi_t[t][0] - 1.0;
}

static inline double mkt_clear_i(const params * p, const eqm * e, uint t, uint i)
{
  return e->ii_t[t][i] - e->is_t[t][i][0] - e->is_t[t][i][1];
}

static inline double prod_m_chk(const params * p, const eqm * e, uint t, uint i, uint s)
{
  double tmp = e->m_t[t][i][s] -
    prod_m(e->m2_t[t][i][s], p->M[i][s], p->mu[i][s], p->zeta[i][s]);

  return tmp;
}

static inline double prod_q_chk(const params * p, const eqm * e, uint t, uint i, uint s)
{
  double tmp = e->q_t[t][i][s] -
    prod_q(e->q2_t[t][i][s], p->H[i][s], p->theta[i][s], p->sig[i][s]);

  return tmp;
}

static inline double foc_m2(const params * p, const eqm * e, uint t, uint i, uint s, uint j)
{
  double tmp = e->pm_t[t][i][s] * p->mu[i][s][j] * pow(p->M[i][s],p->zeta[i][s])
    * pow(e->m_t[t][i][s]/e->m2_t[t][i][s][j],1.0 - p->zeta[i][s]);

  tmp = tmp - (1.0+p->tau_m_ts[t][i][s][j]+p->ntb_m_ts[t][i][s][j]) * e->py_t[t][j][s];
 
  if(j != i)
    {
      if(t>0)
	{
	  tmp = tmp - e->w_t[t][i] * (p->etaM/e->m2_t[t-1][i][s][j]) * (e->m2_t[t][i][s][j]/e->m2_t[t-1][i][s][j] - 1.0);
	}
      else
	{
	  tmp = tmp - e->w_t[t][i] * (p->etaM/p->m02[i][s][j]) * (e->m2_t[t][i][s][j]/p->m02[i][s][j] - 1.0);
	}

      tmp = tmp + e->Q_t[t][i] * e->w_t[t+1][i]
	* ( p->etaM*e->m2_t[t+1][i][s][j] / (e->m2_t[t][i][s][j] * e->m2_t[t][i][s][j]) )
	* (e->m2_t[t+1][i][s][j]/e->m2_t[t][i][s][j] - 1.0);
	}

  return 100.0*tmp;
}

static inline double foc_q2(const params * p, const eqm * e, uint t, uint i, uint s, uint j)
{
  double tmp = e->p_t[t][i][s] * p->theta[i][s][j] * pow(p->H[i][s],p->sig[i][s])
    * pow(e->q_t[t][i][s]/e->q2_t[t][i][s][j],1.0 - p->sig[i][s]);

  tmp = tmp - (1.0+p->tau_f_ts[t][i][s][j]+p->ntb_f_ts[t][i][s][j]) * e->py_t[t][j][s];
 
  if(j != i)
    {
      if(t>0)
	{
	  tmp = tmp - e->w_t[t][i] * (p->etaF/e->q2_t[t-1][i][s][j]) * (e->q2_t[t][i][s][j]/e->q2_t[t-1][i][s][j] - 1.0);
	}
      else
	{
	  tmp = tmp - e->w_t[t][i] * (p->etaF/p->q02[i][s][j]) * (e->q2_t[t][i][s][j]/p->q02[i][s][j] - 1.0);
	}

      tmp = tmp + e->Q_t[t][i] * e->w_t[t+1][i]
	* ( p->etaF*e->q2_t[t+1][i][s][j] / (e->q2_t[t][i][s][j] * e->q2_t[t][i][s][j]) )
	* (e->q2_t[t+1][i][s][j]/e->q2_t[t][i][s][j] - 1.0);
	}

  return 100.0*tmp;
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

#endif
