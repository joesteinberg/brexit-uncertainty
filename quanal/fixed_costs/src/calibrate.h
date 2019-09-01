#ifndef __CALIBRATE_H__
#define __CALIBRATE_H__

#include "globals.h"
#include "exporters.h"
#include "solver.h"

double pi_vote;
double pi_brexit;

typedef struct
{
  // ......................................................................
  // parameters
  // ......................................................................

  // household params
  double rbgp; // steady-state interest rate
  double beta; // discount factor
  double psi; // IER
  double phi[NC]; // consumption share in utility
  double Lbar[NC]; // time endowment
  double delta; // depreciation
  double etaK; // capital adjustment cost
  double tauk[NC]; // tax on capital
  uint cap_adj_cost;

  // aggregators
  double zeta; // Armington elasticity
  double izeta; // 1/zeta
  double fzeta; // (zeta-1)/zeta
  double ifzeta; // zeta/(zeta-1)
  double mu[NC][NC]; // Armington shares
  double mue[NC][NC]; // mu^izeta
  double Ybar[NC]; // top-level scale factor
  double Ybare[NC]; // Ybar^((zeta-1)/zeta)
  double theta; // elast. subst. between varieties
  double Ybar2[NC][NC]; // bottom-level scale factor

  // firms
  double gamma[NC]; // value added share in GO
  double alpha; // capital share in VA
  double sigma[NC]; // std dev of idiosyncratic productivity shocks
  //double std[NC]; // (theta-1)*sigma
  double omega[NC]; // fraction of "type-j" exporters
  double kappa0[NC][NC-1]; // export entry cost
  double kappa1[NC][NC-1]; // export continuation cost
  double Nbar[NC]; // mass of firms
  uint Ji[NC][NC-1]; // set of export destinations for country i
  
  // time series parameters
  double tau_ts[NT+1][NC][NC]; // import tariffs
  double ntb_ts[NT+1][NC][NC]; // non-tariff barriers  
  double a_ts[NT+1][NC]; // aggregate productivity

  // ......................................................................
  // base-period equilibrium values
  // ......................................................................
  double iomat[NC+2][NC*2 + 1]; // input-output matrix
  double K0[NC]; // capital
  double L0[NC]; // labor
  double Y0[NC]; // gross output
  double Y20[NC][NC]; // bilateral expenditures
  double C0[NC]; // consumption
  double X0[NC]; // investment
  double M0[NC]; // intermediates
  double VA0[NC]; // value added
  double n0[NC][NC-1]; // type-j exporters
  double EX0[NC][NC]; // bilateral exports
  double IM0[NC][NC]; // bilateral imports
  double NX0[NC][NC]; // bilateral net exports
  double B0[NC]; // bonds
  double R0[NC]; // rental price of capital
   
}params;

params ppp0[NTH];
params ppp1[NTH];
params ppp2[NTH];

#ifdef CH_SIMPLE
params ppp0_ch[NTH];
params ppp1_ch[NTH];
params ppp2_ch[NTH];
#else
params ppp0_ch[NTH];
params ppp1_ch[NHIST_CH-1][NTH];
#endif

typedef struct
{
  params ppp[NHIST];
}stoch_params;

typedef struct
{
  params ppp[NHIST_CH];
}stoch_params_ch;

stoch_params sppp[NTH];
stoch_params_ch sppp_ch[NTH];
exporter_vars ev[NC][2];

uint load_params_from_file();
uint write_params_to_file();
uint copy_params(params * dest, const params * src);
void set_tariffs(stoch_params * sppp);
void set_tariffs_ch(stoch_params_ch * sppp_ch);
void set_tariffs2(params * p, uint scenario);
uint load_iomat();
void set_nontargeted_params();
uint store_base_period_values();
uint calibrate_agg_params();
uint calibrate_hh_params();
uint calibrate_firm_params();
uint eval_exporter_moments(const double * calpars, double * f);
uint calibrate();

int calfunc_f(const gsl_vector * x, void * data, gsl_vector * f);
int calfunc_df(const gsl_vector * x, void * data, gsl_matrix * J);
int calfunc_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);

static inline double arm_agg(const double Y2[NC], double Ybar, const double mue[NC], double fzeta, double ifzeta)
{
  return Ybar * pow(mue[0]*pow(Y2[0],fzeta) + mue[1]*pow(Y2[1],fzeta) + mue[2]*pow(Y2[2],fzeta),ifzeta);
}

static inline double muc(double c, double l, double Lbar, double phi, double psi)
{
  double leisure;
  if(Lbar-l > 0.0001)
    {
      leisure = Lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(Lbar-l));
    }

  return phi * pow( c, psi*phi - 1.0 ) * pow(leisure,(1.0-phi)*psi);
}

static inline double mul(double c, double l, double Lbar, double phi, double psi)
{
  double leisure;
  if(Lbar-l > 0.0001)
    {
      leisure = Lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(Lbar-l));
    }

  return (1.0-phi) * pow( c, psi*phi )  * 
    pow(leisure,(1.0-phi)*psi - 1.0);
}

static inline double phiK(double x, double delta, double etaK)
{
  return (pow(delta,1.0-etaK) * pow(x,etaK) - (1.0-etaK)*(delta))/etaK;
}

static inline double dphiK(double x, double delta, double etaK)
{
  return pow(delta,1.0-etaK) * pow(x,etaK-1.0);
}

#endif
