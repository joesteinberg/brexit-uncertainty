#ifndef __CALIBRATE_H__
#define __CALIBRATE_H__

#include "globals.h"

double pi_vote;
double pi_brexit;


typedef struct
{
  // ......................................................................
  // final demand
  // ......................................................................

  //combine final demand from different sectors... cons is CES, inv is Cobb-Douglas
  // Cons = [eps * goods^rho + (1-eps) * services^rho]^(1/rho)
  // Inv = G * [goods^eps * services^(1-eps)]
  double rho;
  double eps[NC][NF][NS];
  double G[NC];

  // combine final demand from different countries into sector-specific bundles
  // f_s = H * [sum_{j=1}^{NC} (theta_j * f_j)^sig]^(1/sig)
  double sig[NC][NS];
  double theta[NC][NS][NC];
  double H[NC][NS];

  // ......................................................................
  // gross output parameters
  // ......................................................................

  // combine intermediates from different countries into sector-specific bundles
  // M_s = C * [sum_{j=1}^{NC} (mu_j * m_j)^zeta]^(1/zeta)
  double zeta[NC][NS];
  double mu[NC][NS][NC];
  double M[NC][NS];

  // combine value added and intermediate bundles from different sectors... Leontief
  // Gross output = min[VA/lam_va, M_goods/lam_goods, M_svcs/lam_svcs]
  double lam_va[NC][NS];
  double lam[NC][NS][NS];

  // value added... Cobb-Douglas
  // VA = B * [k^alpha * (gam * ell)^(1-alpha)]
  double alpha[NC][NS];
  double A[NC][NS];
  double ggam[NC][NS];
  double gbgp;
  double sc_speed;

  // ......................................................................
  // households
  // ......................................................................
  
  // capital formation
  double delta; // depreciation rate
  double tauk[NC];  // capital tax rate
  double tauk1[NC];
  double rbgp; // balanced growth path real interest rate
  double etaK; // concavity in capital production with Prescott and Lucas (1971) capital formation
  uint cap_adj_cost;
  
  // household preferences
  double beta[NC]; // discount factors
  double psi; // intertemporal elasticity
  double phi[NC]; // consumption share
  double wedge_speed; // convergence rate for discount factor wedge

  // import tariffs: destination-sector-source
  double tau_m_ts[NT+1][NC][NS][NC];
  double tau_f_ts[NT+1][NC][NS][NC];

  // import NTB: destination-sector-source
  double ntb_m_ts[NT+1][NC][NS][NC];
  double ntb_f_ts[NT+1][NC][NS][NC];

  // ......................................................................
  // time series parameters
  // ......................................................................

  // demographic data
  double popt_ts[NT+1][NC];
  double popw_ts[NT+1][NC];
  double pope_ts[NT+1][NC];
  double lbar_ts[NT+1][NC];
  
  // sector-level productivities
  double gam_ts[NT+1][NC][NS];
  double a_ts[NT+1][NC][NS];

  // ......................................................................
  // base-period equilibrium values
  // ......................................................................
  double iomat[NS*NC+2][NS*NC + NF*NC + 1];
  double kk0[NC];
  double b0[NC];
  double r0[NC];
  double ii0[NC];
  double ll0[NC];
  double lbar0[NC];
  double y0[NC][NS];
  double va0[NC][NS];
  double k0[NC][NS];
  double l0[NC][NS];
  double md0[NC][NS][NS];
  double ex0[NC][NC];
  double im0[NC][NC];
  double nx0[NC][NC];
  double c0[NC][NS];
  double i0[NC][NS];
  double m0[NC][NS];
  double m02[NC][NS][NC];
  double q0[NC][NS];
  double q02[NC][NS][NC];
  double im02[NC][NS][NC];
  double ex02[NC][NS][NC];
  double nx02[NC][NS][NC];
  double lshare0[NC][NS];

  // ......................................................................
  // import adjustment cost parameters
  // ......................................................................
  double etaM;
  double etaF;
 
}params;

params ppp0[NTH];
params ppp1[NTH];
params ppp2[NTH];

typedef struct
{
  params ppp[NHIST];
}stoch_params;

stoch_params sppp[NTH];

void set_nontargeted_params(params * p);
void set_tariffs(stoch_params * sppp);
void set_tariffs2(params * p, uint scenario);
uint load_iomat(params * p);
void load_ts_params(params * p);
uint store_base_period_values(params * p);
uint calibrate_prod_params(params * p);
uint calibrate_fin_params(params * p);
uint calibrate_hh_params(params * p);
uint calibrate(params * p);
uint write_params(const params * p);

static inline double prod_go(double va, const double md[NS], double lam_va, const double lam[NS])
{
  return fmin( va/lam_va, fmin(md[0]/lam[0], md[1]/lam[1]) );
}

static inline double prod_va(double k, double l, double gam, double A, double alpha)
{
  return A * pow(k,alpha) * pow(gam*l,(1.0 - alpha));
}

static inline double prod_inv(const double x[NS], const double eps[NS], double G)
{
  return G * pow(x[0],eps[0]) * pow(x[1],eps[1]);
}

static inline double prod_m(const double m2[NC], double M, const double mu[NC], double zeta)
{
  return M * pow( mu[0]*pow(m2[0],zeta) + mu[1]*pow(m2[1],zeta) + mu[2]*pow(m2[2],zeta), 1.0/zeta );
}

static inline double prod_q(const double q2[NC], double H, const double theta[NC], double sig)
{
  return H * pow( theta[0]*pow(q2[0],sig) + theta[1]*pow(q2[1],sig) + theta[2]*pow(q2[2],sig), 1.0/sig );
}

static inline double muc(const double c[NS], double l, double lbar, double ne, double nw, const double eps[NS], double rho, double phi, double psi, uint s)
{
  double leisure;
  if(lbar-l > 0.0001)
    {
      leisure = lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(lbar-l));
    }

  return phi * eps[s] * pow(c[s]/ne,rho-1.0) * 
    pow( eps[0]*pow(c[0]/ne,rho) + eps[1]*pow(c[1]/ne,rho), psi*phi/rho-1.0 ) * 
    pow(leisure/nw,(1.0-phi)*psi) / ne;
}

static inline double mul(const double c[NS], double l, double lbar, double ne, double nw, const double eps[NS], double rho, double phi, double psi)
{
  double leisure;
  if(lbar-l > 0.0001)
    {
      leisure = lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(lbar-l));
    }

  return (1.0-phi) * 
    pow( eps[0]*pow(c[0]/ne,rho) + eps[1]*pow(c[1]/ne,rho), psi*phi/rho )  * 
    pow(leisure/nw,(1.0-phi)*psi - 1.0) / nw;
}

static inline double phiK(double x, double delta, double gbgp, double etaK)
{
  return (pow(delta+gbgp-1.0,1.0-etaK) * pow(x,etaK) - (1.0-etaK)*(delta+gbgp-1.0))/etaK;
}

static inline double dphiK(double x, double delta, double gbgp, double etaK)
{
  return pow(delta+gbgp-1.0,1.0-etaK) * pow(x,etaK-1.0);
}

#endif
