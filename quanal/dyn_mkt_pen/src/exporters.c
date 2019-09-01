#ifndef __EXPORTERS_C__
#define __EXPORTERS_C__

#include "exporters.h"

#define max_root_iter 1000
#define root_tol_rel 1.0e-8
#define root_tol_abs 1.0e-8
#define policy_fn_tol 1.0e-8
#define dist_tol 1.0e-8
#define max_policy_fn_iter 500
#define max_dist_iter 1000

///////////////////////////////////////////////////////////////////////////////
// Simple utilities
///////////////////////////////////////////////////////////////////////////////

// copy exporter parameters
void copy_exporter_params(exporter_params * dest, const exporter_params * src)
{
  dest->kappa = src->kappa;
  dest->alpha = src->alpha;
  dest->beta = src->beta;
  dest->psi = src->psi;
  dest->omega = src->omega;
  dest->delta = src->delta;
  dest->xi = src->xi;
  dest->theta = src->theta;
  dest->Zi = src->Zi;

  memcpy((double *)(dest->phi_grid),(const double *)(src->phi_grid),sizeof(double)*NPHI);
  memcpy((double *)(dest->phi_probs),(const double *)(src->phi_probs),sizeof(double)*NPHI);
  memcpy((double *)(dest->phi_cumprobs),(const double *)(src->phi_cumprobs),sizeof(double)*NPHI);
  memcpy((double *)(dest->phi_pow_grid),(const double *)(src->phi_pow_grid),sizeof(double)*NPHI);

  return;
}

// copy vars
void copy_exporter_vars(exporter_vars * dest, const exporter_vars * src)
{
  dest->W = src->W;
  dest->B = src->B;
  dest->Q = src->Q;
  dest->L = src->L;
  dest->La = src->La;
  dest->MC = src->MC;
  dest->Bii=src->Bii;
  dest->Lii=src->Lii;
  dest->tau=src->tau;
  dest->psi_mult=src->psi_mult;

  dest->exports=src->exports;
  dest->export_participation_rate=src->export_participation_rate;
  dest->theoretical_export_participation_rate=src->theoretical_export_participation_rate;
  dest->top5_share = src->top5_share;
  dest->avg_rel_entrant_size = src->avg_rel_entrant_size;
  dest->lf = src->lf;
  dest->exit_rate = src->exit_rate;
  dest->entrant_exit_rate = src->entrant_exit_rate;
  dest->rel_entrant_growth_rate = src->rel_entrant_growth_rate;
  dest->cutoff = src->cutoff;
  dest->Z = src->Z;
  dest->n = src->n;
  dest->startup_cost = src->startup_cost;
  dest->continuation_cost = src->continuation_cost;

  memcpy((double *)(dest->hne),(const double *)(src->hne),sizeof(double)*NPHI);
  memcpy((double *)(dest->diste),(const double *)(src->diste),sizeof(double)*NPHI);
  memcpy((double *)(dest->disti),(const double *)(src->disti),sizeof(double)*NPHI);
  memcpy((double *)(dest->n_grid),(const double *)(src->n_grid),sizeof(double)*NPHI*NN);
  memcpy((double *)(dest->hn),(const double *)(src->hn),sizeof(double)*NPHI*NN);
  memcpy((double *)(dest->dist),(const double *)(src->dist),sizeof(double)*NPHI*NN);
  memcpy((double *)(dest->vx),(const double *)(src->vx),sizeof(double)*NPHI);
  memcpy((double *)(dest->vn),(const double *)(src->vn),sizeof(double)*NPHI);

  return;
}

// set constants (W, Q, B)
void set_exporter_consts1(exporter_vars * ev, double W, double Q, double MC, double psi_mult)
{
  ev->W = W;
  ev->Q = Q;
  ev->MC = MC;
  ev->psi_mult = psi_mult;
}

void set_exporter_consts2(exporter_vars * ev, double L, double P2, double Y2, double Ybar2, double theta, double alpha, double tau)
{
  ev->L = L;
  ev->La = pow(L,alpha);
  ev->B = (1.0/theta) * pow(theta/(theta-1.0),1.0-theta) * 
    pow(ev->MC*tau,1.0-theta) * 
    pow(Ybar2, theta-1.0) * pow(P2,theta) * Y2;
}

void set_exporter_consts3(exporter_vars * ev, double Lii, double P2ii, double Y2ii, double Ybar2ii, double theta)
{
  ev->Lii = Lii;
  ev->Bii = (1.0/theta) * pow(theta/(theta-1.0),1.0-theta) * 
    pow(ev->MC,1.0-theta) * 
    pow(Ybar2ii, theta-1.0) * pow(P2ii,theta) * Y2ii;
}

// discretize productivity grid
void discretize_phi(exporter_params * ep)
{
  // simple evenly spaced method for loglinear
  double inprob = 1.0e-8;
  double lo = gsl_cdf_ugaussian_Pinv(inprob)*ep->kappa*1.5;
  double hi = -gsl_cdf_ugaussian_Pinv(inprob)*ep->kappa*1.5;
  double ome = (hi-lo)/(NPHI-1.0);
  linspace(lo,hi,NPHI,ep->phi_grid);

  double sum=0.0;
  int i;
  for(i=0; i<NPHI; i++)
    {
      ep->phi_probs[i] = (gsl_cdf_ugaussian_P((ep->phi_grid[i]+ome/2.0)/ep->kappa)-
			  gsl_cdf_ugaussian_P((ep->phi_grid[i]-ome/2.0)/ep->kappa));
      sum+=ep->phi_probs[i];
    }

  for(i=0; i<NPHI; i++)
    {
      ep->phi_probs[i] = ep->phi_probs[i]/sum;
    }

  for(i=0; i<NPHI; i++)
    {
      ep->phi_grid[i] = exp(ep->phi_grid[i]);
      ep->phi_pow_grid[i] = pow(ep->phi_grid[i],ep->theta-1.0);
    }

  sum = 0.0;
  for(i=0; i<NPHI; i++)
    {
      sum+=ep->phi_probs[i];
      ep->phi_cumprobs[i] = sum;
    }

  return;
}

// marketing cost
double mp_f(const exporter_params * ep, double np, double n, double L)
{
  double top = 1.0-np;
  double bottom = 1.0 - n*(1.0-ep->delta);
  double frac = top/bottom;

  double retval=(pow(L,ep->alpha)/(1.0-ep->beta)/ep->psi) * (1.0-pow(frac,1.0-ep->beta));
  
  return retval;
}

// derivative of marketing cost with respect to n_{t+1}
double mp_f1(const exporter_params * ep, double np, double n, double La)
{
  double retval=0.0;
  if(full_mkt_pen==1)
    {
      retval=La/ep->psi;
    }
  else
    {
      retval = (La/ep->psi) * pow(1.0-n*(1.0-ep->delta),ep->beta-1.0) * pow(1.0-np,-ep->beta);
    }
  
  return retval;
}

// derivative of marketing cost with respect to n_t
double mp_f2(const exporter_params * ep, double np, double n, double La)
{
  double retval=0.0;
  if(full_mkt_pen==1)
    {
      retval=0.0;
    }
  else
    {
      retval = -(La/ep->psi)*(1.0-ep->delta)*pow(1.0-n*(1.0-ep->delta),ep->beta)*pow(1.0-np,1.0-ep->beta);
    }

  return retval;
}

///////////////////////////////////////////////////////////////////////////////
// Policy function
///////////////////////////////////////////////////////////////////////////////

// analytical entry cutoff
int set_entry_cutoff(const exporter_params * ep, exporter_vars * ev)
{
  double tmp=0.0;
  if(full_mkt_pen==1)
    {
      tmp = ev->W*ev->La*ev->psi_mult/(ep->psi*ev->B)*(1.0-ev->Q*ep->xi*ep->omega);
    }
  else
    {
      tmp = ev->W*ev->La*ev->psi_mult/(ep->psi*ev->B)*(1.0-ev->Q*ep->xi*ep->omega*(1.0-ep->delta));
    }
  ev->cutoff = pow(tmp,1.0/(ep->theta-1.0));
  
  return 0;
}

// convert marketing efficiency from static to dynamic environment
int convert_psi(exporter_params * ep, exporter_vars * ev)
{
  if(full_mkt_pen==1)
    {
      ep->psi = ep->psi*(1.0-ev->Q*ep->xi*ep->omega);
    }
  else
    {
      ep->psi = ep->psi*(1.0-ev->Q*ep->xi*ep->omega*(1.0-ep->delta));
    }
  return 0;
}

// cutoff resulting from solving firm's problem numerically
int check_cutoff(const exporter_params * ep, exporter_vars  * ev)
{
  int ip;
  for(ip=0; ip<NPHI; ip++)
    {
      if(ev->hn[ip][0]>1e-12)
	{
	  break;
	}
    }
  return ip;
}

// initialize policy and value functions
void init_policy(exporter_vars * ev)
{
  int ip;
  for(ip=0; ip<NPHI; ip++)
    {
      expspace(0.0,0.9999,NN,3.0,ev->n_grid[ip]);
      ev->hne[ip] = 0.0;
      set_all_v(ev->hn[ip],NN,0.0);
      ev->vx[ip] = 0.0;
      ev->vn[ip] = 0.0;
    }
}

double interp(gsl_interp_accel * acc, const double *xa, const double *ya, int n, double x)
{
  double x0=0.0;
  double x1=0.0;
  double xd=0.0;
  double q0=0.0;
  double q1=0.0;
  double retval=0.0;

  int ix = gsl_interp_accel_find(acc, xa, n, x);

  if(ix==0)
    {
      x0 = xa[0];
      x1 = xa[1];
      xd = x1-x0;
      q0 = ya[0];
      q1 = ya[1];
    }
  else if(ix==n-1)
    {
      x0 = xa[n-2];
      x1 = xa[n-1];
      xd = x1-x0;
      q0 = ya[n-2];
      q1 = ya[n-1];
    }
  else
   {
      x0 = xa[ix];
      x1 = xa[ix+1];
      xd = x1-x0;
      q0 = ya[ix];
      q1 = ya[ix+1];
    }

  retval = ( q0*(x1-x) + q1*(x-x0) ) / xd;  
  return retval;
}

// Euler equation
double entrant_foc(double np, void * params)
{
  foc_params * params2 = (foc_params *)params;

  const exporter_params * ep = params2->ep;
  exporter_vars * ev = params2->ev;
  exporter_vars * evp = params2->evp;
  gsl_interp_accel * acc = params2->acc;
  double W = ev->W;
  double B = ev->B;
  double Q = ev->Q;
  double La = ev->La;
  int ip = params2->ip;

  double phi_pow = ep->phi_pow_grid[ip];
  double n = 0.0;

  double npp = interp(acc,evp->n_grid[ip],evp->hn[ip],NN,np);
  double lhs = W*ev->psi_mult*mp_f1(ep,np,n,La);
  double rhs = B*phi_pow;

  //rhs = rhs - Q*ep->xi*evp->W*evp->psi_mult*mp_f2(ep,npp,np,La);
  rhs -= Q*ep->xi*ep->omega*evp->W*evp->psi_mult*mp_f2(ep,npp,np,La);
  /*if(evp->hne[ip]>1.0e-8)
    {
      rhs -= Q*ep->xi*(1.0-ep->omega)*evp->W*evp->psi_mult*mp_f2(ep,evp->hne[ip],0.0,La);
      }*/

  return lhs-rhs;
}

double entrant_foc_stoch(double np, void * params)
{
  foc_stoch_params * params2 = (foc_stoch_params *)params;

  const exporter_params * ep = params2->ep;
  exporter_vars * ev = params2->ev;
  exporter_vars * evp1 = params2->evp1;
  exporter_vars * evp2 = params2->evp2;
  double prob = params2->prob;
  gsl_interp_accel * acc1 = params2->acc1;
  gsl_interp_accel * acc2 = params2->acc2;
  double W = ev->W;
  double B = ev->B;
  double Q = ev->Q;
  double La = ev->La;
  int ip = params2->ip;

  double phi_pow = ep->phi_pow_grid[ip];
  double n = 0.0;
  
  double npp1 = interp(acc1,evp1->n_grid[ip],evp1->hn[ip],NN,np);
  double npp2 = interp(acc2,evp2->n_grid[ip],evp2->hn[ip],NN,np);
  double lhs = W*ev->psi_mult*mp_f1(ep,np,n,La);
  double rhs = B*phi_pow;

  rhs -= (prob*Q*ep->xi*evp1->W*evp1->psi_mult*ep->omega*mp_f2(ep,npp1,np,La) -
	  (1.0-prob)*Q*ep->xi*evp2->W*evp2->psi_mult*ep->omega*mp_f2(ep,npp2,np,La));

  return lhs-rhs;
}

double entrant_foc_stoch_idio(double np, void * params)
{
  foc_stoch_idio_params * params2 = (foc_stoch_idio_params *)params;

  const exporter_params * ep = params2->ep;
  exporter_vars * ev = params2->ev;
  exporter_vars ** evplist = params2->evplist;
  double * probs = params2->probs;
  gsl_interp_accel ** acclist = params2->acclist;
  double W = ev->W;
  double B = ev->B;
  double Q = ev->Q;
  double La = ev->La;
  int ip = params2->ip;

  double phi_pow = ep->phi_pow_grid[ip];
  double n = 0.0;
  
  double npplist[6];
  int k;
  for(k=0; k<6; k++)
    {
      npplist[k] = interp(acclist[k],evplist[k]->n_grid[ip],evplist[k]->hn[ip],NN,np);
    }
  double lhs = W*ev->psi_mult*mp_f1(ep,np,n,La);
  double rhs = B*phi_pow;

  for(k=0; k<6; k++)
    {
      rhs -= (probs[k]*Q*ep->xi*evplist[k]->W*evplist[k]->psi_mult*ep->omega*mp_f2(ep,npplist[k],np,La));
    }
  
  return lhs-rhs;
}

// root finder
int find_root_1d(gsl_function * f, double xlo, double xhi, double * x)
{
  int status = 0;
  int iter = 0;
  const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;
  gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
  
  status = gsl_root_fsolver_set(s,f,xlo,xhi);
  if(status)
    {
      printf("Error initializing root-finder!\n");
    }
  else
    {
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate(s);
	  if(status)
	    {
	      printf("Error iterating root-finder!\n");
	      break;
	    }
	  *x = gsl_root_fsolver_root(s);
	  xlo = gsl_root_fsolver_x_lower(s);
	  xhi = gsl_root_fsolver_x_upper(s);
	  status = gsl_root_test_interval(xlo,xhi,root_tol_abs,root_tol_rel);
	}while(status==GSL_CONTINUE && iter<max_root_iter);
    }

  gsl_root_fsolver_free(s);

  return status;
  
}

double invert_foc(const exporter_params * ep,
		  exporter_vars * ev,
		  exporter_vars * evp,
		  int ip,
		  double np,
		  double npp)
{
  double phi_pow = ep->phi_pow_grid[ip];
  
  double rhs = (ev->B*phi_pow - 
		ev->Q*ep->xi*evp->W*evp->psi_mult*ep->omega*mp_f2(ep,npp,np,evp->La));

  double n = (1.0/(1.0-ep->delta))*(1.0 - pow( (pow(1.0-np,ep->beta)*ep->psi/ev->W/ev->psi_mult/ev->La)*rhs,1.0/(ep->beta-1.0) ));

  return n;
}

double invert_foc_stoch(const exporter_params * ep,
			exporter_vars * ev,
			exporter_vars * evp1,
			exporter_vars * evp2,
			int ip,
			double np,
			double npp1,
			double npp2,
			double prob)
{
  double phi_pow = ep->phi_pow_grid[ip];

  double rhs = (ev->B*phi_pow - 
		prob*ev->Q*ep->xi*evp1->W*evp1->psi_mult*ep->omega*mp_f2(ep,npp1,np,evp1->La) - 
		(1.0-prob)*ev->Q*ep->xi*evp2->W*evp2->psi_mult*ep->omega*mp_f2(ep,npp2,np,evp2->La));

  double n = (1.0/(1.0-ep->delta))*(1.0-pow( (pow(1.0-np,ep->beta)*ep->psi/ev->W/ev->psi_mult/ev->La)*rhs,1.0/(ep->beta-1.0) ));

  return n;
}

double invert_foc_stoch_idio(const exporter_params * ep,
			     exporter_vars * ev,
			     exporter_vars * evplist[6],
			     int ip,
			     double np,
			     double npplist[6],
			     double probs[6])
{
  double phi_pow = ep->phi_pow_grid[ip];

  double rhs = (ev->B*phi_pow);
  int k;
  for(k=0; k<6; k++)
    {
      rhs -= probs[k]*ev->Q*ep->xi*evplist[k]->W*evplist[k]->psi_mult*ep->omega*mp_f2(ep,npplist[k],np,evplist[k]->La);
    }

  double n = (1.0/(1.0-ep->delta))*(1.0-pow( (pow(1.0-np,ep->beta)*ep->psi/ev->W/ev->psi_mult/ev->La)*rhs,1.0/(ep->beta-1.0) ));

  return n;
}

int iterate_entrant_policy_fn(const exporter_params * ep,
			      exporter_vars * ev,
			      exporter_vars * evp,
			      double * supnorm)
{

  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  int ip;

  for(ip=0; ip<NPHI; ip++)
    {
      double hne_old = evp->hne[ip];

      //gsl_spline_init(evp->spline,evp->n_grid[ip],evp->hn[ip],NN);
      gsl_interp_accel_reset(acc);

      foc_params params = {ip, ep, ev, evp, acc};
      double lower_bound=1.0e-8;
      double upper_bound = 0.9999-1.0e-8;
	  
      // if FOC evaluated at 0 is positive, marginal cost is too high to enter
      if(entrant_foc(lower_bound,&params)>0.0)
	{
	  ev->hne[ip] = 0.0;
	}
      // otherwise, solve the entrant's Euler equation
      else
	{
	  if(full_mkt_pen==1)
	    {
	      ev->hne[ip]=1.0;
	    }
	  else
	    {
	      if(ip>0)
		{
		  lower_bound = fmax(lower_bound,ev->hne[ip-1]-1.0e-4);
		}

	      gsl_function f;
	      f.function = &entrant_foc;
	      f.params=&params;
	      if(find_root_1d(&f,lower_bound,upper_bound,&(ev->hne[ip])))
		{
		  fprintf(logfile,
			  KRED"\nError solving entrant's first-order condition! ip = %d\n" RESET,ip);
		  gsl_interp_accel_free(acc);
		  return 1;
		}
	    }
	}

      if(fabs(hne_old - ev->hne[ip])> *supnorm)
	{
	  *supnorm = fabs(hne_old - ev->hne[ip]);
	}
    }

  gsl_interp_accel_free(acc);
  return 0;
}

int iterate_entrant_policy_fn_stoch(const exporter_params * ep,
				    exporter_vars * ev,
				    exporter_vars * evp1,
				    exporter_vars * evp2,
				    double prob)
{

  gsl_interp_accel * acc1 = gsl_interp_accel_alloc();
  gsl_interp_accel * acc2 = gsl_interp_accel_alloc();
  int ip;

  for(ip=0; ip<NPHI; ip++)
    {
      gsl_interp_accel_reset(acc1);
      gsl_interp_accel_reset(acc2);
      
      foc_stoch_params params = {ip, ep, ev, evp1, evp2, prob, acc1, acc2};
      double lower_bound=1.0e-8;
      double upper_bound = 0.9999-1.0e-8;
	  
      // if FOC evaluated at 0 is positive, marginal cost is too high to enter
      if(entrant_foc_stoch(lower_bound,&params)>0.0)
	{
	  ev->hne[ip] = 0.0;
	}
      // otherwise, solve the entrant's Euler equation
      else
	{
	  if(full_mkt_pen==1)
	    {
	      ev->hne[ip]=1.0;
	    }
	  else
	    {
	      if(ip>0)
		{
		  lower_bound = fmax(lower_bound,ev->hne[ip-1]-1.0e-4);
		}

	      gsl_function f;
	      f.function = &entrant_foc_stoch;
	      f.params=&params;
	      if(find_root_1d(&f,lower_bound,upper_bound,&(ev->hne[ip])))
		{
		  fprintf(logfile,
			  KRED"\nError solving entrant's first-order condition! ip = %d\n" RESET,ip);
		  gsl_interp_accel_free(acc1);
		  gsl_interp_accel_free(acc2);
		  return 1;
		}
	    }
	}
    }


  gsl_interp_accel_free(acc1);
  gsl_interp_accel_free(acc2);
  return 0;
}

int iterate_entrant_policy_fn_stoch_idio(const exporter_params * ep,
					 exporter_vars * ev,
					 exporter_vars * evplist[6],
					 double probs[6])
{
  gsl_interp_accel * acclist[6];
  int k;
  for(k=0; k<6; k++)
    {
      acclist[k] = gsl_interp_accel_alloc();
    }

  int ip;

  for(ip=0; ip<NPHI; ip++)
    {
      for(k=0; k<6; k++)
	{
	  gsl_interp_accel_reset(acclist[k]);
	}
      
      foc_stoch_idio_params params = {ip, ep, ev, evplist, probs, acclist};
      double lower_bound=1.0e-8;
      double upper_bound = 0.9999-1.0e-8;
	  
      // if FOC evaluated at 0 is positive, marginal cost is too high to enter
      if(entrant_foc_stoch_idio(lower_bound,&params)>0.0)
	{
	  ev->hne[ip] = 0.0;
	}
      // otherwise, solve the entrant's Euler equation
      else
	{
	  if(full_mkt_pen==1)
	    {
	      ev->hne[ip]=1.0;
	    }
	  else
	    {
	      if(ip>0)
		{
		  lower_bound = fmax(lower_bound,ev->hne[ip-1]-1.0e-4);
		}

	      gsl_function f;
	      f.function = &entrant_foc_stoch_idio;
	      f.params=&params;
	      if(find_root_1d(&f,lower_bound,upper_bound,&(ev->hne[ip])))
		{
		  fprintf(logfile,
			  KRED"\nError solving entrant's first-order condition! ip = %d\n" RESET,ip);
		  for(k=0; k<6; k++)
		    {
		      gsl_interp_accel_free(acclist[k]);
		    }
		  return 1;
		}
	    }
	}
    }

  for(k=0; k<6; k++)
    {
      gsl_interp_accel_free(acclist[k]);
    }

  return 0;
}

int iterate_incumbent_policy_fn(const exporter_params * ep,
				exporter_vars * ev,
				exporter_vars * evp,
				double * supnorm)
{
  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  int ip, in, inp;

  for(ip=0; ip<NPHI; ip++)
    {
      // use the entrant's choice as the lower bound on the exogenous grid formation
      // to do this, we need to interpolate to get the current policies defined in the new
      // grid (they were defined on last iteration/period's grid)

      // first, reset the spline on the old grid and old policies (or future ones!)
      gsl_interp_accel_reset(acc);
      
      // now form the new grid based on the curren entrant policy...
      // to avoid extrapolation issues, the lower bound must be no lower than the old one
      expspace(ev->hne[ip],0.9999,NN,2.0,ev->n_grid[ip]);

      // now interpolate the old policies on the new grid
      double hn_tmp[NN];
      for(in=0; in<NN; in++)
	{
	  hn_tmp[in] = interp(acc,evp->n_grid[ip],evp->hn[ip],NN,ev->n_grid[ip][in]);
	}
      
      // now form the endogenous grid implied by interpreting the points on the exogenous
      // grid as optimal policies
      double endo_grid[NN];
      
      // find endogenous grid points by inverting FOC for all np>n0
      for(inp=0; inp<NN; inp++)
	{
	  double n = invert_foc(ep,ev,evp,ip,ev->n_grid[ip][inp],hn_tmp[inp]);
	  
	  // check if FOC does not hold with equality (i.e. firm does not want to to any marketing)
	  // this happens if the new grid point, np, is less than (1-delta)*n
	  if(n*(1.0-ep->delta)>ev->n_grid[ip][inp])
	    {
	      n = ev->n_grid[ip][inp]/(1.0-ep->delta);
	    }
	  
	  // set endogenous grid point
	  endo_grid[inp] = n;
	}    

      // endogenous grid defines a new policy function, hn', which maps endo_grid to exo grid
      // to update hn, which defines policies on the exo grid, we need to interpolate again
      //gsl_spline_init(ev->spline,endo_grid,ev->n_grid[ip],NN);
      gsl_interp_accel_reset(acc);
      for(in=0; in<NN; in++)
	{
	  double hn_new = interp(acc,endo_grid,ev->n_grid[ip],NN,ev->n_grid[ip][in]);
	  if(fabs(hn_new - evp->hn[ip][in])> *supnorm)
	    {
	      *supnorm = fabs(hn_new - evp->hn[ip][in]);
	    }
	  ev->hn[ip][in] = hn_new;
	}
    }

  gsl_interp_accel_free(acc);
  return 0;
}

int iterate_incumbent_policy_fn_stoch(const exporter_params * ep,
				      exporter_vars * ev,
				      exporter_vars * evp1,
				      exporter_vars * evp2,
				      double prob)
{
  gsl_interp_accel * acc1 = gsl_interp_accel_alloc();
  gsl_interp_accel * acc2 = gsl_interp_accel_alloc();
  int ip, in, inp;

  for(ip=0; ip<NPHI; ip++)
    {
      // use the entrant's choice as the lower bound on the exogenous grid formation
      // to do this, we need to interpolate to get the current policies defined in the new
      // grid (they were defined on last iteration/period's grid)

      // first, reset the spline on the old grid and old policies (or future ones!)
      gsl_interp_accel_reset(acc1);
      gsl_interp_accel_reset(acc2);
      
      // now form the new grid based on the curren entrant policy...
      // to avoid extrapolation issues, the lower bound must be no lower than the old one
      expspace(ev->hne[ip],0.9999,NN,2.0,ev->n_grid[ip]);

      // now interpolate the old policies on the new grid
      double hn_tmp1[NN];
      double hn_tmp2[NN];
      for(in=0; in<NN; in++)
	{
	  hn_tmp1[in] = interp(acc1,evp1->n_grid[ip],evp1->hn[ip],NN,ev->n_grid[ip][in]);
	  hn_tmp2[in] = interp(acc2,evp2->n_grid[ip],evp2->hn[ip],NN,ev->n_grid[ip][in]);
	}
      
      // now form the endogenous grid implied by interpreting the points on the exogenous
      // grid as optimal policies
      double endo_grid[NN];
      
      // find endogenous grid points by inverting FOC for all np>n0
      for(inp=0; inp<NN; inp++)
	{
	  double n = invert_foc_stoch(ep,ev,evp1,evp2,ip,ev->n_grid[ip][inp],hn_tmp1[inp],hn_tmp2[inp],prob);
	  
	  // check if FOC does not hold with equality (i.e. firm does not want to to any marketing)
	  // this happens if the new grid point, np, is less than (1-delta)*n
	  if(n*(1.0-ep->delta)>ev->n_grid[ip][inp])
	    {
	      n = ev->n_grid[ip][inp]/(1.0-ep->delta);
	    }
	  
	  // set endogenous grid point
	  endo_grid[inp] = n;
	}    

      // endogenous grid defines a new policy function, hn', which maps endo_grid to exo grid
      // to update hn, which defines policies on the exo grid, we need to interpolate again
      gsl_interp_accel_reset(acc1);
      for(in=0; in<NN; in++)
	{
	  //double hn_new = gsl_spline_eval(ev->spline,ev->n_grid[ip][in],ev->acc);
	  double hn_new = interp(acc1,endo_grid,ev->n_grid[ip],NN,ev->n_grid[ip][in]);
	  ev->hn[ip][in] = hn_new;
	}
    }

  gsl_interp_accel_free(acc1);
  gsl_interp_accel_free(acc2);
  return 0;
}

int iterate_incumbent_policy_fn_stoch_idio(const exporter_params * ep,
					   exporter_vars * ev,
					   exporter_vars * evplist[6],
					   double probs[6])
{
  gsl_interp_accel * acclist[6];
  int k;
  for(k=0; k<6; k++)
    {
      acclist[k] = gsl_interp_accel_alloc();
    }

  int ip, in, inp;

  for(ip=0; ip<NPHI; ip++)
    {
      // use the entrant's choice as the lower bound on the exogenous grid formation
      // to do this, we need to interpolate to get the current policies defined in the new
      // grid (they were defined on last iteration/period's grid)

      // first, reset the spline on the old grid and old policies (or future ones!)
      for(k=0; k<6; k++)
	{
	  gsl_interp_accel_reset(acclist[k]);
	}
      
      // now form the new grid based on the curren entrant policy...
      // to avoid extrapolation issues, the lower bound must be no lower than the old one
      expspace(ev->hne[ip],0.9999,NN,2.0,ev->n_grid[ip]);

      // now interpolate the old policies on the new grid
      double hn_tmp[NN][6];
      for(k=0; k<6; k++)
	{
	  for(in=0; in<NN; in++)
	    {
	      hn_tmp[in][k] = interp(acclist[k],evplist[k]->n_grid[ip],evplist[k]->hn[ip],NN,ev->n_grid[ip][in]);
	    }
	}
      
      // now form the endogenous grid implied by interpreting the points on the exogenous
      // grid as optimal policies
      double endo_grid[NN];
      
      // find endogenous grid points by inverting FOC for all np>n0
      for(inp=0; inp<NN; inp++)
	{
	  double n = invert_foc_stoch_idio(ep,ev,evplist,ip,ev->n_grid[ip][inp],hn_tmp[inp],probs);
	  
	  // check if FOC does not hold with equality (i.e. firm does not want to to any marketing)
	  // this happens if the new grid point, np, is less than (1-delta)*n
	  if(n*(1.0-ep->delta)>ev->n_grid[ip][inp])
	    {
	      n = ev->n_grid[ip][inp]/(1.0-ep->delta);
	    }
	  
	  // set endogenous grid point
	  endo_grid[inp] = n;
	}    

      // endogenous grid defines a new policy function, hn', which maps endo_grid to exo grid
      // to update hn, which defines policies on the exo grid, we need to interpolate again
      gsl_interp_accel_reset(acclist[0]);
      for(in=0; in<NN; in++)
	{
	  //double hn_new = gsl_spline_eval(ev->spline,ev->n_grid[ip][in],ev->acc);
	  double hn_new = interp(acclist[0],endo_grid,ev->n_grid[ip],NN,ev->n_grid[ip][in]);
	  ev->hn[ip][in] = hn_new;
	}
    }

  for(k=0; k<6; k++)
    {
      gsl_interp_accel_free(acclist[k]);
    }

  return 0;
}

// Full market penetration version
int iterate_vf_policy_full_mkt_pen(const exporter_params * ep, exporter_vars * ev, exporter_vars * evp, double * supnorm)
{
  double entry_cost = ev->W*ev->psi_mult*ev->La/ep->psi;
  int ip;
  for(ip=0; ip<NPHI; ip++)
    {
      double vx_old = evp->vx[ip];
      double vn_old = evp->vn[ip];
      ev->vx[ip] = (ev->B*ep->phi_pow_grid[ip] + 
		    ev->Q*ep->xi*ep->omega*evp->vx[ip] +
		    ev->Q*ep->xi*(1.0-ep->omega)*evp->vn[ip]);

      if(ev->vx[ip]-entry_cost > ev->Q*ep->xi*evp->vn[ip])
	{
	  ev->hne[ip]=1.0;
	  ev->vn[ip] = ev->vx[ip]-entry_cost;
	}
      else
	{
	  ev->hne[ip] = 0.0;
	  ev->vn[ip] = ev->Q*ep->xi*evp->vn[ip];
	}
            
      if(fabs(vx_old-ev->vx[ip])>*supnorm)
	{
	  *supnorm=fabs(vx_old-ev->vx[ip]);
	}
      if(fabs(vn_old-ev->vn[ip])>*supnorm)
	{
	  *supnorm=fabs(vn_old-ev->vn[ip]);
	}

    }

  return 0;
}

int iterate_vf_policy_full_mkt_pen_stoch(const exporter_params * ep,
					 exporter_vars * ev,
					 exporter_vars * evp1,
					 exporter_vars * evp2,
					 double prob)
{
  double entry_cost = ev->W*ev->psi_mult*ev->La/ep->psi;
  int ip;
  for(ip=0; ip<NPHI; ip++)
    {
      ev->vx[ip] = (ev->B*ep->phi_pow_grid[ip] + 
		    prob*ev->Q*ep->xi*ep->omega*evp1->vx[ip] + 
		    prob*ev->Q*ep->xi*(1.0-ep->omega)*evp1->vn[ip] + 
		    (1.0-prob)*ev->Q*ep->xi*ep->omega*evp2->vx[ip] + 
		    (1.0-prob)*ev->Q*ep->xi*(1.0-ep->omega)*evp2->vn[ip]);

      double vn_continuation = ev->Q*ep->xi*(prob*evp1->vn[ip] + (1.0-prob)*evp2->vn[ip]);
      if(ev->vx[ip]-entry_cost>vn_continuation)
	{
	  ev->hne[ip]=1.0;
	  ev->vn[ip] = ev->vx[ip]-entry_cost;
	}
      else
	{
	  ev->hne[ip] = 0.0;
	  ev->vn[ip]  =vn_continuation;
	}
    }

  return 0;
}

int iterate_vf_policy_full_mkt_pen_stoch_idio(const exporter_params * ep,
					      exporter_vars * ev,
					      exporter_vars * evplist[6],
					      double probs[6])
{
  double entry_cost = ev->W*ev->psi_mult*ev->La/ep->psi;
  int ip;
  for(ip=0; ip<NPHI; ip++)
    {
      ev->vx[ip] = (ev->B*ep->phi_pow_grid[ip]);
      double vn_continuation = 0.0;
      int k;
      for(k=0; k<6; k++)
	{
	  ev->vx[ip] += probs[k]*ev->Q*ep->xi*ep->omega*evplist[k]->vx[ip] + 
	    probs[k]*ev->Q*ep->xi*(1.0-ep->omega)*evplist[k]->vn[ip];
	  
	  vn_continuation += ev->Q*ep->xi*probs[k]*evplist[k]->vn[ip];
	}
      
      if(ev->vx[ip]-entry_cost>vn_continuation)
	{
	  ev->hne[ip]=1.0;
	  ev->vn[ip] = ev->vx[ip]-entry_cost;
	}
      else
	{
	  ev->hne[ip] = 0.0;
	  ev->vn[ip]  =vn_continuation;
	}
    }

  return 0;
}

// iteration loop
int solve_steady_state_policies(const exporter_params * ep,
				exporter_vars * ev)
{
  time_t start, stop;
  time(&start);

  init_policy(ev);

  int status = 0;
  double supnorm = 999;

  int iter=0;
  do
    {
      iter++;
      supnorm=0.0;
      
      if(full_mkt_pen==0)
	{
	  status = iterate_entrant_policy_fn(ep,ev,ev,&supnorm);
	  if(status)
	    {
	      fprintf(logfile, KRED"\nError iterating entrant policy function!\n RESET");
	      break;
	    }

	  status = iterate_incumbent_policy_fn(ep,ev,ev,&supnorm);
	  if(status)
	    {
	      fprintf(logfile, KRED"\nError iterating incumbent policy function!\n RESET");
	      break;
	    }
	}
      else
	{
	  status = iterate_vf_policy_full_mkt_pen(ep,ev,ev,&supnorm);
	  if(status)
	    {
	      fprintf(logfile, KRED"\nError iterating entrant policy function!\n RESET");
	      break;
	    }
	}
    }
  while(supnorm>policy_fn_tol && iter < max_policy_fn_iter);

  time(&stop);
  
  if(iter==max_policy_fn_iter)
    {
      status=1;
      fprintf(logfile, KRED "\nIteration failed! ||h1-h0|| = %0.4g\n" RESET,supnorm);
    } 

  return status;
}

///////////////////////////////////////////////////////////////////////////////
// Distribution
///////////////////////////////////////////////////////////////////////////////

// initialize distribution
void init_dist(const exporter_params * ep, exporter_vars * ev)
{
  double sum=0.0;
  int ip, in;
  for(ip=0; ip<NPHI; ip++)
    {
      ev->diste[ip]=0.0;
      ev->disti[ip]=0.0;
      sum+=ev->diste[ip];

      for(in=0; in<NN; in++)
	{
	  ev->dist[ip][in]=0.0;
	  sum+=ev->dist[ip][in];
	}
    }

  for(ip=0; ip<NPHI; ip++)
    {
      ev->diste[ip] = ep->phi_probs[ip];
      sum+=ev->diste[ip];
    }
  if(fabs(sum-1.0)>1.0e-8)
    {
      printf("\nNew distribution does not sum to one! sum = %0.4g\n",sum);
    }
}

// distribution iteration driver
int update_dist(const exporter_params * ep,
		exporter_vars * ev,
		exporter_vars * evp,
		double diste[NPHI],
		double dist[NPHI][NN],
		double * supnorm)
{
  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  int in, ip;
  double death_measure=0.0;

  for(ip=0; ip<NPHI; ip++)
    {
      diste[ip]=0.0;

      for(in=0; in<NN; in++)
	{
	  dist[ip][in]=0.0;
	}
    }

  for(ip=0; ip<NPHI; ip++)
    {
      gsl_interp_accel_reset(acc);

      // first update the distributions using the entrant's policy function
      diste[ip] += (1.0-ep->xi)*ev->diste[ip];
      death_measure += ev->diste[ip]*(1.0-ep->xi);
      
      if(ev->hne[ip]<1.0e-8)
	{
	  diste[ip] += ep->xi*ev->diste[ip];
	}
      else
	{
	  // placeholder for extrapolation
	  int ihn = gsl_interp_accel_find(acc,evp->n_grid[ip],NN,ev->hne[ip]);
	  double m1 = (ev->hne[ip]-evp->n_grid[ip][ihn])/(evp->n_grid[ip][ihn+1]-evp->n_grid[ip][ihn]);
	  double m0 = 1.0-m1;
	  dist[ip][ihn] +=  m0*ev->diste[ip]*ep->xi*ep->omega;
	  dist[ip][ihn+1] += m1*ev->diste[ip]*ep->xi*ep->omega;
	  diste[ip] += ep->xi*(1.0-ep->omega)*ev->diste[ip];
	}

      // now do it for incumbents
      for(in=0; in<NN; in++)
	{
	  diste[ip] += (1.0-ep->xi)*ev->dist[ip][in];
	  death_measure += ev->dist[ip][in]*(1.0-ep->xi);

	  int ihn = gsl_interp_accel_find(acc,evp->n_grid[ip],NN,ev->hn[ip][in]);
	  double m1 = (ev->hn[ip][in]-evp->n_grid[ip][ihn])/(evp->n_grid[ip][ihn+1]-evp->n_grid[ip][ihn]);
	  double m0 = 1.0-m1;
	  dist[ip][ihn] += m0*ev->dist[ip][in]*ep->xi*ep->omega;
	  dist[ip][ihn+1]+= m1*ev->dist[ip][in]*ep->xi*ep->omega;
	  diste[ip]+=ep->xi*(1.0-ep->omega)*ev->dist[ip][in];
	}
    }

  /*for(ip=0; ip<NPHI; ip++)
    {
      diste[ip] += death_measure*ep->phi_probs[ip];
      }*/

  gsl_interp_accel_free(acc);

  double sum_old=0.0;
  double sum = 0.0;
  for(ip=0; ip<NPHI; ip++)
    {
      sum_old+=ev->diste[ip];
      sum += diste[ip];
      if(fabs(diste[ip]-ev->diste[ip])>*supnorm)
	{
	  *supnorm = fabs(diste[ip]-ev->diste[ip]);
	}

      for(in=0; in<NN; in++)
	{
	  sum_old += dist[ip][in];
	  sum += dist[ip][in];
	  if(fabs(dist[ip][in]-ev->dist[ip][in])>*supnorm)
	    {
	      *supnorm = fabs(dist[ip][in]-ev->dist[ip][in]);
	    }
	}
    }

  if(fabs(sum-1.0)>1.0e-8)
    {
      printf("\nUpdated distribution does not sum to one! sum = %0.4g, sum_old = %0.4g\n",sum,sum_old);
      return 1;
    }

  return 0;
}

int update_dist_full_mkt_pen(const exporter_params * ep,
			     exporter_vars * ev,
			     exporter_vars * evp,
			     double diste[NPHI],
			     double disti[NPHI],
			     double * supnorm)
{
  int ip;
  double death_measure=0.0;
  for(ip=0; ip<NPHI; ip++)
    {
      diste[ip]=0.0;
      disti[ip]=0.0;
    }
  for(ip=0; ip<NPHI; ip++)
    {
      // a fraction (1-xi) of potential entrants stay there exogenously
      diste[ip] += (1.0-ep->xi)*ev->diste[ip];
      death_measure += ev->diste[ip]*(1.0-ep->xi);
      
      // if the potential entrant chooses to enter, a fraction xi of the
      // potential entrant distribution shifts to the incumbent distribution
      if(ev->hne[ip]<1.0e-8)
	{
	  diste[ip] += ep->xi*ev->diste[ip];
	}
      else
	{
	  disti[ip] += ev->diste[ip]*ep->xi*ep->omega;
	  diste[ip] += ev->diste[ip]*ep->xi*(1.0-ep->omega);
	}

      // a fraction (1-xi) of incumbents exit, and become potential entrants, while
      // the remainder stay incumbents
      diste[ip] += (1.0-ep->xi)*ev->disti[ip];
      death_measure += ev->disti[ip]*(1.0-ep->xi);
      disti[ip] += ep->xi * ep->omega * ev->disti[ip];
      diste[ip] += ep->xi * (1.0-ep->omega) * ev->disti[ip];
    }

  /*for(ip=0; ip<NPHI; ip++)
    {
      diste[ip] += death_measure * ep->phi_probs[ip];
      }*/


  double sum=0.0;
  double sum_old=0.0;
  for(ip=0; ip<NPHI; ip++)
    {
      if(fabs(diste[ip]-ev->diste[ip])>*supnorm)
	{
	  *supnorm=fabs(diste[ip]-ev->diste[ip]);
	}
      if(fabs(disti[ip]-ev->disti[ip])>*supnorm)
	{
	  *supnorm=fabs(disti[ip]-ev->disti[ip]);
	}
      sum+=diste[ip];
      sum+=disti[ip];
      sum_old+=diste[ip];
      sum_old+=disti[ip];
    }

  if(fabs(sum-1.0)>1.0e-8)
    {
      printf("\nUpdated distribution does not sum to one! sum = %0.4g\n",sum);
      return 1;
    }
  else
    {
      return 0;
    }
}

// distribution iteration loop
int solve_steady_state_dist(const exporter_params * ep, exporter_vars * ev)
{
  time_t start, stop;
  int iter=0;
  double tmp_diste[NPHI];
  double tmp_disti[NPHI];
  double tmp_dist[NPHI][NN];
  double supnorm=999;
  int status=0;

  time(&start);
  init_dist(ep,ev);

  do
    {
      iter++;
      supnorm=0.0;

      if(full_mkt_pen==1)
	{
	  status = update_dist_full_mkt_pen(ep,ev,ev,tmp_diste,tmp_disti,&supnorm);
	  memcpy(ev->diste,tmp_diste,NPHI*sizeof(double));
	  memcpy(ev->disti,tmp_disti,NPHI*sizeof(double));
	}
      else
	{
	  status = update_dist(ep,ev,ev,tmp_diste,tmp_dist,&supnorm);
	  memcpy(ev->diste,tmp_diste,NPHI*sizeof(double));
	  memcpy(ev->dist,tmp_dist,NPHI*NN*sizeof(double));
	}

      if(status)
	{
	  printf("\nError iterating distribution!\n");
	  break;
	}      
    }
  while(supnorm>dist_tol && iter < max_dist_iter);

  time(&stop);

  if(iter==max_policy_fn_iter)
    {
      status=1;
      printf("\nIteration failed! ||h1-h0|| = %0.4g\n",supnorm);
    }

  return status;

}

///////////////////////////////////////////////////////////////////////////////
// Computing cross-sectional and lifecycle moments
///////////////////////////////////////////////////////////////////////////////

// qsort comparison (used to find top 5% share)
int compare_obs(const void *a, const void *b)
{
  export_obs *x = (export_obs *)a;
  export_obs *y = (export_obs *)b;

  if (x->v < y->v)
    return -1;
  else if(x->v > y->v)
    return 1;
  else
    return 0;
}

int compute_moments(const exporter_params * ep, exporter_vars * ev)
{
  ev->export_participation_rate=0.0;
  ev->exports=0.0;
  ev->top5_share=0.0;
  ev->avg_rel_entrant_size=0.0;
  ev->lf=0.0;
  ev->Z=0.0;
  ev->n=0.0;
  ev->n2=0.0;
  ev->lf=0.0;
  ev->exit_rate=0.0;
  ev->entrant_exit_rate=0.0;
  ev->rel_entrant_growth_rate = 0.0;
  ev->startup_cost=0.0;
  ev->continuation_cost=0.0;
  set_entry_cutoff(ep,ev);  

  double avg_exporter_size = 0.0;
  double avg_entrant_size = 0.0;
  double entrant_mass=0.0;
  double incumbent_survivor_mass = 0.0;
  double entrant_survivor_mass = 0.0;
  double incumbent_growth_rate = 0.0;
  double entrant_growth_rate = 0.0;
  double f0 = 0.0;
  double f1 = 0.0;
  double f02 = 0.0;
  double f12 = 0.0;
  
  int export_dist_n=0;
  export_obs export_dist[NPHI + NPHI*NN];

  gsl_interp_accel * acc;
  if(calc_all_moments)
    {
      acc = gsl_interp_accel_alloc();
    }
 
  double w2 = 0.0;
  int ip,in;
  for(ip=0; ip<NPHI; ip++)
    {      
      // new entrants
      if(ev->hne[ip]>1.0e-8)
	{
	  ev->lf += ev->psi_mult * mp_f(ep,ev->hne[ip],0.0,ev->L) * ev->diste[ip];
	  ev->export_participation_rate += ev->diste[ip];
	  ev->Z += ev->diste[ip] * ev->hne[ip] * ep->phi_pow_grid[ip];
	  ev->n += ev->diste[ip] * ev->hne[ip];
	  ev->n2 += ev->diste[ip] * ep->phi_pow_grid[ip] * ev->hne[ip];
	  w2 += ev->diste[ip] * ep->phi_pow_grid[ip];
	  
	  if(calc_all_moments)
	    {
	      double py = ev->hne[ip] * ep->theta * ev->B * ep->phi_pow_grid[ip];
	      ev->exports += ev->diste[ip] * (py);
	      ev->exit_rate += ev->diste[ip]*((1.0-ep->xi) + ep->xi*(1.0-ep->omega));
      
	      gsl_interp_accel_reset(acc);
	      //double npp =   gsl_spline_eval(ev->spline,ev->hne[ip],ev->acc);
	      double npp =   interp(acc,ev->n_grid[ip],ev->hn[ip],NN,ev->hne[ip]);
	      double growth = npp/ev->hne[ip]-1.0;

	      export_dist[export_dist_n].v = py;
	      export_dist[export_dist_n].w = ev->diste[ip];
	      export_dist_n++;      

	      entrant_mass += ev->diste[ip];
	      avg_entrant_size += ev->diste[ip] * (py);
	      ev->entrant_exit_rate += ev->diste[ip]*((1.0-ep->xi) + ep->xi*(1.0-ep->omega));
		  
	      entrant_survivor_mass += ev->diste[ip] * ep->xi * ep->omega;
	      entrant_growth_rate += ev->diste[ip] * ep->xi * ep->omega * growth;
	      f0 += (ev->diste[ip] * ev->W *
		     ev->psi_mult * mp_f(ep,ev->hne[ip],0.0,ev->La)
		     );
	      f02 += (ev->diste[ip] * ev->W *
		     ev->psi_mult * mp_f(ep,ev->hne[ip],0.0,ev->La) / py
		     );

	    }
	}

      // incumbents
      for(in=0; in<NN; in++)
	{
	  if(ev->hn[ip][in]>1.0e-8)
	    {
	      ev->lf += ev->psi_mult * mp_f(ep,ev->hn[ip][in],ev->n_grid[ip][in],ev->L) * ev->dist[ip][in];
	      
	      // all exporters
	      ev->export_participation_rate += ev->dist[ip][in];
	      ev->Z += ev->dist[ip][in] * ev->hn[ip][in] * ep->phi_pow_grid[ip];
	      ev->n += ev->dist[ip][in] * ev->hn[ip][in];
	      ev->n2 += ev->dist[ip][in] * ep->phi_pow_grid[ip] * ev->hne[ip];
	      w2 += ev->dist[ip][in] * ep->phi_pow_grid[ip];
	      
	      if(calc_all_moments)
		{
		  double py = ev->hn[ip][in] * ep->theta * ev->B * ep->phi_pow_grid[ip];
		  ev->exports += ev->dist[ip][in] * (py);
		  ev->exit_rate += ev->dist[ip][in]*((1.0-ep->xi) + ep->xi*(1.0-ep->omega));
	      
		  //double npp =   gsl_spline_eval(ev->spline,ev->hn[ip][in],ev->acc);
		  double npp =   interp(acc,ev->n_grid[ip],ev->hn[ip],NN,ev->hn[ip][in]);
		  double growth = npp/ev->hn[ip][in]-1.0;

		  export_dist[export_dist_n].v = py;
		  export_dist[export_dist_n].w = ev->dist[ip][in];
		  export_dist_n++;      

		  incumbent_survivor_mass += ev->dist[ip][in] * ep->xi * ep->omega;
		  incumbent_growth_rate += growth * ev->dist[ip][in] * ep->xi * ep->omega;

		  f1 += (ev->dist[ip][in] * ev->W * ev->psi_mult *
			 mp_f(ep,ev->hn[ip][in],ev->n_grid[ip][in],ev->La)
			 );
		  f12 += (ev->dist[ip][in] * ev->W * ev->psi_mult *
			 mp_f(ep,ev->hn[ip][in],ev->n_grid[ip][in],ev->La) / py
			 );


		}
	    }
	}
    }

  // theoretical_export participation rate (to smooth out calibration function)
  double theoretical_cutoff = ev->cutoff;
  double theoretical_export_participation_rate = 1.0-gsl_cdf_ugaussian_P(log(theoretical_cutoff)/ep->kappa);
  ev->theoretical_export_participation_rate = theoretical_export_participation_rate;
  
  // finalize main cross-sectional moments
  ev->n = ev->n/ev->export_participation_rate;
  ev->n2 = ev->n2/w2;

  if(calc_all_moments)
    {
      avg_exporter_size = avg_exporter_size / ev->export_participation_rate;
      avg_entrant_size = avg_entrant_size / entrant_mass;
      double avg_exports = ev->exports/ev->export_participation_rate;

      // relative entrant size
      ev->avg_rel_entrant_size = avg_entrant_size / avg_exports;

      // exit rates
      ev->exit_rate = ev->exit_rate/ev->export_participation_rate;
      ev->entrant_exit_rate = ev->entrant_exit_rate/entrant_mass;

      // entrant growth
      ev->rel_entrant_growth_rate = ((entrant_growth_rate/entrant_survivor_mass) - 
				     (incumbent_growth_rate/incumbent_survivor_mass));
  
      // top5 share
      int i;
      for(i=0; i<export_dist_n; i++)
	{
	  export_dist[i].w = export_dist[i].w/ev->export_participation_rate;
	}
      qsort(export_dist, export_dist_n, sizeof(export_obs), compare_obs);

      double cw = 0.0;
      double top5_sum = 0.0;
      for(i=0; i<export_dist_n; i++)
	{
	  cw += export_dist[i].w;
	  export_dist[i].cw = cw;
	  if(export_dist[i].cw>=0.95)
	    top5_sum += export_dist[i].w * ev->export_participation_rate * export_dist[i].v;
	}
      ev->top5_share = top5_sum/ev->exports;

      // sunkc costs
      ev->startup_cost = f0/entrant_mass;
      ev->continuation_cost = f1/(ev->export_participation_rate-entrant_mass);
      ev->startup_cost2 = f02/entrant_mass;
      ev->continuation_cost2 = f12/(ev->export_participation_rate-entrant_mass); 

      gsl_interp_accel_free(acc);
    }
  
  return 0;
}

int compute_moments_full_mkt_pen(const exporter_params * ep, exporter_vars * ev)
{
  ev->export_participation_rate=0.0;
  ev->exports=0.0;
  ev->top5_share=0.0;
  ev->avg_rel_entrant_size=0.0;
  ev->lf=0.0;
  ev->Z=0.0;
  ev->n=0.0;
  ev->lf=0.0;
  ev->exit_rate=0.0;
  ev->entrant_exit_rate=0.0;
  ev->rel_entrant_growth_rate = 0.0;
  
  int export_dist_n=0;
  export_obs export_dist[NPHI + NPHI];
 
  int ip;
  for(ip=0; ip<NPHI; ip++)
    {      
      // new entrants
      if(ev->hne[ip]>1.0e-8)
	{
	  ev->lf += (ev->La * ev->psi_mult/ep->psi) * ev->diste[ip];
	  ev->export_participation_rate += ev->diste[ip];
	  ev->Z += ev->diste[ip] * ep->phi_pow_grid[ip];
	  
	  if(calc_all_moments)
	    {
	      double py = ep->theta * ev->B * ep->phi_pow_grid[ip];
	      ev->exports += ev->diste[ip] * (py);
      
	      export_dist[export_dist_n].v = py;
	      export_dist[export_dist_n].w = ev->diste[ip];
	      export_dist_n++;
	    }
	}

      // incumbents
      ev->export_participation_rate += ev->disti[ip];
      ev->Z += ev->disti[ip] * ep->phi_pow_grid[ip];
	      
      if(calc_all_moments)
	{
	  double py = ep->theta * ev->B * ep->phi_pow_grid[ip];
	  ev->exports += ev->disti[ip] * (py);

	  export_dist[export_dist_n].v = py;
	  export_dist[export_dist_n].w = ev->disti[ip];
	  export_dist_n++;
	}
    }

  // theoretical_export participation rate (to smooth out calibration function)
  double theoretical_cutoff = log(ev->cutoff);
  double theoretical_export_participation_rate = 1.0-gsl_cdf_ugaussian_P(theoretical_cutoff/ep->kappa);
  ev->theoretical_export_participation_rate = theoretical_export_participation_rate;

  // finalize main cross-sectional moments
  ev->n = 1.0;

  if(calc_all_moments)
    {
      // relative entrant size
      ev->avg_rel_entrant_size = 1.0;

      // exit rates
      ev->exit_rate = 1.0-ep->xi;
      ev->entrant_exit_rate = 1.0-ep->xi;

      // entrant growth
      ev->rel_entrant_growth_rate = 0.0;
  
      // top5 share
      int i;
      for(i=0; i<export_dist_n; i++)
	{
	  export_dist[i].w = export_dist[i].w/ev->export_participation_rate;
	}
      qsort(export_dist, export_dist_n, sizeof(export_obs), compare_obs);

      double cw = 0.0;
      double top5_sum = 0.0;
      for(i=0; i<export_dist_n; i++)
	{
	  cw += export_dist[i].w;
	  export_dist[i].cw = cw;
	  if(export_dist[i].cw>=0.95)
	    top5_sum += export_dist[i].w * ev->export_participation_rate * export_dist[i].v;
	}
      ev->top5_share = top5_sum/ev->exports;
    }
  
  return 0;
}

int compute_moments_static_version(const exporter_params * ep,
				   exporter_vars * ev)
{
  if(!allexport)
    {  
      // compute entry cutoff using equation (9) in Arkolakis (2010)
      double phi_pow = ev->W * ev->La * ev->psi_mult / (ep->psi * ev->B);
      double entry_cutoff = pow(phi_pow,1.0/(ep->theta-1.0));
      ev->cutoff=entry_cutoff;

      // compute export participation rate
      double log_entry_cutoff = log(entry_cutoff);
      ev->export_participation_rate = 1.0-gsl_cdf_ugaussian_P(log_entry_cutoff/ep->kappa);
      ev->theoretical_export_participation_rate = ev->export_participation_rate;

      if(full_mkt_pen==1)
	{
	  double tmp1 = exp(ep->kappa * ep->kappa * (ep->theta-1.0)*(ep->theta-1.0) / 2.0) * 
	    gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*(ep->theta-1.0) - log_entry_cutoff)/ ep->kappa );

	  // compute total exports using above in equation (14)
	  ev->Z = tmp1;
	  ev->exports = ev->B * ep->theta * ev->Z;
	  ev->n = 1.0;
	  ev->lf = (ev->La*ev->psi_mult/ep->psi) * ev->export_participation_rate;
  
	  // compute top 5 share
	  double log_top5_cutoff = gsl_cdf_ugaussian_Pinv(1.0-0.05*ev->export_participation_rate)*ep->kappa;  
	  double tmp1_5 = exp(ep->kappa * ep->kappa * (ep->theta-1.0)*(ep->theta-1.0) / 2.0) * 
	    gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*(ep->theta-1.0) - log_top5_cutoff)/ ep->kappa );
	  double Z_5 = tmp1_5;
	  double top5_exports = ev->B * ep->theta * Z_5;
	  ev->top5_share = top5_exports/ev->exports;
      
	  if(ev->export_participation_rate<1.0e-6)
	    {
	      ev->top5_share = 1.0;
	    }
	}
      else
	{
	  // compute integral of n*phi^(theta-1) from entry_cutoff to infinity
	  //[
	  // tmp1 = int_{phi*}^{infty} phi^(theta-1) g(phi) dphi
	  double tmp1 = exp(ep->kappa * ep->kappa * (ep->theta-1.0)*(ep->theta-1.0) / 2.0) * 
	    gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*(ep->theta-1.0) - log_entry_cutoff)/ ep->kappa );
	  
	  // tmp2 = int_{phi*}^{infty} (phi*)^{(theta-1)/beta} phi^{(1-theta)/beta + theta-1} g(phi) dphi
	  double k = (1.0-ep->theta)/ep->beta + ep->theta-1.0;
	  double tmp2 = exp(ep->kappa * ep->kappa * k * k / 2.0) * 
	    pow(phi_pow,1.0/ep->beta) * 
	    gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*k - log_entry_cutoff)/ ep->kappa );
	  //]

	  // compute total exports using above in equation (14)
	  ev->Z = tmp1-tmp2;
	  ev->exports = ev->B * ep->theta * ev->Z;

	  // compute average market penetration and labor cost of marketing
	  // int_{phi*}^{infty} [1-(phi*/phi)^{(theta-1)/beta}] g(phi) dphi
	  k = (1.0-ep->theta)/ep->beta;
	  ev->n = exp(ep->kappa * ep->kappa * k * k / 2.0) * 
	    pow(phi_pow,1.0/ep->beta) * 
	    gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*k - log_entry_cutoff)/ ep->kappa );
	  ev->n = 1.0-ev->n;

	  k = (1.0-ep->theta)*(1.0-ep->beta)/ep->beta;
	  ev->lf = (ev->La*ev->psi_mult/(ep->psi*(1.0-ep->beta))) * ev->export_participation_rate;
  
	  ev->lf = ev->lf - ( (ev->La*ev->psi_mult/(ep->psi*(1.0-ep->beta))) * 
			      pow(phi_pow,(1.0-ep->beta)/ep->beta) * 
			      exp(ep->kappa * ep->kappa * k * k / 2.0) *
			      gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*k - log_entry_cutoff)/ ep->kappa ) );

	  // compute top 5 share
	  //[
	  // first compute cutoff for top 5 exporters
	  double log_top5_cutoff = gsl_cdf_ugaussian_Pinv(1.0-0.05*ev->export_participation_rate)*ep->kappa;
  
	  // now compute total exports of firms above the cutoff
	  double tmp1_5 = exp(ep->kappa * ep->kappa * (ep->theta-1.0)*(ep->theta-1.0) / 2.0) * 
	    gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*(ep->theta-1.0) - log_top5_cutoff)/ ep->kappa );

	  k = (1.0-ep->theta)/ep->beta + ep->theta-1.0;
	  double tmp2_5 = exp(ep->kappa * ep->kappa * k * k / 2.0) * 
	    pow(phi_pow,1.0/ep->beta) * 
	    gsl_cdf_ugaussian_P( (ep->kappa*ep->kappa*k - log_top5_cutoff)/ ep->kappa );

	  double Z_5 = tmp1_5-tmp2_5;
	  // compute total exports using above in equation (14)
	  double top5_exports = ev->B * ep->theta * Z_5;
	  ev->top5_share = top5_exports/ev->exports;
	  //]

	  if(ev->export_participation_rate<1.0e-6)
	    {
	      ev->top5_share = 1.0;
	    }
	}

      ev->avg_rel_entrant_size = 1.0;
      ev->exit_rate = 0.0;
      ev->entrant_exit_rate = 0.0;
      ev->rel_entrant_growth_rate = 0.0;
    }
  else
    {
      ev->Z=ep->Zi;
      ev->n=1.0;
      ev->export_participation_rate=1.0;
    }

  return 0;
}

#endif
