#ifndef __EXPORTERS_C__

#define __EXPORTERS_C__

#include "exporters.h"

void copy_exporter_vars(exporter_vars * dest, const exporter_vars * src)
{
  dest->zp = src->zp;
  dest->zm = src->zm;
  dest->Zi = src->Zi;
  dest->Zp = src->Zp;
  dest->Zm = src->Zm;
  dest->Z = src->Z;
  dest->Fzp = src->Fzp;
  dest->Fzm = src->Fzm;
  dest->dV = src->dV;
  dest->n = src->n;
}

void calc_z_Z_Fz(double sigma,
		 double theta,
		 double kappa0,
		 double kappa1,
		 double W,
		 double Dtj,
		 double dVp,
		 exporter_vars * ev)
{
  // now calculate conditional productivities
  if(nokappa || W*kappa0<dVp)
    {
      ev->zp=-HUGE_VAL;
    }
  else
    {
      ev->zp = (1.0/(theta-1.0)) * (log(W*kappa0 - dVp) - log(Dtj));
    }

  if(nokappa || W*kappa1 < dVp)
    {
      ev->zm = -HUGE_VAL;
    }
  else
    {
      ev->zm = (1.0/(theta-1.0)) * (log(W*kappa1 - dVp) - log(Dtj));
    }

  // now calculate conditional productivities
  //double std = sigma*(theta-1.0);
  //double std2 = std*std;
  if(nokappa)
    {
      ev->Zi = 1.0;
      }
  else
    {
      ev->Zi = exp(sigma*sigma*(theta-1.0)*(theta-1.0)/2.0);
    }
  
  if(nokappa || gsl_isinf(ev->zp)==-1)
    {
      ev->Zp = ev->Zi;
    }
  else if(gsl_isinf(ev->zp)==1)
    {
      ev->Zp = 0.0;
    }
  else
    {
      //ev->Zp = ev->Zi * gsl_cdf_ugaussian_P( (std2 - ev->zp/(theta-1.0))/std );
      ev->Zp = ev->Zi * gsl_cdf_ugaussian_P( ((theta-1.0)*sigma*sigma - ev->zp)/sigma );
    }

  if(nokappa || gsl_isinf(ev->zm)==-1)
    {
      ev->Zm = ev->Zi;
    }
  else if(gsl_isinf(ev->zm)==1)
    {
      ev->Zm = 0.0;
    }
  else
    {
      //ev->Zm = ev->Zi * gsl_cdf_ugaussian_P( (std2 - ev->zm/(theta-1.0))/std );
      ev->Zm = ev->Zi * gsl_cdf_ugaussian_P( ((theta-1.0)*sigma*sigma - ev->zm)/sigma );
    }

  // now calculate probabilities
  if(nokappa)
    {
      ev->Fzp = 0.0;
      ev->Fzm = 0.0;
    }
  else
    {
      ev->Fzp = gsl_cdf_ugaussian_P(ev->zp/sigma);
      ev->Fzm = gsl_cdf_ugaussian_P(ev->zm/sigma);
    }

  return;
}

void calc_stoch_z_Z_Fz(double sigma,
		       double theta,
		       double kappa0,
		       double kappa1,
		       double W_good,
		       double W_bad,
		       double Dtj_good,
		       double Dtj_bad,
		       double dVp_good,
		       double dVp_bad,
		       double pi,
		       exporter_vars * ev)
{
  // now calculate conditional productivities
  if(nokappa || pi*(W_good*kappa0-dVp_good) + (1.0-pi)*(W_bad*kappa0-dVp_bad) < 0.0)
    {
      ev->zp=-HUGE_VAL;
    }
  else
    {
      ev->zp = (1.0/(theta-1.0)) * ( log( pi*(W_good*kappa0 - dVp_good) + 
					  (1.0-pi)*(W_bad*kappa0 - dVp_bad) )  
				     - log( pi*Dtj_good + (1.0-pi)*Dtj_bad ) );
    }

  if(nokappa || pi*(W_good*kappa1-dVp_good) + (1.0-pi)*(W_bad*kappa1-dVp_bad) < 0.0)
    {
      ev->zm = -HUGE_VAL;
    }
  else
    {
      ev->zm = (1.0/(theta-1.0)) * ( log( pi*(W_good*kappa1 - dVp_good) + 
					  (1.0-pi)*(W_bad*kappa1 - dVp_bad) )  
				     - log( pi*Dtj_good + (1.0-pi)*Dtj_bad ) );
    }

  // now calculate conditional productivities
  //double std = sigma*(theta-1.0);
  //double std2 = std*std;
  if(nokappa)
    {
      ev->Zi = 1.0;
      }
  else
    {
      ev->Zi = exp(sigma*sigma*(theta-1.0)*(theta-1.0)/2.0);
    }
  
  if(nokappa || gsl_isinf(ev->zp)==-1)
    {
      ev->Zp = ev->Zi;
    }
  else if(gsl_isinf(ev->zp)==1)
    {
      ev->Zp = 0.0;
    }
  else
    {
      //ev->Zp = ev->Zi * gsl_cdf_ugaussian_P( (std2 - ev->zp/(theta-1.0))/std );
      ev->Zp = ev->Zi * gsl_cdf_ugaussian_P( ((theta-1.0)*sigma*sigma - ev->zp)/sigma );
    }

  if(nokappa || gsl_isinf(ev->zm)==-1)
    {
      ev->Zm = ev->Zi;
    }
  else if(gsl_isinf(ev->zm)==1)
    {
      ev->Zm = 0.0;
    }
  else
    {
      //ev->Zm = ev->Zi * gsl_cdf_ugaussian_P( (std2 - ev->zm/(theta-1.0))/std );
      ev->Zm = ev->Zi * gsl_cdf_ugaussian_P( ((theta-1.0)*sigma*sigma - ev->zm)/sigma );
    }

  // now calculate probabilities
  if(nokappa)
    {
      ev->Fzp = 0.0;
      ev->Fzm = 0.0;
    }
  else
    {
      ev->Fzp = gsl_cdf_ugaussian_P(ev->zp/sigma);
      ev->Fzm = gsl_cdf_ugaussian_P(ev->zm/sigma);
    }

  return;
}


void init_n(double n, exporter_vars * ev)
{
  ev->n = n;
}

void update_n(double nm,
	      exporter_vars * ev)
{
  if(nokappa)
    {
      ev->Z = ev->Zi;
      ev->n = 1.0;
    }
  else
    {
      ev->Z = nm * ev->Zm + (1.0-nm) * ev->Zp;
      ev->n = nm * (1.0-ev->Fzm) + (1.0-nm) * (1.0-ev->Fzp);
    }
}

void init_dV(double sigma,
	     double theta,
	     double kappa1,
	     double W,
	     double Lambda,
	     double Dtj,
	     exporter_vars * ev)
{
  if(eqkappa)
    {
      ev->dV = 0.0;
    }
  else
    {
      double Zi = exp(sigma*sigma*(theta-1.0)*(theta-1.0)/2.0);
      //double V0 = Dti*Zi/(1.0-Lambda);
      //double V1 = (Dti*Zi + Dtj*Zi)/(1.0-Lambda);
      ev->dV = Dtj*Zi/(1.0-Lambda);
    }

  return;
}

void update_dV(double kappa0,
	       double kappa1,
	       double W,
	       double Lambda,
	       double Dtj,
	       double dVm,
	       exporter_vars * ev)
{	
  if(nokappa || eqkappa)
    {
      ev->dV=0.0;
    }
  else
    {
      ev->dV = Lambda * (Dtj*(ev->Zm - ev->Zp) 
	+ (ev->Fzp - ev->Fzm) * dVm
			 - (1.0-ev->Fzm)*W*kappa1 + (1.0-ev->Fzp)*W*kappa0);
    }
  
}

void update_stoch_dV(double kappa0,
		     double kappa1,
		     double W_good,
		     double W_bad,
		     double Lambda_good,
		     double Lambda_bad,
		     double Dtj_good,
		     double Dtj_bad,
		     double dVm_good,
		     double dVm_bad,
		     double pi,
		     exporter_vars * ev)
{	
  if(nokappa || eqkappa)
    {
      ev->dV=0.0;
    }
  else
    {
      ev->dV = pi * Lambda_good * ( Dtj_good*(ev->Zm-ev->Zp) + 
				    (ev->Fzp-ev->Fzm)*dVm_good - 
				    (1.0-ev->Fzm)*W_good*kappa1 + 
				    (1.0-ev->Fzp)*W_good*kappa0 ) + 
	(1.0-pi) * Lambda_bad * ( Dtj_bad*(ev->Zm-ev->Zp) + 
				  (ev->Fzp-ev->Fzm)*dVm_bad - 
				  (1.0-ev->Fzm)*W_bad*kappa1 + 
				  (1.0-ev->Fzp)*W_bad*kappa0 );
    }
  
}


uint ev_steady_state(double sigma,
		     double theta,
		     double kappa0,
		     double kappa1,
		     double W,
		     double Lambda,
		     double Dtj,
		     double n0,
		     exporter_vars * ev)
{
  init_n(n0, ev);
  init_dV(sigma,theta,kappa1,W,Lambda,Dtj,ev);

  uint t=0;
  double ddV = +HUGE_VAL;
  double dVm=0.0;
  while(t<10000 && ddV > TINYSQ)
    {
      t++;
      dVm = ev->dV;

      calc_z_Z_Fz(sigma,theta,kappa0,kappa1,W,Dtj,dVm,ev);
      update_n(ev->n,ev);
      update_dV(kappa0,kappa1,W,Lambda,Dtj,ev->dV,ev);
      ddV = fabs(ev->dV-dVm);
    }

  if(t==10000)
    {
      fprintf(logfile,KRED "Steady state value function failed to converge! ddV = %0.3g\n" KNRM,ev->dV - dVm);
      return 1;    
    }

  return 0;
}


#endif
