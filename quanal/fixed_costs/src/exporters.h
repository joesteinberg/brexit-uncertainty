#ifndef __EXPORTERS_H__
#define __EXPORTERS_H__

#include "globals.h"

typedef struct
{
  double zp;
  double zm;

  double Zi;
  double Zp;
  double Zm;
  double Z;

  double Fzp;
  double Fzm;

  double dV;

  double n;
}exporter_vars;

void copy_exporter_vars(exporter_vars * dest, const exporter_vars * src);

void calc_z_Z_Fz(double sigma,
		 double theta,
		 double kappa0,
		 double kappa1,
		 double W,
		 double Dtj,
		 double dVp,
		 exporter_vars * ev);

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
		       exporter_vars * ev);

void init_n(double n0,exporter_vars * ev);

void update_n(double nm,
	      exporter_vars * ev);

void init_dV(double sigma,
		double theta,
		double kappa1,
		double W,
		double Lambda,
		double Dtj,
		exporter_vars * ev);

void update_dV(double kappa0,
	       double kappa1,
	       double W,
	       double Lambda,
	       double Dtj,
	       double dVm,
	       exporter_vars * ev);

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
		     exporter_vars * ev);

uint ev_steady_state(double sigma,
		     double theta,
		     double kappa0,
		     double kappa1,
		     double W,
		     double Lambda,
		     double Dtj,
		     double n0,
		     exporter_vars * ev);

#endif
