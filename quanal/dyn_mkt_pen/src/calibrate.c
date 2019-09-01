#ifndef __CALIBRATE_C__
#define __CALIBRATE_C__

#include "calibrate.h"

double export_participation_rate_target[NC][2] = {{0.55,0.43},{0.0572,0.37},{0.041,0.1}};
double top5_share_target[NC] = {0.58,0.58,0.58};
double rel_entrant_growth_rate_target[NC] = {0.132280,0.132280,0.132280};
//double new_exporter_size_ratio_target[NC] = {0.43,0.43,0.43};
//double exit_rate_target[NC] = {0.45,0.45,0.45};
//double entrant_exit_rate_target[NC] = {0.64,0.64,0.64};

uint homotopy_times = 30;
double soft_tfp_loss = 0.025;
double hard_tfp_loss = 0.05;
double ch_hi = 5.6;

uint load_params_from_file()
{
  uint i,j;
  params * p = &(ppp0[0]);

  FILE * file = fopen("output/params.txt","rb");
  if(!file)
    {
      return 1;
    }
  
  uint cnt=0;
  uint matched=0;

  cnt += 12;
  matched += fscanf(file,"%lf",&(p->rbgp));
  matched += fscanf(file,"%lf",&(p->beta));
  matched += fscanf(file,"%lf",&(p->psi));
  matched += fscanf(file,"%lf",&(p->delta));
  matched += fscanf(file,"%lf",&(p->etaK));
  matched += fscanf(file,"%lf",&(p->zeta));
  matched +=  fscanf(file,"%lf",&(p->izeta));
  matched += fscanf(file,"%lf",&(p->fzeta));
  matched += fscanf(file,"%lf",&(p->ifzeta));
  matched += fscanf(file,"%lf",&(p->theta));
  matched += fscanf(file,"%lf",&(p->alpha));
  matched += fscanf(file,"%u",&(p->cap_adj_cost));

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->phi[i]));
    }

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->Lbar[i]));
    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->mu[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->mue[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->Ybar[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->Ybare[i]));
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->Ybar2[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->gamma[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 2;
      matched += fscanf(file,"%lf %lf",&(p->ep[i][0].kappa),&(p->ep[i][1].kappa));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->Nbar[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 2;
      matched += fscanf(file,"%lf %lf",&(p->ep[i][0].alpha),&(p->ep[i][1].alpha));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 2;
      matched += fscanf(file,"%lf %lf",&(p->ep[i][0].beta),&(p->ep[i][1].beta));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 2;
      matched += fscanf(file,"%lf %lf",&(p->ep[i][0].psi),&(p->ep[i][1].psi));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 2;
      matched += fscanf(file,"%lf %lf",&(p->ep[i][0].omega),&(p->ep[i][1].omega));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 2;
      matched += fscanf(file,"%lf %lf",&(p->ep[i][0].delta),&(p->ep[i][1].delta));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 2;
      matched += fscanf(file,"%lf %lf",&(p->ep[i][0].xi),&(p->ep[i][1].xi));
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%u",&(p->Ji[i][j]));
	}
    } 

  for(i=0; i<(NC+2); i++)
    {
      for(j=0; j<(NC*2+1); j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->iomat[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->K0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->L0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->Y0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->Y20[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->C0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->X0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->M0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->VA0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->EX0[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->IM0[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->NX0[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->B0[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->R0[i]));
    }

  fclose(file);

  if(matched != cnt)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

uint write_params_to_file()
{
  const params * p = &(ppp0[0]);
  uint i,j;

  FILE * file = fopen("output/params.txt","wb");

  if(!file)
    {
      return 1;
    }

  fprintf(file,"%0.15f\n",p->rbgp);
  fprintf(file,"%0.15f\n",p->beta);
  fprintf(file,"%0.15f\n",p->psi);
  fprintf(file,"%0.15f\n",p->delta);
  fprintf(file,"%0.15f\n",p->etaK);
  fprintf(file,"%0.15f\n",p->zeta);
  fprintf(file,"%0.15f\n",p->izeta);
  fprintf(file,"%0.15f\n",p->fzeta);
  fprintf(file,"%0.15f\n",p->ifzeta);
  fprintf(file,"%0.15f\n",p->theta);
  fprintf(file,"%0.15f\n",p->alpha);
  fprintf(file,"%d\n",p->cap_adj_cost);

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->phi[i]);
    }

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->Lbar[i]);
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file,"%0.15f\n",p->mu[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file,"%0.15f\n",p->mue[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->Ybar[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->Ybare[i]);
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file,"%0.15f\n",p->Ybar2[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->gamma[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f %0.15f\n",p->ep[i][0].kappa,p->ep[i][1].kappa);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->Nbar[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f %0.15f\n",p->ep[i][0].alpha,p->ep[i][1].alpha);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f %0.15f\n",p->ep[i][0].beta,p->ep[i][1].beta);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f %0.15f\n",p->ep[i][0].psi,p->ep[i][1].psi);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f %0.15f\n",p->ep[i][0].omega,p->ep[i][1].omega);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f %0.15f\n",p->ep[i][0].delta,p->ep[i][1].delta);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f %0.15f\n",p->ep[i][0].xi,p->ep[i][1].xi);
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  fprintf(file,"%d\n",p->Ji[i][j]);
	}
    } 

  for(i=0; i<(NC+2); i++)
    {
      for(j=0; j<(NC*2+1); j++)
	{
	  fprintf(file,"%0.15f\n",p->iomat[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->K0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->L0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->Y0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file,"%0.15f\n",p->Y20[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->C0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->X0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->M0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->VA0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file,"%0.15f\n",p->EX0[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file,"%0.15f\n",p->IM0[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file,"%0.15f\n",p->NX0[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->B0[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->R0[i]);
    }

  fclose(file);

  /////////////////////////////////////////////////

  
  file = fopen("output/params_text.txt","wb");

  if(!file)
    {
      return 1;
    }

  fprintf(file,"beta: %0.3f\n",p->beta);
  fprintf(file,"psi: %0.3f\n",p->psi);
  fprintf(file,"delta: %0.3f\n",p->delta);
  fprintf(file,"varphi: %0.3f\n",p->etaK);
  fprintf(file,"zeta: %0.3f\n",p->zeta);
  fprintf(file,"theta: %0.3f\n",p->theta);
  fprintf(file,"alpha: %0.3f\n",p->alpha);
  fprintf(file,"cap_adj_cost: %d\n\n",p->cap_adj_cost);

  fprintf(file,"phi:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->phi[i]);
    }
  fprintf(file,"\n\n");

  fprintf(file,"Lbar:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->Lbar[i]);
    } 
  fprintf(file,"\n\n");

  fprintf(file,"mu:");
  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file," %0.3f",p->mu[i][j]);
	}
      fprintf(file,"\n");
    } 
  fprintf(file,"\n");

  fprintf(file,"Ybar:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->Ybar[i]);
    } 
  fprintf(file,"\n\n");

  fprintf(file,"Ybare:");
  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.3f\n",p->Ybare[i]);
    } 

  fprintf(file,"Ybar2:");
  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  fprintf(file," %0.3f",p->Ybar2[i][j]);
	}
      fprintf(file,"\n");
    } 
  fprintf(file,"\n");

  fprintf(file,"gamma:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->gamma[i]);
    } 
  fprintf(file,"\n\n");

  fprintf(file,"kappa:");
  for(i=0; i<NC; i++)
    {
      fprintf(file,"% 0.4f",p->ep[i][0].kappa);
    } 
  fprintf(file,"\n\n");

  fprintf(file,"lambda:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->ep[i][0].beta);
    } 
  fprintf(file,"\n");

  fprintf(file,"psi:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f %0.3f",p->ep[i][0].psi,p->ep[i][1].psi);
    } 
  fprintf(file,"\n");

  fprintf(file,"omega:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->ep[i][0].omega);
    } 
  fprintf(file,"\n");

  fprintf(file,"chi:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->ep[i][0].delta);
    } 
  fprintf(file,"\n");

  fprintf(file,"xi:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->ep[i][0].xi);
    } 
  fprintf(file,"\n");


  fclose(file);

  return 0;
}


uint copy_params(params * dest, const params * src)
{
  // household params
  dest->rbgp = src->rbgp;
  dest->beta = src->beta;
  dest->psi = src->psi;
  memcpy((double *)(dest->phi),(const double *)(src->phi),sizeof(double)*NC);
  memcpy((double *)(dest->Lbar),(const double *)(src->Lbar),sizeof(double)*NC);
  dest->delta = src->delta;
  dest->etaK = src->etaK;
  memcpy((double *)(dest->tauk),(const double *)(src->tauk),sizeof(double)*NC);
  dest->cap_adj_cost = src->cap_adj_cost;

  // aggregators
  dest->zeta = src->zeta;
  dest->izeta = src->izeta;
  dest->fzeta = src->fzeta;
  dest->ifzeta = src->ifzeta;
  memcpy((double *)(dest->mu),(const double *)(src->mu),sizeof(double)*NC*NC);
  memcpy((double *)(dest->mue),(const double *)(src->mue),sizeof(double)*NC*NC);
  memcpy((double *)(dest->Ybar),(const double *)(src->Ybar),sizeof(double)*NC);
  memcpy((double *)(dest->Ybare),(const double *)(src->Ybare),sizeof(double)*NC);
  dest->theta = src->theta;
  memcpy((double *)(dest->Ybar2),(const double *)(src->Ybar2),sizeof(double)*NC*NC);

  // firms
  memcpy((double *)(dest->gamma),(const double *)(src->gamma),sizeof(double)*NC);
  dest->alpha = src->alpha;
  memcpy((double *)(dest->Nbar),(const double *)(src->Nbar),sizeof(double)*NC);
  memcpy((unsigned int *)(dest->Ji),(const unsigned int *)(src->Ji),sizeof(unsigned int)*NC*(NC-1));

  int ic;
  for(ic=0; ic<NC; ic++)
    {
      copy_exporter_params(&(dest->ep[ic][0]),&(src->ep[ic][0]));
      copy_exporter_params(&(dest->ep[ic][1]),&(src->ep[ic][1]));
    }
  
  // time series parameters
  memcpy((double *)(dest->tau_ts),(const double *)(src->tau_ts),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(dest->ntb_ts),(const double *)(src->ntb_ts),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(dest->psi_mult_ts),(const double *)(src->psi_mult_ts),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(dest->a_ts),(const double *)(src->a_ts),sizeof(double)*(NT+1)*NC);

  // ......................................................................
  // base-period equilibrium values
  // ......................................................................
  memcpy((double *)(dest->iomat),(const double *)(src->iomat),sizeof(double)*(NC+2)*(NC*2+1));
  memcpy((double *)(dest->K0),(const double *)(src->K0),sizeof(double)*NC);
  memcpy((double *)(dest->L0),(const double *)(src->L0),sizeof(double)*NC);
  memcpy((double *)(dest->Y0),(const double *)(src->Y0),sizeof(double)*NC);
  memcpy((double *)(dest->Y20),(const double *)(src->Y20),sizeof(double)*NC*NC);
  memcpy((double *)(dest->C0),(const double *)(src->C0),sizeof(double)*NC);
  memcpy((double *)(dest->X0),(const double *)(src->X0),sizeof(double)*NC);
  memcpy((double *)(dest->M0),(const double *)(src->M0),sizeof(double)*NC);
  memcpy((double *)(dest->VA0),(const double *)(src->VA0),sizeof(double)*NC);
  memcpy((double *)(dest->EX0),(const double *)(src->EX0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->IM0),(const double *)(src->IM0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->NX0),(const double *)(src->NX0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->B0),(const double *)(src->B0),sizeof(double)*NC);
  memcpy((double *)(dest->R0),(const double *)(src->R0),sizeof(double)*NC);
  memcpy((double *)(dest->diste0),(const double *)(src->diste0),sizeof(double)*NC*(NC-1)*NPHI);
  memcpy((double *)(dest->disti0),(const double *)(src->disti0),sizeof(double)*NC*(NC-1)*NPHI);
  memcpy((double *)(dest->dist0),(const double *)(src->dist0),sizeof(double)*NC*(NC-1)*NPHI*NN);

  return 0;
}

void set_nontargeted_params()
{
  params * p = &(ppp0[0]);

  // parameters common across countries

#ifdef _QUARTERLY
  p->rbgp = 0.02/4.0; // quarterly!!
  p->delta = 0.06/4.0;
  p->etaK = 0.8;
#else
  if(high_interest_rate)
    {
      p->rbgp=0.1;
    }
  else
    {
      p->rbgp = 0.02;
    }
  p->delta = 0.06;
  p->etaK = 0.8;
#endif

  p->beta = 1.0 / (1.0 + p->rbgp);

  if(allexport)
    {
      p->zeta = 7.2;
    }
  else if(no_exporter_dyn && full_mkt_pen)
    {
      p->zeta = 6.2;
    }
  else if(no_exporter_dyn)
    {
      p->zeta = 5.0;
    }
  else if(full_mkt_pen)
    {
      p->zeta = 6.2;
    }
  else if(zeta_sens)
    {
      p->zeta = 1.01;
    }
  else
    {
      p->zeta = 5.0;
    }

  p->psi = -1.0;
  if(psi_sens==1)
    {
      p->psi = -4.0;
    }

  p->cap_adj_cost = 1;
  p->alpha = 1.0/3.0;
  p->theta = 5.0;

  SET_ALL_V(p->R0,NC,p->rbgp + p->delta);
  SET_ALL_V(p->tauk,NC,0.0);
  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);

  p->izeta=1.0/p->zeta;
  p->fzeta=(p->zeta-1.0)/p->zeta;
  p->ifzeta=1.0/p->fzeta;

  uint i;
  for(i=0; i<NC; i++)
    {
      if(i==0)
	{
	  p->Ji[i][0]=1;
	  p->Ji[i][1]=2;
	}
      else if(i==1)
	{
	  p->Ji[i][0]=0;
	  p->Ji[i][1]=2;
	}
      else if(i==2)
	{
	  p->Ji[i][0]=0;
	  p->Ji[i][1]=1;
	}
    }

  return;
}

//Average EU-USA reducible NTBs for WIOD industries weighted by EU-UK trade flows
//UK imports from EU:8.70797371111
//EU imports from UK:6.94973604997

//Average EU MFN tariffs for 6-digit HS industries weighted by EU-UK trade flows
//		      UK --> EU: 3.2935061316
//		      EU --> UK: 4.23013777143
//Scaled by goods share of total trade flows
//		      UK --> EU: 2.12371862378
//		      EU --> UK: 3.57662801726

void set_tariffs(stoch_params * sppp_)
{
  uint t, ih;

  // optimistic scenario: no goods tariffs, small non-tariff barriers
  // EU-USA NTBs calculated using same apprroach is in Dhingra et al. (2016)
  // "small" defined as 50% of US-EU NTBs
  ih=1;
  params * p = &((sppp_->ppp)[ih]);
  double ntb_mult = 0.25;
  
  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	}
      if(psi_uncertainty_sens >= 1)
	{
	  p->psi_mult_ts[t][0][1] = 1.27;
	  p->psi_mult_ts[t][1][0] = 1.41;
	}

      if(tfp_flag==1)
	{
	  uint i;
	  for(i=0; i<NC; i++)
	    {
	      p->a_ts[t][i] = 1.0-soft_tfp_loss;
	    }
	}
    }

  // pessimistic scenario: tariffs calculated from EU MFN tariffs as in Dhingra et al. (2016), plus
  // 100% of EU-USA NTBs from above
  ih=2;
  p = &((sppp_->ppp)[ih]);
  ntb_mult = 0.75;

  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	}
      if(psi_uncertainty_sens >= 1)
	{
	  p->psi_mult_ts[t][0][1] = 1.875;
	  p->psi_mult_ts[t][1][0] = 2.5;
	}

      p->tau_ts[t][0][1] = 3.57662801726/100.0;
      p->tau_ts[t][1][0] = 2.12371862378/100.0;

      if(tfp_flag==1)
	{
	  uint i;
	  for(i=0; i<NC; i++)
	    {
	      p->a_ts[t][i] = 1.0-hard_tfp_loss;
	    }
	}
    }
}

void set_tariffs2(params * p, uint scenario)
{
  uint t;

  if(scenario==2 || scenario==4)
    {
      double ntb_mult=0.25;
      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	    {
	      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	    }
	  if(psi_uncertainty_sens >= 1)
	    {
	      p->psi_mult_ts[t][0][1] = 1.27;
	      p->psi_mult_ts[t][1][0] = 1.41;
	    }

	  if(tfp_flag==1)
	    {
	      uint i;
	      for(i=0; i<NC; i++)
		{
		  p->a_ts[t][i] = 1.0-soft_tfp_loss;
		}
	    }
	}
    }
  else if(scenario==3 || scenario==5)
    {
      double ntb_mult=0.75;
      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2) 
	    {
	      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	    }
	  if(psi_uncertainty_sens>=1)
	    {
	      p->psi_mult_ts[t][0][1] = 1.875;
	      p->psi_mult_ts[t][1][0] = 2.5;
	    }

	  p->tau_ts[t][0][1] = 3.57662801726/100.0;
	  p->tau_ts[t][1][0] = 2.12371862378/100.0;

	  if(tfp_flag==1)
	    {
	      uint i;
	      for(i=0; i<NC; i++)
		{
		  p->a_ts[t][i] = 1.0-hard_tfp_loss;
		}
	    }
	}
    }
  else if(scenario==7)
    {
      double ntb_mult=0.25;
      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=TBREXIT; t<TREV; t++)
	{
	  if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	    {
	      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	    }
	  if(psi_uncertainty_sens >= 1)
	    {
	      p->psi_mult_ts[t][0][1] = 1.27;
	      p->psi_mult_ts[t][1][0] = 1.41;
	    }

	  if(tfp_flag==1)
	    {
	      uint i;
	      for(i=0; i<NC; i++)
		{
		  p->a_ts[t][i] = 1.0-soft_tfp_loss;
		}
	    }
	}
    }
  else if(scenario==8)
    {
      double ntb_mult=0.75;
      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=TBREXIT; t<TREV; t++)
	{
	  if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2) 
	    {
	      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	    }
	  if(psi_uncertainty_sens>=1)
	    {
	      p->psi_mult_ts[t][0][1] = 1.875;
	      p->psi_mult_ts[t][1][0] = 2.5;
	    }

	  p->tau_ts[t][0][1] = 3.57662801726/100.0;
	  p->tau_ts[t][1][0] = 2.12371862378/100.0;

	  if(tfp_flag==1)
	    {
	      uint i;
	      for(i=0; i<NC; i++)
		{
		  p->a_ts[t][i] = 1.0-hard_tfp_loss;
		}
	    }
	}
    }
}

void set_tariffs_rb(stoch_params_rb * sppp_)
{
  uint t, ih;

  // permanent soft Brexit
  ih=1;
  params * p = &((sppp_->ppp)[ih]);
  double ntb_mult = 0.25;
  
  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	}
      if(psi_uncertainty_sens >= 1)
	{
	  p->psi_mult_ts[t][0][1] = 1.27;
	  p->psi_mult_ts[t][1][0] = 1.41;
	}

      if(tfp_flag==1)
	{
	  uint i;
	  for(i=0; i<NC; i++)
	    {
	      p->a_ts[t][i] = 1.0-soft_tfp_loss;
	    }
	}
    }

  // permanent hard Brexit
  ih=2;
  p = &((sppp_->ppp)[ih]);
  ntb_mult = 0.75;

  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	}
      if(psi_uncertainty_sens >= 1)
	{
	  p->psi_mult_ts[t][0][1] = 1.875;
	  p->psi_mult_ts[t][1][0] = 2.5;
	}

      p->tau_ts[t][0][1] = 3.57662801726/100.0;
      p->tau_ts[t][1][0] = 2.12371862378/100.0;

      if(tfp_flag==1)
	{
	  uint i;
	  for(i=0; i<NC; i++)
	    {
	      p->a_ts[t][i] = 1.0-hard_tfp_loss;
	    }
	}
    }

  // temporary soft Brexit
  ih=3;
  p = &((sppp_->ppp)[ih]);
  ntb_mult = 0.25;
  
  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<TREV; t++)
    {
      if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	}
      if(psi_uncertainty_sens >= 1)
	{
	  p->psi_mult_ts[t][0][1] = 1.27;
	  p->psi_mult_ts[t][1][0] = 1.41;
	}

      if(tfp_flag==1)
	{
	  uint i;
	  for(i=0; i<NC; i++)
	    {
	      p->a_ts[t][i] = 1.0-soft_tfp_loss;
	    }
	}
    }

  // temporary hard Brexit
  ih=4;
  p = &((sppp_->ppp)[ih]);
  ntb_mult = 0.75;

  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->psi_mult_ts,(NT+1)*NC*NC,1.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<TREV; t++)
    {
      if(psi_uncertainty_sens==0 || psi_uncertainty_sens==2)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;
	}
      if(psi_uncertainty_sens >= 1)
	{
	  p->psi_mult_ts[t][0][1] = 1.875;
	  p->psi_mult_ts[t][1][0] = 2.5;
	}

      p->tau_ts[t][0][1] = 3.57662801726/100.0;
      p->tau_ts[t][1][0] = 2.12371862378/100.0;

      if(tfp_flag==1)
	{
	  uint i;
	  for(i=0; i<NC; i++)
	    {
	      p->a_ts[t][i] = 1.0-hard_tfp_loss;
	    }
	}
    }
}

uint load_iomat()
{
  params * p = &(ppp0[0]);

  uint i, j, got;
  double tmp;
  FILE * file = fopen("input/iomat_balanced2011.txt","rb");
  if(file)
    {
      for(i=0; i<(NC+2); i++)
	{
	  for(j=0; j<(NC+NC+1); j++)
	    {
	      got = fscanf(file,"%lf",&tmp);
	      if(got != 1)
		{
		  fprintf(logfile,KRED "Error reading IO matrix!\n" RESET);
		  fclose(file);
		  return 1;
		}
	      
	      p->iomat[i][j] = tmp;
	    }
	}
      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error loading IO matrix!\n" RESET);
      return 1;
    }
}

uint store_base_period_values()
{
  params * p = &(ppp0[0]);

  double mkt_clear_tol = 1.0e-7;
  uint varow = NC;
  uint gorow = NC+1;

  // initialize all variables we will read from IO matrix to zero
  SET_ALL_V(p->B0,NC,0.0);
  SET_ALL_V(p->K0,NC,0.0);
  SET_ALL_V(p->L0,NC,0.0);
  SET_ALL_V(p->Y0,NC,0.0);
  SET_ALL_V(p->Y20,NC*NC,0.0);
  SET_ALL_V(p->X0,NC,0.0);
  SET_ALL_V(p->C0,NC,0.0);
  SET_ALL_V(p->M0,NC,0.0);
  SET_ALL_V(p->VA0,NC,0.0);
  SET_ALL_V(p->EX0,NC*NC,0.0);
  SET_ALL_V(p->IM0,NC*NC,0.0);
  SET_ALL_V(p->NX0,NC*NC,0.0);

  // first get value added and factors, compute investment and intermediates
  // note: some of these will have to get reshuffled later on... we do this now to
  // check consistency of the IO table
  uint i, j;
  for(i=0; i<NC; i++)
    {
      uint ycol = i;
      uint fcol = NC+i;

      p->Y0[i] = p->iomat[gorow][ycol];
      p->VA0[i] = p->iomat[varow][ycol];
      p->L0[i] = (1.0 - p->alpha) * p->VA0[i];
      p->K0[i] = p->alpha * p->VA0[i] / p->R0[i];
      p->X0[i] = p->delta*p->K0[i];
      p->M0[i] = p->Y0[i]-p->VA0[i];

      // final demand and consumption
      double f = p->iomat[gorow][fcol];
      p->C0[i] = f - p->X0[i];

      // international trade
      for(j=0; j<NC; j++)
	{
	  p->Y20[i][j] = p->iomat[j][ycol] + p->iomat[j][fcol]; 
	  if(j!=i)
	    {
	      p->IM0[i][j] = p->Y20[i][j];
	    }
	}
    }

  // calculate value-added shares in production function and calculate trade variables
  for(i=0; i<NC; i++)
    {
      
      if(leontief)
	{
	  double sq = (p->theta/(p->theta-1.0)) * p->M0[i]/(p->Y20[0][i] + p->Y20[1][i] + p->Y20[2][i]);
	  double br = pow(p->R0[i]/p->alpha,p->alpha) *  
	    pow(1.0/(1.0-p->alpha),1.0-p->alpha);

	  p->gamma[i] = (1.0-sq)/(1.0-sq+br*sq);
	}
      else
	{
	  p->gamma[i] = (p->theta/(p->theta-1.0)) * p->M0[i]/(p->Y20[0][i] + p->Y20[1][i] + p->Y20[2][i]);
	  p->gamma[i] = 1.0-p->gamma[i];
	}
      
      for(j=0; j<NC; j++)
	{
	  p->EX0[i][j] = p->IM0[j][i];
	  p->NX0[i][j] = p->EX0[i][j] - p->IM0[i][j];
	}
    }

  double tmp=0.0;
  for(i=0; i<NC; i++)
    {
      tmp = (p->VA0[i] - (p->C0[i]+p->X0[i]+sum(p->NX0[i],NC)))/p->VA0[i];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "GDP != C+I+NX for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      tmp = (p->Y0[i] - p->Y20[0][i] - p->Y20[1][i] - p->Y20[2][i])/p->Y0[i];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "supply != demand for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      tmp = (p->Y0[i] - (p->VA0[i] + p->M0[i]))/p->Y0[i];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "go != va + m for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      tmp = (p->Y0[i] - p->Y20[i][0] - p->Y20[i][1] - p->Y20[i][2])/p->Y0[i];
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "y != sum(y2) for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}
    }

  return 0;

}

uint calibrate_agg_params()
{
  params * p = &(ppp0[0]);

  uint i,j,cnt,idx;
  double tmp;
  double tmp1[NC];
  
  SET_ALL_V(p->mu,NC*NC,0.0);
  SET_ALL_V(p->Ybar,NC,0.0);

  for(i=0; i<NC; i++)
    {
      idx=2;
      for(j=0; j<NC; j++)
	{
	  tmp1[j] = p->Y20[i][j]/p->Y20[i][idx];
	}
      p->mu[i][idx] = 1.0/sum(tmp1,NC);
      p->mue[i][idx] = pow(p->mu[i][idx],p->izeta);
      cnt=0;
      for(j=0; j<NC; j++)
	{
	  cnt=cnt+1;
	  if(j!=idx)
	    {
	      if(cnt<NC)
		{
		  p->mu[i][j] = p->mu[i][idx]*tmp1[j];
		}
	      else
		{
		  p->mu[i][j] = 1.0-sum(p->mu[i],NC-1);
		}
	      p->mue[i][j] = pow(p->mu[i][j],p->izeta);
	    }
	}

      tmp = pow(DOT_PROD_EX2(p->mu[i],p->Y20[i],NC,p->izeta,p->fzeta),p->ifzeta);
      p->Ybar[i] = p->Y0[i]/tmp;
      p->Ybare[i] = pow(p->Ybar[i],p->fzeta);
    }

  for(i=0; i<NC; i++)
    {
      tmp = p->Y0[i] - arm_agg(p->Y20[i],p->Ybar[i],p->mue[i],p->fzeta,p->ifzeta);
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "Armington aggregator for country %d, error = %f" RESET,i,tmp);
	  return 1;
	}

      for(j=0; j<NC; j++)
	{
	  tmp = p->Ybare[i] * p->mue[i][j] * pow(p->Y0[i]/p->Y20[i][j],p->izeta) - 1.0;
	  if(fabs(tmp)>TINYSQ)
	    {
	      fprintf(logfile,KRED "Armington FOC for country/country %d/%d, , error = %f" RESET,i,j,tmp);
	      return 1;
	    }
	}
    }

  return 0;
}

uint calibrate_hh_params()
{
  params * p = &(ppp0[0]);

  uint i;
  double tmp;
  
  for(i=0; i<NC; i++)
    {
      p->Lbar[i] = 3.0 * p->L0[i];

      if(fixl==1)
	{
	  p->phi[i]=1.0;
	}
      else
	{
	  tmp = (p->Lbar[i] - p->L0[i]) / p->C0[i];
	  p->phi[i] = tmp/(1.0+tmp);
	}

      if(fixl==0)
	{
	  tmp = muc(p->C0[i],p->L0[i],p->Lbar[i],p->phi[i], p->psi) /
	    mul(p->C0[i],p->L0[i],p->Lbar[i],p->phi[i],p->psi) - 1.0;
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "HH intratemp FOC for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }
	}
    }

  return 0;
}

uint calibrate_firm_params(int reset_params)
{
  params * p = &(ppp0[0]);

  uint i,j,ii;
  double MC[NC];
  double D[NC][NC];
  double Db[NC][NC];
  
  // initialize memory for exporter dynamics objects
  /*if(no_exporter_dyn==0 && allexport==0)
    {
      for(i=0; i<NC; i++)
	{
	  for(j=0; j<NC-1; j++)
	    {
	      alloc_splines(&(ev[i][j]));
	    }
	}
	}*/

  // normalize mass of firms to one
  SET_ALL_V(p->Nbar,NC,1.0);

  if(reset_params)
    {
      // initialize exporter params
      for(i=0; i<NC; i++)
	{
	  for(j=0; j<2; j++)
	    {
	      p->ep[i][j].kappa=0.4;
	      p->ep[i][j].alpha=0.0;
	      p->ep[i][j].beta=5.0;
	      p->ep[i][j].delta=0.75;
	      //p->ep[i][j].xi=0.55;
	      if(low_exit_rate_sens==1)
		{
		  p->ep[i][j].omega=0.98;
		  p->ep[i][j].xi=0.96;
		}
	      else
		{
		  p->ep[i][j].omega=0.85;
		  p->ep[i][j].xi=0.6470588;
		}
	      p->ep[i][j].theta=p->theta;
	      p->ep[i][j].Zi = exp(p->ep[i][j].kappa * p->ep[i][j].kappa * (p->theta-1.0)*(p->theta-1.0) / 2.0);

	    }
	}

      p->ep[0][0].psi = 1.0;
      p->ep[0][1].psi = 1.0;
      p->ep[1][0].psi = 0.05;
      p->ep[1][1].psi = 0.25;
      p->ep[2][0].psi = 0.05;
      p->ep[2][1].psi = 0.25;
    }

  if(reset_params)
    {
      // first calculate the domestic demand stuff
      for(i=0; i<NC; i++)
	{
	  if(leontief)
	    {
	      MC[i] = p->gamma[i]*(
				   pow(p->R0[i]/p->alpha,p->alpha) *  
				   pow(1.0/(1.0-p->alpha),1.0-p->alpha)
				   )
		+ (1.0-p->gamma[i])*1.0;
	    }
	  else
	    {
	      MC[i] = pow(p->R0[i]/(p->alpha*p->gamma[i]),p->gamma[i]*p->alpha) *
		pow(1.0/(p->gamma[i]*(1.0-p->alpha)),p->gamma[i]*(1.0-p->alpha)) * 
		pow(1.0/(1.0-p->gamma[i]),1.0-p->gamma[i]);
	    }
	  
	  p->Ybar2[i][i] = pow(p->theta/(p->theta-1.0),p->theta-1.0) * 
	    pow(MC[i],p->theta-1.0) * (1.0/p->Nbar[i]) * (1.0/p->ep[i][0].Zi);
	  p->Ybar2[i][i] = pow(p->Ybar2[i][i],1.0/(p->theta-1.0));

	  D[i][i] = pow(p->Ybar2[i][i], p->theta-1.0) * p->Y20[i][i];
	  Db[i][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * D[i][i] * pow(MC[i],1.0-p->theta);
      
	  double Phome = (1.0/p->Ybar2[i][i]) * pow(p->Nbar[i],1.0/(1.0-p->theta)) * 
	    (p->theta/(p->theta-1.0))*MC[i] * pow(p->ep[i][0].Zi,1.0/(1.0-p->theta));

	  if(fabs(Phome-1.0)>TINYSQ)
	    {
	      fprintf(logfile,"Error! P[%d][%d] = %0.2f\n\n",i,i,Phome);
	      return 1;
	    }
      
	  double Yhome = p->Nbar[i] * Db[i][i] * p->ep[i][0].Zi;
	  if(fabs(Yhome-p->Y20[i][i])>TINYSQ)
	    {
	      fprintf(logfile,"Error! P*Y2[%d][%d] = %0.2f, P*Y20[%d][%d] = %0.2f\n\n",i,i,Yhome,i,i,p->Y20[i][i]);
	      return 1;
	    }

	  Yhome = p->Ybar2[i][i] * pow( MC[i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	    pow(p->Nbar[i],p->theta/(p->theta-1.0)) * D[i][i] * pow(p->ep[i][0].Zi,p->theta/(p->theta-1.0));
	  if(fabs(Yhome-p->Y20[i][i])>TINYSQ)
	    {
	      fprintf(logfile,"Error! Y2[%d][%d] = %0.2f, Y20[%d][%d] = %0.2f\n\n",i,i,Yhome,i,i,p->Y20[i][i]);
	      return 1;
	    }
	}
    }

  if(reset_params)
    {
      // calculate export demand stuff in "approximate" model
      // be equal
      for(i=0; i<NC; i++)
	{
	  for(ii=0; ii<2; ii++)
	    {
	      j=p->Ji[i][ii];
	      double tmp = gsl_cdf_ugaussian_Pinv(1.0-export_participation_rate_target[i][ii])*p->ep[i][ii].kappa;

	      if(allexport)
		{
		  ev[i][ii].Z = p->ep[i][ii].Zi;
		  ev[i][ii].n = 1.0;
		}
	      else if(full_mkt_pen)
		{
		  ev[i][ii].Z = p->ep[i][ii].Zi * 
		    gsl_cdf_ugaussian_P( (p->ep[i][ii].kappa*p->ep[i][ii].kappa*(p->theta-1.0) - tmp)/
					 p->ep[i][ii].kappa );
		}
	      else
		{      
		  // assume all firms that export have 50% market penetration
		  ev[i][ii].Z = 0.5*p->ep[i][ii].Zi * 
		    gsl_cdf_ugaussian_P( (p->ep[i][ii].kappa*p->ep[i][ii].kappa*(p->theta-1.0) - tmp)/
					 p->ep[i][ii].kappa );
		}

	      p->Ybar2[j][i] = pow(p->theta/(p->theta-1.0),p->theta-1.0) * 
		pow(MC[i],p->theta-1.0) * (1.0/(p->Nbar[i])) * (1.0/ev[i][ii].Z);
	      p->Ybar2[j][i] = pow(p->Ybar2[j][i],1.0/(p->theta-1.0));
	  
	      D[j][i] = pow(p->Ybar2[j][i], p->theta-1.0) * p->Y20[j][i];
	      Db[j][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * D[j][i] * pow(MC[i],1.0-p->theta);
      
	      double Yex = p->Nbar[i] * Db[j][i] * ev[i][ii].Z;
	      if(fabs(Yex-p->Y20[j][i])>TINYSQ)
		{
		  fprintf(logfile,"Error! P*Y2[%d][%d] = %0.2f, P*Y20[%d][%d] = %0.2f\n\n",j,i,Yex,j,i,p->Y20[j][i]);
		  return 1;
		}

	      Yex = p->Ybar2[j][i] * pow( MC[i] * (p->theta/(p->theta-1.0)),-p->theta) * 
		pow(p->Nbar[i],p->theta/(p->theta-1.0)) * D[j][i] * pow(ev[i][ii].Z,p->theta/(p->theta-1.0));
	      if(fabs(Yex-p->Y20[j][i])>TINYSQ)
		{
		  fprintf(logfile,"Error! Y2[%d][%d] = %0.2f, Y20[%d][%d] = %0.2f\n\n",j,i,Yex,j,i,p->Y20[j][i]);
		  return 1;
		}
	    }
	}

      for(i=0; i<NC; i++)
	{
	  for(ii=0; ii<2; ii++)
	    {
	      // set constants
	      set_exporter_consts1(&(ev[i][ii]), 1.0, p->beta, MC[i], 1.0);
	      set_exporter_consts2(&(ev[i][ii]), p->L0[j]/p->L0[0], 1.0, p->Y20[j][i], p->Ybar2[j][i], p->theta, p->ep[i][ii].alpha, 1.0);
	      set_exporter_consts3(&(ev[i][ii]), p->L0[i]/p->L0[0], 1.0, p->Y20[i][i], p->Ybar2[i][i], p->theta);
	    }
	}
    } 

  if(no_exporter_dyn || allexport || full_mkt_pen)
    {
      solver_n = 15;
      alloc_solver_mem();
      solver_x->data[0] = log(p->ep[0][0].kappa);
      solver_x->data[1] = log(p->ep[1][0].kappa);
      solver_x->data[2] = log(p->ep[2][0].kappa);
      solver_x->data[3] = log(p->Ybar2[0][1]);
      solver_x->data[4] = log(p->Ybar2[0][2]);
      solver_x->data[5] = log(p->Ybar2[1][0]);
      solver_x->data[6] = log(p->Ybar2[1][2]);
      solver_x->data[7] = log(p->Ybar2[2][0]);
      solver_x->data[8] = log(p->Ybar2[2][1]);
      solver_x->data[9] = log(p->ep[0][0].psi);
      solver_x->data[10] = log(p->ep[0][1].psi);
      solver_x->data[11] = log(p->ep[1][0].psi);
      solver_x->data[12] = log(p->ep[1][1].psi);
      solver_x->data[13] = log(p->ep[2][0].psi);
      solver_x->data[14] = log(p->ep[2][1].psi);
    }
  else
    {
      solver_n = 18;
      alloc_solver_mem();
      solver_x->data[0] = log(p->ep[0][0].kappa);
      solver_x->data[1] = log(p->ep[1][0].kappa);
      solver_x->data[2] = log(p->ep[2][0].kappa);
      solver_x->data[3] = log(p->Ybar2[0][1]);
      solver_x->data[4] = log(p->Ybar2[0][2]);
      solver_x->data[5] = log(p->Ybar2[1][0]);
      solver_x->data[6] = log(p->Ybar2[1][2]);
      solver_x->data[7] = log(p->Ybar2[2][0]);
      solver_x->data[8] = log(p->Ybar2[2][1]);
      solver_x->data[9] = log(p->ep[0][0].psi);
      solver_x->data[10] = log(p->ep[0][1].psi);
      solver_x->data[11] = log(p->ep[1][0].psi);
      solver_x->data[12] = log(p->ep[1][1].psi);
      solver_x->data[13] = log(p->ep[2][0].psi);
      solver_x->data[14] = log(p->ep[2][1].psi);
      solver_x->data[15] = log(-p->ep[0][0].delta/(p->ep[0][0].delta-1.0));
      solver_x->data[16] = log(-p->ep[1][0].delta/(p->ep[1][0].delta-1.0));
      solver_x->data[17] = log(-p->ep[2][0].delta/(p->ep[2][0].delta-1.0));
    }

  uint status = 0;


  if(allexport==1 || eval_calfn_once)
    {
      status = calfunc_f(solver_x,NULL,f0[0]);
    }
  else
    {
      //small_gnewton_step_flag=1;
      solver_verbose=1;
      par=0;
      gsl_multiroot_function_fdf f = {&calfunc_f,&calfunc_df,&calfunc_fdf,solver_n,NULL};

      status = find_root_deriv_mkl(&f);
      if(status)
	{
	  fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
	  status=1;
	}
    
      free_solver_mem();
    }

  if(!allexport)
    {
      fprintf(logfile, KBLU "\n\tKappa:   (%0.4f,%0.4f,%0.4f)\n",p->ep[0][0].kappa,p->ep[1][0].kappa,p->ep[2][0].kappa);
      fprintf(logfile, KBLU "\tBeta:    (%0.4f,%0.4f,%0.4f)\n",p->ep[0][0].beta,p->ep[1][0].beta,p->ep[2][0].beta);
      fprintf(logfile, KBLU "\tPsi:     (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
	      p->ep[0][0].psi,p->ep[0][1].psi,
	      p->ep[1][0].psi,p->ep[1][1].psi,
	      p->ep[2][0].psi,p->ep[2][1].psi);
      if(no_exporter_dyn==0 && full_mkt_pen==0)
	{
	  fprintf(logfile, KBLU "\tDelta:   (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
		  p->ep[0][0].delta,p->ep[0][1].delta,
		  p->ep[1][0].delta,p->ep[1][1].delta,
		  p->ep[2][0].delta,p->ep[2][1].delta);
	}
      fprintf(logfile, KBLU "\tYB2:     (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
	      p->Ybar2[1][0],p->Ybar2[2][0],
	      p->Ybar2[0][1],p->Ybar2[2][1],
	      p->Ybar2[0][2],p->Ybar2[1][2]);
      fprintf(logfile, KBLU "\n\tEPR:     (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
	      ev[0][0].export_participation_rate,ev[0][1].export_participation_rate,
	      ev[1][0].export_participation_rate,ev[1][1].export_participation_rate,
	      ev[2][0].export_participation_rate,ev[2][1].export_participation_rate);
      fprintf(logfile, KBLU "\tTEPR:    (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
	      ev[0][0].theoretical_export_participation_rate,ev[0][1].theoretical_export_participation_rate,
	      ev[1][0].theoretical_export_participation_rate,ev[1][1].theoretical_export_participation_rate,
	      ev[2][0].theoretical_export_participation_rate,ev[2][1].theoretical_export_participation_rate);
      fprintf(logfile, KBLU "\tn:       (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
	      ev[0][0].n,ev[0][1].n,
	      ev[1][0].n,ev[1][1].n,
	      ev[2][0].n,ev[2][1].n);
      fprintf(logfile, KBLU "\tt5:      (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
	      ev[0][0].top5_share,ev[0][1].top5_share,
	      ev[1][0].top5_share,ev[1][1].top5_share,
	      ev[2][0].top5_share,ev[2][1].top5_share);

      if(no_exporter_dyn==0 && full_mkt_pen==0)
	{
	  fprintf(logfile, KBLU "\tEGR:     (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
		  ev[0][0].rel_entrant_growth_rate,ev[0][1].rel_entrant_growth_rate,
		  ev[1][0].rel_entrant_growth_rate,ev[1][1].rel_entrant_growth_rate,
		  ev[2][0].rel_entrant_growth_rate,ev[2][1].rel_entrant_growth_rate);

	  fprintf(logfile, KBLU "\tF0/F1:   (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
		  ev[0][0].startup_cost/ev[0][0].continuation_cost,
		  ev[0][1].startup_cost/ev[0][1].continuation_cost,
		  ev[1][0].startup_cost/ev[1][0].continuation_cost,
		  ev[1][1].startup_cost/ev[1][1].continuation_cost,
		  ev[2][0].startup_cost/ev[2][0].continuation_cost,
		  ev[2][1].startup_cost/ev[2][1].continuation_cost
		  );
	  
	  fprintf(logfile, KBLU "\tF0/F1[2]: (%0.4f,%0.4f), (%0.4f,%0.4f), (%0.4f,%0.4f)\n",
		  ev[0][0].startup_cost/ev[0][0].continuation_cost,
		  ev[0][1].startup_cost/ev[0][1].continuation_cost,
		  ev[1][0].startup_cost/ev[1][0].continuation_cost,
		  ev[1][1].startup_cost/ev[1][1].continuation_cost,
		  ev[2][0].startup_cost/ev[2][0].continuation_cost,
		  ev[2][1].startup_cost/ev[2][1].continuation_cost
		  );

	}
    }

  return status;
}

uint calibrate()
{
  if(recal==1 || allexport || psi_sens || zeta_sens || no_exporter_dyn || full_mkt_pen)
    { 
      BREAK1;
      fprintf(logfile, KBLU "\tAssigning non-targeted parameter values...\n");
      set_nontargeted_params();

      BREAK1;
      fprintf(logfile, KBLU "\tLoading input-output matrix...\n");     
      if(load_iomat())
	{
	  fprintf(logfile, KRED "\nFailed to load iomat!\n" RESET);
	  return 1;
	}

      BREAK1;
      fprintf(logfile, KBLU "\tStoring base-period values...\n");     
      if(store_base_period_values())
	{
	  fprintf(logfile, KRED "\nFailed to store base period values!\n" RESET);
	  return 1;
	}

      BREAK1;
      fprintf(logfile, KBLU "\tCalibrating aggregate production parameters...\n");
      if(calibrate_agg_params())
	{
	  fprintf(logfile, KRED "\nFailed to calibrate agg params!\n" RESET);
	  return 1;
	}

      BREAK1;
      if(no_exporter_dyn==0)
	{
	  no_exporter_dyn=1;
	  fprintf(logfile, KBLU "\tCalibrating exporter-level parameters in simpler no-exporter-dynamics model...\n");
	  if(calibrate_firm_params(1))
	    {
	      fprintf(logfile, KRED "\nFailed to calibrate exporter-level params!\n" RESET);
	      return 1;
	    }
	  no_exporter_dyn=0;
	  int i,ii;
	  for(i=0; i<NC; i++)
	    {
	      for(ii=0; ii<2; ii++)
		{
		  discretize_phi(&(ppp0[0].ep[i][ii]));
		  convert_psi(&(ppp0[0].ep[i][ii]),&ev[i][ii]);
		}
	    }
	  BREAK1;
	}
      fprintf(logfile, KBLU "\tCalibrating exporter-level parameters...\n");
      if(calibrate_firm_params(no_exporter_dyn))
	{
	  fprintf(logfile, KRED "\nFailed to calibrate exporter-level params!\n" RESET);
	  return 1;
	}

      BREAK1;
      fprintf(logfile, KBLU "\tCalibrating household parameters...\n");
      if(calibrate_hh_params())
	{
	  fprintf(logfile, KRED "\nFailed to calibrate hh params!\n" RESET);
	  return 1;
	}

      if(allexport==0 && no_exporter_dyn==0 && psi_sens==0 && zeta_sens==0)
	{
	  BREAK1;
	  fprintf(logfile, KBLU "\tWriting parameters to disk...\n");
	  if(write_params_to_file())
	    {
	      fprintf(logfile, KRED "\nFailed to write parameters to file!\n" RESET);
	      return 1;
	    }
	}
    }
  else
    {
      BREAK1;
      fprintf(logfile, KBLU "\tReading pre-calibrated parameters from disk...\n");
      if(load_params_from_file())
	{
	  fprintf(logfile, KRED "\nFailed to read parameters from file!\n" RESET);
	  return 1;
	}
    }
  
  uint it;
  for(it=0; it<NTH; it++)
    {
      if(it>0 && copy_params(&(ppp0[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp0!\n" RESET);
	  return 1;
	}
      if(copy_params(&(ppp1[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp1!\n" RESET);
	  return 1;
	}
      if(copy_params(&(ppp2[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp2!\n" RESET);
	  return 1;
	}
      if(copy_params(&(ppp3[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp1!\n" RESET);
	  return 1;
	}
      if(copy_params(&(ppp4[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp2!\n" RESET);
	  return 1;
	}

      stoch_params * sp = &(sppp[it]);
      uint ih;
      for(ih=0; ih<NHIST; ih++)
	{
	  if(copy_params(&((sp->ppp)[ih]),&(ppp0[0])))
	    {
	      fprintf(logfile, KRED "\nFailed to copy stoch_params!\n" RESET);
	      return 1;
	    }
	}

      stoch_params_rb * sp_rb = &(sppp_rb[it]);
      for(ih=0; ih<NHIST_RB; ih++)
	{
	  if(copy_params(&((sp_rb->ppp)[ih]),&(ppp0[0])))
	    {
	      fprintf(logfile, KRED "\nFailed to copy stoch_params!\n" RESET);
	      return 1;
	    }
	}
    }

  return 0;
  }

///////////////////////////////

uint eval_exporter_moments(const double * calpars, double * myf)
{
  params * p = &(ppp0[0]);

  if(allexport==1 || no_exporter_dyn==1 || full_mkt_pen==1)
    {
      p->ep[0][0].kappa = exp(calpars[0]);
      p->ep[1][0].kappa = exp(calpars[1]);
      p->ep[2][0].kappa = exp(calpars[2]);
      p->Ybar2[0][1] = exp(calpars[3]);
      p->Ybar2[0][2] = exp(calpars[4]);
      p->Ybar2[1][0] = exp(calpars[5]);
      p->Ybar2[1][2] = exp(calpars[6]);
      p->Ybar2[2][0] = exp(calpars[7]);
      p->Ybar2[2][1] = exp(calpars[8]);
      p->ep[0][0].psi = exp(calpars[9]);
      p->ep[0][1].psi = exp(calpars[10]);
      p->ep[1][0].psi = exp(calpars[11]);
      p->ep[1][1].psi = exp(calpars[12]);
      p->ep[2][0].psi = exp(calpars[13]);
      p->ep[2][1].psi = exp(calpars[14]);
    }
  else
    {
      p->ep[0][0].kappa = exp(calpars[0]);
      p->ep[1][0].kappa = exp(calpars[1]);
      p->ep[2][0].kappa = exp(calpars[2]);
      p->Ybar2[0][1] = exp(calpars[3]);
      p->Ybar2[0][2] = exp(calpars[4]);
      p->Ybar2[1][0] = exp(calpars[5]);
      p->Ybar2[1][2] = exp(calpars[6]);
      p->Ybar2[2][0] = exp(calpars[7]);
      p->Ybar2[2][1] = exp(calpars[8]);
      p->ep[0][0].psi = exp(calpars[9]);
      p->ep[0][1].psi = exp(calpars[10]);
      p->ep[1][0].psi = exp(calpars[11]);
      p->ep[1][1].psi = exp(calpars[12]);
      p->ep[2][0].psi = exp(calpars[13]);
      p->ep[2][1].psi = exp(calpars[14]);
      p->ep[0][0].delta = 1.0/(1.0+exp(-calpars[15]));
      p->ep[1][0].delta = 1.0/(1.0+exp(-calpars[16]));
      p->ep[2][0].delta = 1.0/(1.0+exp(-calpars[17]));
    }
  
  uint i, ii, j;
  
  for(i=0; i<NC; i++)
    {
      p->ep[i][1].kappa = p->ep[i][0].kappa;
      p->ep[i][0].Zi = exp(p->ep[i][0].kappa * p->ep[i][0].kappa * (p->theta-1.0)*(p->theta-1.0) / 2.0);
      p->ep[i][1].Zi = exp(p->ep[i][1].kappa * p->ep[i][1].kappa * (p->theta-1.0)*(p->theta-1.0) / 2.0);
      if(allexport==0 && no_exporter_dyn==0 && full_mkt_pen==0)
	{
	  p->ep[i][1].delta = p->ep[i][0].delta;
	}
    }

  double MC[NC];
  double D[NC][NC];
  double Db[NC][NC];
  double Dh[NC][NC];

  // recalculate export demand scale factors
  for(i=0; i<NC; i++)
    {
      if(leontief)
	{
	  MC[i] = p->gamma[i]*(
			    pow(p->R0[i]/p->alpha,p->alpha) *
			    pow(1.0/(1.0-p->alpha),1.0-p->alpha)
			    )
	    + (1.0-p->gamma[i])*1.0;
 	}
      else
	{
	  MC[i] = pow(p->R0[i]/(p->alpha*p->gamma[i]),p->gamma[i]*p->alpha) * 
	    pow(1.0/(p->gamma[i]*(1.0-p->alpha)),p->gamma[i]*(1.0-p->alpha)) * 
	    pow(1.0/(1.0-p->gamma[i]),1.0-p->gamma[i]);
	}

      for(j=0; j<NC; j++)
	{
	  D[j][i] = pow(p->Ybar2[j][i], p->theta-1.0) * p->Y20[j][i];
	  Db[j][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * D[j][i] * pow(MC[i],1.0-p->theta);
	  Dh[j][i] = pow(p->theta/(p->theta-1.0),-p->theta) * D[j][i] * pow(MC[i],1.0-p->theta);
      	}
    }

  for(i=0; i<NC; i++)
    {
      if(allexport)
	{
	  ev[i][0].Z=p->ep[i][0].Zi;
	  ev[i][1].Z=p->ep[i][1].Zi;
	  ev[i][0].n=1.0;
	  ev[i][1].n=1.0;
	}

      p->Ybar2[i][i] = pow(p->theta/(p->theta-1.0),p->theta-1.0) * 
	pow(MC[i],p->theta-1.0) * (1.0/p->Nbar[i]) * (1.0/p->ep[i][0].Zi);
      p->Ybar2[i][i] = pow(p->Ybar2[i][i],1.0/(p->theta-1.0));
      
      D[i][i] = pow(p->Ybar2[i][i], p->theta-1.0) * p->Y20[i][i];
      Db[i][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * D[i][i] * pow(MC[i],1.0-p->theta);
      
      double Phome = (1.0/p->Ybar2[i][i]) * pow(p->Nbar[i],1.0/(1.0-p->theta)) * 
	(p->theta/(p->theta-1.0))*MC[i] * pow(p->ep[i][0].Zi,1.0/(1.0-p->theta));

      if(fabs(Phome-1.0)>TINYSQ)
	{
	  fprintf(logfile,"Error! P[%d][%d] = %0.2f\n\n",i,i,Phome);
	  return 1;
	}
      
      double Yhome = p->Nbar[i] * Db[i][i] * p->ep[i][0].Zi;
      if(fabs(Yhome-p->Y20[i][i])>TINYSQ)
	{
	  fprintf(logfile,"Error! P*Y2[%d][%d] = %0.2f, P*Y20[%d][%d] = %0.2f\n\n",i,i,Yhome,i,i,p->Y20[i][i]);
	  return 1;
	}

      Yhome = p->Ybar2[i][i] * pow( MC[i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	pow(p->Nbar[i],p->theta/(p->theta-1.0)) * D[i][i] * pow(p->ep[i][0].Zi,p->theta/(p->theta-1.0));
      if(fabs(Yhome-p->Y20[i][i])>TINYSQ)
	{
	  fprintf(logfile,"Error! Y2[%d][%d] = %0.2f, Y20[%d][%d] = %0.2f\n\n",i,i,Yhome,i,i,p->Y20[i][i]);
	  return 1;
	}
    }

  // steady state exporter dynamics
  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];
	  set_exporter_consts1(&(ev[i][ii]), 1.0, p->beta, MC[i], 1.0);
	  set_exporter_consts2(&(ev[i][ii]), p->L0[j]/p->L0[0], 1.0, p->Y20[j][i], p->Ybar2[j][i], p->theta, p->ep[i][ii].alpha, 1.0);
	  set_exporter_consts3(&(ev[i][ii]), p->L0[i]/p->L0[0], 1.0, p->Y20[i][i], p->Ybar2[i][i], p->theta);
	}

      for(ii=0; ii<2; ii++)
	{
	  if(no_exporter_dyn)
	    {
	      if(compute_moments_static_version(&(p->ep[i][ii]), &(ev[i][ii])))
		{
		  fprintf(logfile,"Error! Failed to compute steady-state moments! (i,ii)=(%d,%d)\n\n",i,ii);
		  return 1;
		}
	    }
	  else
	    {
	      discretize_phi(&(p->ep[i][ii]));

	      if(set_entry_cutoff(&(p->ep[i][ii]), &(ev[i][ii])))
		{
		  fprintf(logfile,"Error! Failed to solve for analytical entry cutoff (i,ii)=(%d,%d)\n\n",i,ii);
		  return 1;
		}

	      if(solve_steady_state_policies(&(p->ep[i][ii]), &(ev[i][ii])))
		{
		  fprintf(logfile,"Error! Failed to solve for steady-state policy function! (i,ii)=(%d,%d)\n\n",i,ii);
		  return 1;
		}

	      if(solve_steady_state_dist(&(p->ep[i][ii]), &(ev[i][ii])))
		{
		  fprintf(logfile,"Error! Failed to solve for stationary distribution! (i,ii)=(%d,%d)\n\n",i,ii);
		  return 1;
		}

	      if(full_mkt_pen==0)
		{
		  if(compute_moments(&(p->ep[i][ii]), &(ev[i][ii])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (i,ii)=(%d,%d)\n\n",i,ii);
		      return 1;
		    }
		}
	      else
		{
		  if(compute_moments_full_mkt_pen(&(p->ep[i][ii]), &(ev[i][ii])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (i,ii)=(%d,%d)\n\n",i,ii);
		      return 1;
		    }
		}
	    }
	}
    }


  // factor demand
  double exrate[NC][NC-1];
  double exportdiff[NC][NC-1];
  double t5[NC];
  double grdiff[NC];
  double Kd[NC];
  double Ld[NC];
  double Md[NC];

  for(i=0; i<NC; i++)
    {      
      if(leontief)
	{
	  Kd[i] = (p->gamma[i]/MC[i]) *
	    pow(p->R0[i]*(1.0-p->alpha)/1.0/p->alpha,p->alpha-1.0)*p->Nbar[i]*Dh[i][i]*p->ep[i][0].Zi;

	  Ld[i] = (p->R0[i]*(1.0-p->alpha)/1.0/p->alpha)*Kd[i];

	  Md[i] = ((1.0-p->gamma[i])/MC[i]) * p->Nbar[i]*Dh[i][i]*p->ep[i][0].Zi;
	}
      else
	{
	  Kd[i] = (p->gamma[i]*p->alpha/p->R0[i]) * p->Nbar[i] * Dh[i][i] * p->ep[i][0].Zi;
	  Ld[i] = (p->gamma[i]*(1.0-p->alpha)/1.0) * p->Nbar[i] * Dh[i][i] * p->ep[i][0].Zi;
	  Md[i] = ((1.0-p->gamma[i])/1.0) * p->Nbar[i] * Dh[i][i] * p->ep[i][0].Zi;
	}

      double lf = 0.0;

      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];
	  
	  exrate[i][ii] = ev[i][ii].theoretical_export_participation_rate;
	  //exrate[i][ii] = ev[i][ii].export_participation_rate;

	  if(leontief)
	    {
	      Kd[i] += (p->gamma[i]/MC[i]) *
		pow(p->R0[i]*(1.0-p->alpha)/1.0/p->alpha,p->alpha-1.0)*p->Nbar[i]*Dh[j][i]*ev[i][ii].Z;
	      
	      Ld[i] += (p->gamma[i]/MC[i]) *
		pow(p->R0[i]*(1.0-p->alpha)/1.0/p->alpha,p->alpha)*p->Nbar[i]*Dh[j][i]*ev[i][ii].Z;
	      
	      Md[i] += ((1.0-p->gamma[i])/MC[i]) * p->Nbar[i]*Dh[j][i]*ev[i][ii].Z;

	    }
	  else
	    {
	      Kd[i] += (p->gamma[i]*p->alpha/p->R0[i]) * p->Nbar[i] * Dh[j][i] * ev[i][ii].Z;
	      Ld[i] += (p->gamma[i]*(1.0-p->alpha)/1.0) * p->Nbar[i] * Dh[j][i] * ev[i][ii].Z;
	      Md[i] += ((1.0-p->gamma[i])/1.0) * p->Nbar[i] * Dh[j][i] * ev[i][ii].Z;
	    }

	  if(!allexport)
	    {
	      lf += p->Nbar[i] * ev[i][ii].lf;
	    }

	  double PYj = Db[j][i] * p->Nbar[i] * ev[i][ii].Z;
	  exportdiff[i][ii] = (PYj - p->Y20[j][i])/p->Y20[j][i];
	}

      t5[i] = ((export_participation_rate_target[i][0]*ev[i][0].top5_share + 
		export_participation_rate_target[i][1]*ev[i][1].top5_share)/
	       (export_participation_rate_target[i][0]+export_participation_rate_target[i][1]));

      grdiff[i] = ((export_participation_rate_target[i][0]*ev[i][0].rel_entrant_growth_rate + 
		    export_participation_rate_target[i][1]*ev[i][1].rel_entrant_growth_rate)/
		   (export_participation_rate_target[i][0]+export_participation_rate_target[i][1]));
      
      // recompute factors, investment, and consumption as necessary
      p->K0[i] = Kd[i];
      Ld[i] = Ld[i] + lf;
      p->L0[i] = Ld[i];
      p->Lbar[i] = 3.0*p->L0[i];
      p->X0[i] = p->delta*Kd[i];
      p->C0[i] = p->VA0[i] - sum(p->NX0[i],2) - p->X0[i];
      
      for(ii=0; ii<NC-1; ii++)
	{
	  if(no_exporter_dyn==0)
	    {
	      if(full_mkt_pen==0)
		{
		  memcpy((double *)(p->diste0[i][ii]),
			 (const double *)(ev[i][ii].diste),
			 NPHI*sizeof(double));
		  
		  memcpy((double *)(p->dist0[i][ii]),
			 (const double *)(ev[i][ii].dist),
			 NPHI*NN*sizeof(double));
		}
	      else
		{
		  memcpy((double *)(p->diste0[i][ii]),
			 (const double *)(ev[i][ii].diste),
			 NPHI*sizeof(double));
		  memcpy((double *)(p->disti0[i][ii]),
			 (const double *)(ev[i][ii].disti),
			 NPHI*sizeof(double));
		}
	    }
	}
    }

  if(!allexport)
    {
      uint cnt=0;
      for(i=0; i<NC; i++)
	{
	  myf[cnt] = (exrate[i][0] - export_participation_rate_target[i][0])/export_participation_rate_target[i][0];
	  cnt = cnt+1;

	  myf[cnt] = (exrate[i][1] - export_participation_rate_target[i][1])/export_participation_rate_target[i][1];
	  cnt = cnt+1;
	  
	  if(no_exporter_dyn==1 || full_mkt_pen==1)
	    {
	      myf[cnt] = t5[i] - top5_share_target[i];
	      cnt = cnt+1;
	    }
	  else
	    {
	      myf[cnt] = t5[i] - top5_share_target[i];
	      cnt = cnt+1;

	      if(low_dep)
		{
		  myf[cnt] = grdiff[i] - 0.3;
		}
	      else
		{
		  myf[cnt] = grdiff[i] - rel_entrant_growth_rate_target[i];
		}
	      cnt=cnt+1;
	    }
	  
	  for(ii=0; ii<2; ii++)
	    {
	      myf[cnt] = exportdiff[i][ii];
	      cnt = cnt+1;
	    }
	}
	
      for(i=0; i<cnt; i++)
	{
	  if(gsl_isnan(myf[i]) || gsl_isinf(myf[i]))
	    {
	      fprintf(logfile,KRED "NaN or Inf detected in calibration function!\n" RESET);
	      return 1;
	    }
	}
    }

  return 0;
}


int calfunc_f(const gsl_vector * x, void * data, gsl_vector * f)
{
  if(eval_exporter_moments(x->data,f->data))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int calfunc_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&calfunc_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int calfunc_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(calfunc_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&calfunc_f, x, J, 0))
	{
	  return GSL_EFAILED;
	}
      else
	{
	  return GSL_SUCCESS;
	}
    }
}

#endif
