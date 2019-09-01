#ifndef __CALIBRATE_C__
#define __CALIBRATE_C__

#include "calibrate.h"

double export_participation_rate_target[NC][2] = {{0.55,0.43},{0.0572,0.37},{0.041,0.1}};
double top5_share_target[NC] = {0.58,0.58,0.58};
double exit_rate_target[NC][NC-1] = {{0.45,0.45},{0.45,0.45},{0.45,0.45}};
double exit_rate_target2[NC][NC-1];

uint homotopy_times = 5;
//double exit_rate_target = 0.75;
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
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->sigma[i]));
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->omega[i]));
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->kappa0[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->kappa1[i][j]));
	}
    } 

  for(i=0; i<NC; i++)
    {
      cnt += 1;
      matched += fscanf(file,"%lf",&(p->Nbar[i]));
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
      for(j=0; j<2; j++)
	{
	  cnt += 1;
	  matched += fscanf(file,"%lf",&(p->n0[i][j]));
	}
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
      fprintf(file,"%0.15f\n",p->sigma[i]);
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->omega[i]);
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  fprintf(file,"%0.15f\n",p->kappa0[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  fprintf(file,"%0.15f\n",p->kappa1[i][j]);
	}
    } 

  for(i=0; i<NC; i++)
    {
      fprintf(file,"%0.15f\n",p->Nbar[i]);
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
      for(j=0; j<2; j++)
	{
	  fprintf(file,"%0.15f\n",p->n0[i][j]);
	}
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
  fprintf(file,"%d\n\n",p->cap_adj_cost);

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

  fprintf(file,"sigma:");
  for(i=0; i<NC; i++)
    {
      fprintf(file,"% 0.4f",p->sigma[i]);
    } 
  fprintf(file,"\n\n");

  fprintf(file,"omega:");
  for(i=0; i<NC; i++)
    {
      fprintf(file," %0.3f",p->omega[i]);
    } 
  fprintf(file,"\n\n");

  fprintf(file,"kappa0:");
  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  fprintf(file," %0.3f",p->kappa0[i][j]);
	}
      fprintf(file,"\n");
    } 
  fprintf(file,"\n");

  fprintf(file,"kappa1:");
  for(i=0; i<NC; i++)
    {
      for(j=0; j<2; j++)
	{
	  fprintf(file," %0.3f",p->kappa1[i][j]);
	}
      fprintf(file,"\n");
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
  memcpy((double *)(dest->sigma),(const double *)(src->sigma),sizeof(double)*NC);
  memcpy((double *)(dest->omega),(const double *)(src->omega),sizeof(double)*NC);
  memcpy((double *)(dest->kappa0),(const double *)(src->kappa0),sizeof(double)*NC*(NC-1));
  memcpy((double *)(dest->kappa1),(const double *)(src->kappa1),sizeof(double)*NC*(NC-1));
  memcpy((double *)(dest->Nbar),(const double *)(src->Nbar),sizeof(double)*NC);
  memcpy((unsigned int *)(dest->Ji),(const unsigned int *)(src->Ji),sizeof(unsigned int)*NC*(NC-1));
  
  // time series parameters
  memcpy((double *)(dest->tau_ts),(const double *)(src->tau_ts),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(dest->ntb_ts),(const double *)(src->ntb_ts),sizeof(double)*(NT+1)*NC*NC);
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
  memcpy((double *)(dest->n0),(const double *)(src->n0),sizeof(double)*NC*(NC-1));
  memcpy((double *)(dest->EX0),(const double *)(src->EX0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->IM0),(const double *)(src->IM0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->NX0),(const double *)(src->NX0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->B0),(const double *)(src->B0),sizeof(double)*NC);
  memcpy((double *)(dest->R0),(const double *)(src->R0),sizeof(double)*NC);

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
  p->rbgp = 0.02;
  p->delta = 0.06;
  p->etaK = 0.8;
#endif

  p->beta = 1.0 / (1.0 + p->rbgp);

  if(nokappa)
    {
      p->zeta = 7.2;
    }
  else if(eqkappa)
    {
      p->zeta = 3.7;
    }
  else if(zeta_sens)
    {
      p->zeta = 1.01;
    }
  else
    {
      p->zeta = 3.65;
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
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

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
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

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
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

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
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

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
  else if(scenario==6)
    {
      double ntb_mult=0.75;
      if(higher_ch_trd_costs)
	{
	  ntb_mult = ntb_mult*1.0;
	}

      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=0; t<(NT+1); t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	  if(higher_ch_trd_costs)
	    {
	      p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	      p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	    }
	  else
	    {
	      p->tau_ts[t][0][1] = 3.57662801726/100.0;
	      p->tau_ts[t][1][0] = 2.12371862378/100.0;
	    }
	}
    }
  else if(scenario==7 || scenario==10)
    {
      double ntb_mult=0.75;
      if(higher_ch_trd_costs)
	{
	  ntb_mult = ntb_mult*1.0;
	}

      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=0; t<TCH1; t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	  if(higher_ch_trd_costs)
	    {
	      p->tau_ts[t][0][1] = 27.0/100.0;
	      p->tau_ts[t][1][0] = 27.0/100.0;
	    }
	  else
	    {
	      p->tau_ts[t][0][1] = 3.57662801726/100.0;
	      p->tau_ts[t][1][0] = 2.12371862378/100.0;
	    }
	}
    }
  else if(scenario==8 || scenario==11)
    {
      double ntb_mult=0.75;
      if(higher_ch_trd_costs)
	{
	  ntb_mult = ntb_mult*1.0;
	}

      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=0; t<TCH1; t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	  if(higher_ch_trd_costs)
	    {
	      p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	      p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	    }
	  else
	    {
	      p->tau_ts[t][0][1] = 3.57662801726/100.0;
	      p->tau_ts[t][1][0] = 2.12371862378/100.0;
	    }
	}
      for(t=TCH2; t<(NT+1); t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	  if(higher_ch_trd_costs)
	    {
	      p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	      p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	    }
	  else
	    {
	      p->tau_ts[t][0][1] = 3.57662801726/100.0;
	      p->tau_ts[t][1][0] = 2.12371862378/100.0;
	    }
	}
    }
  else if(scenario>=1000 && scenario<2000)
    {
      double ntb_mult=0.75;
      if(higher_ch_trd_costs)
	{
	  ntb_mult = ntb_mult*1.0;
	}

      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=0; t<TCH1; t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	  if(higher_ch_trd_costs)
	    {
	      p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	      p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	    }
	  else
	    {
	      p->tau_ts[t][0][1] = 3.57662801726/100.0;
	      p->tau_ts[t][1][0] = 2.12371862378/100.0;
	    }
	}
      if(TCH1+1+(scenario-1000)<=TCH2)
	{
	  for(t=TCH1+1+(scenario-1000); t<(NT+1); t++)
	    {
	      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	      if(higher_ch_trd_costs)
		{
		  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
		  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
		}
	      else
		{
		  p->tau_ts[t][0][1] = 3.57662801726/100.0;
		  p->tau_ts[t][1][0] = 2.12371862378/100.0;
		}
	    }
	}
    }
  else if(scenario>=2000)
    {
      double ntb_mult=0.75;
      if(higher_ch_trd_costs)
	{
	  ntb_mult = ntb_mult*1.0;
	}

      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=0; t<TCH1; t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	  if(higher_ch_trd_costs)
	    {
	      p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	      p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	    }
	  else
	    {
	      p->tau_ts[t][0][1] = 3.57662801726/100.0;
	      p->tau_ts[t][1][0] = 2.12371862378/100.0;
	    }
	}
      if(TCH1+1+(scenario-2000)<=TCH2)
	{
	  for(t=TCH1+1+(scenario-2000); t<(NT+1); t++)
	    {
	      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	      if(higher_ch_trd_costs)
		{
		  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
		  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
		}
	      else
		{
		  p->tau_ts[t][0][1] = 3.57662801726/100.0;
		  p->tau_ts[t][1][0] = 2.12371862378/100.0;
		}
	    }
	}
    }
}

#ifdef CH_SIMPLE
void set_tariffs_ch(stoch_params_ch * sppp_ch_)
{
  uint t, ih;
  double ntb_mult = 0.75;
  if(higher_ch_trd_costs)
    {
      ntb_mult = ntb_mult*1.0;
    }

  // no-reform scenario: tariffs never fall
  ih=0;
  params * p = &((sppp_ch_->ppp)[ih]); 

  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=0; t<(NT+1); t++)
    {
      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

      if(higher_ch_trd_costs)
	{
	  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	}
      else
	{
	  p->tau_ts[t][0][1] = 3.57662801726/100.0;
	  p->tau_ts[t][1][0] = 2.12371862378/100.0;
	}
    }

  // pessimistic scenario: tariffs rise again in TCH1
  ih=1;
  p = &((sppp_ch_->ppp)[ih]);

  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=0; t<TCH1; t++)
    {
      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

      if(higher_ch_trd_costs)
	{
	  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	}
      else
	{
	  p->tau_ts[t][0][1] = 3.57662801726/100.0;
	  p->tau_ts[t][1][0] = 2.12371862378/100.0;
	}
    }
  for(t=TCH2; t<(NT+1); t++)
    {
      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

      if(higher_ch_trd_costs)
	{
	  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	}
      else
	{
	  p->tau_ts[t][0][1] = 3.57662801726/100.0;
	  p->tau_ts[t][1][0] = 2.12371862378/100.0;
	}
    }

  // optimistic scenario: tariffs stay low forever
  ih=2;
  p = &((sppp_ch_->ppp)[ih]);
  
  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=0; t<TCH1; t++)
    {
      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

      if(higher_ch_trd_costs)
	{
	  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	}
      else
	{
	  p->tau_ts[t][0][1] = 3.57662801726/100.0;
	  p->tau_ts[t][1][0] = 2.12371862378/100.0;
	}
    }
}
#else
void set_tariffs_ch(stoch_params_ch * sppp_ch_)
{
  uint t, ih;
  double ntb_mult = 0.75;
  if(higher_ch_trd_costs)
    {
      ntb_mult = ntb_mult*1.0;
    }

  // no-reform scenario: tariffs never fall
  ih=0;
  params * p = &((sppp_ch_->ppp)[ih]); 

  SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

  for(t=0; t<(NT+1); t++)
    {
      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

      if(higher_ch_trd_costs)
	{
	  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	}
      else
	{
	  p->tau_ts[t][0][1] = 3.57662801726/100.0;
	  p->tau_ts[t][1][0] = 2.12371862378/100.0;
	}
    }

  for(ih=1; ih<NHIST_CH; ih++)
    {
      p = &((sppp_ch_->ppp)[ih]);

      SET_ALL_V(p->tau_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->ntb_ts,(NT+1)*NC*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC,1.0);

      for(t=0; t<TCH1; t++)
	{
	  p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	  p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	  if(higher_ch_trd_costs)
	    {
	      p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
	      p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
	    }
	  else
	    {
	      p->tau_ts[t][0][1] = 3.57662801726/100.0;
	      p->tau_ts[t][1][0] = 2.12371862378/100.0;
	    }
	}
      if(TCH1+ih<=TCH2)
	{
	  for(t=TCH1+ih; t<(NT+1); t++)
	    {
	      p->ntb_ts[t][0][1] = ntb_mult*8.708/100.0;
	      p->ntb_ts[t][1][0] = ntb_mult*6.950/100.0;

	      if(higher_ch_trd_costs)
		{
		  p->tau_ts[t][0][1] = ch_hi*3.57662801726/100.0;
		  p->tau_ts[t][1][0] = ch_hi*2.12371862378/100.0;
		}
	      else
		{
		  p->tau_ts[t][0][1] = 3.57662801726/100.0;
		  p->tau_ts[t][1][0] = 2.12371862378/100.0;
		}
	    }
	}
    }
}
#endif

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

  // export participation:
  // - assume type-j fractions are equal to share of exports to country j
  // - assume type-j and type-k exporters have same propensity for multi-market participation
  // - non-exporter shares to 75% for all countries
  // - 65% of exporters export to one market only
  //SET_ALL_V(p->omega,NC,0.5);
  //p->omega[0] = 0.489; // UK: type-EU = 45.3%, type-RW = 54.7%
  //p->omega[1] = 0.109; // EU: type-UK = 11.8%, type-RW = 88.2%
  //p->omega[2] = 0.130; // ROW; type-UK = 12.1%, type-EU = 87.9%
  p->omega[0] = 0.5;
  p->omega[1] = 0.5;
  p->omega[2] = 0.5;

  if(nokappa)
    {
      SET_ALL_V(p->n0,NC*(NC-1),1.0);
    }
  else
    {
      uint i, ii;
      for(i=0; i<NC; i++)
	{
	  for(ii=0; ii<2; ii++)
	    {
	      p->n0[i][ii]=export_participation_rate_target[i][ii];
	    }
	}
    }

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

uint calibrate_firm_params()
{
  params * p = &(ppp0[0]);

  uint i,j,ii;
  double MC[NC];
  double D[NC][NC];
  double Db[NC][NC];
  //double Dh[NC][NC];
  double Dt[NC][NC];
  double Zi[NC];

  // normalize mass of firms to one
  SET_ALL_V(p->Nbar,NC,1.0);

  // initial guess for dispersion (George's China paper)
  if(nokappa)
    {
      SET_ALL_V(p->sigma,NC,0.25);
    }
  else
    {
      SET_ALL_V(p->sigma,NC,0.05);
    }

  // initialize some other stuff
  //SET_ALL_V(p->Ybar2,NC*NC,1.0);
  SET_ALL_V(p->kappa0,NC*(NC-1),0.0);
  SET_ALL_V(p->kappa1,NC*(NC-1),0.0);
 
  // first calculate the domestic demand stuff
  for(i=0; i<NC; i++)
    {
      // normalize number of firms to employment (so avg firm size = 1 worker in all countries)
      //p->Nbar[i] = p->L0[i];

      if(nokappa)
	{
	  Zi[i]=1.0;
	}
      else
	{
	  Zi[i] = exp(p->sigma[i]*p->sigma[i]*(p->theta-1.0)*(p->theta-1.0)/2.0);
	}

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
      
      //p->Nbar[i] = pow(p->theta/(p->theta-1.0),p->theta-1.0) * 
      //pow(MC[i],p->theta-1.0) * pow(p->Ybar2[i][i],1.0-p->theta) * (1.0/Zi[i]);
      
      p->Ybar2[i][i] = pow(p->theta/(p->theta-1.0),p->theta-1.0) * 
	pow(MC[i],p->theta-1.0) * (1.0/p->Nbar[i]) * (1.0/Zi[i]);
      p->Ybar2[i][i] = pow(p->Ybar2[i][i],1.0/(p->theta-1.0));

      D[i][i] = pow(p->Ybar2[i][i], p->theta-1.0) * p->Y20[i][i];
      Db[i][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * D[i][i] * pow(MC[i],1.0-p->theta);
      //Dh[i][i] = pow(p->theta/(p->theta-1.0),-p->theta) * D[i][i] * pow(MC[i],1.0-p->theta);
      Dt[i][i] = Db[i][i]/p->theta;
      
      double Phome = (1.0/p->Ybar2[i][i]) * pow(p->Nbar[i],1.0/(1.0-p->theta)) * 
	(p->theta/(p->theta-1.0))*MC[i] * pow(Zi[i],1.0/(1.0-p->theta));

      if(fabs(Phome-1.0)>TINYSQ)
	{
	  fprintf(logfile,"Error! P[%d][%d] = %0.2f\n\n",i,i,Phome);
	  return 1;
	}
      
      double Yhome = p->Nbar[i] * Db[i][i] * Zi[i];
      if(fabs(Yhome-p->Y20[i][i])>TINYSQ)
	{
	  fprintf(logfile,"Error! P*Y2[%d][%d] = %0.2f, P*Y20[%d][%d] = %0.2f\n\n",i,i,Yhome,i,i,p->Y20[i][i]);
	  return 1;
	}

      Yhome = p->Ybar2[i][i] * pow( MC[i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	pow(p->Nbar[i],p->theta/(p->theta-1.0)) * D[i][i] * pow(Zi[i],p->theta/(p->theta-1.0));
      if(fabs(Yhome-p->Y20[i][i])>TINYSQ)
	{
	  fprintf(logfile,"Error! Y2[%d][%d] = %0.2f, Y20[%d][%d] = %0.2f\n\n",i,i,Yhome,i,i,p->Y20[i][i]);
	  return 1;
	}
    }

  // calculate export demand stuff and fixed costs in static model where kappa0 and kappa1 are assumed to
  // be equal
  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  if(nokappa)
	    {
	      ev[i][ii].zp = -HUGE_VAL;
	      ev[i][ii].zm = -HUGE_VAL;
	      ev[i][ii].Zp = Zi[i];
	      ev[i][ii].Zm = Zi[i];
	      ev[i][ii].Z = Zi[i];
	      ev[i][ii].n = 1.0;
	    }
	  else
	    {
	      ev[i][ii].zp = p->sigma[i]*gsl_cdf_ugaussian_Pinv(1.0-p->n0[i][ii]);
	      ev[i][ii].zm = ev[i][ii].zp;
	      
	      ev[i][ii].Zp = Zi[i] * 
		gsl_cdf_ugaussian_P( (p->sigma[i]*p->sigma[i]*(p->theta-1.0) - ev[i][ii].zp)/p->sigma[i] );
	      
	      ev[i][ii].Zm = Zi[i] * 
		gsl_cdf_ugaussian_P( (p->sigma[i]*p->sigma[i]*(p->theta-1.0) - ev[i][ii].zm)/p->sigma[i]);
	      
	      ev[i][ii].Z = p->n0[i][ii] * ev[i][ii].Zm + (1.0-p->n0[i][ii]) * ev[i][ii].Zp;
	    }

	  double om = p->omega[i];
	  if(ii>0)
	    {
	      om = 1.0-p->omega[i];
	    }

	  p->Ybar2[j][i] = pow(p->theta/(p->theta-1.0),p->theta-1.0) * 
	    pow(MC[i],p->theta-1.0) * (1.0/(p->Nbar[i]*om)) * (1.0/ev[i][ii].Z);
	  p->Ybar2[j][i] = pow(p->Ybar2[j][i],1.0/(p->theta-1.0));
	  
	  D[j][i] = pow(p->Ybar2[j][i], p->theta-1.0) * p->Y20[j][i];
	  Db[j][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * D[j][i] * pow(MC[i],1.0-p->theta);
	  //Dh[j][i] = pow(p->theta/(p->theta-1.0),-p->theta) * D[j][i] * pow(MC[i],1.0-p->theta);
	  Dt[j][i] = Db[j][i]/p->theta;
      
	  double Yex = p->Nbar[i] * om * Db[j][i] * ev[i][ii].Z;
	  if(fabs(Yex-p->Y20[j][i])>TINYSQ)
	    {
	      fprintf(logfile,"Error! P*Y2[%d][%d] = %0.2f, P*Y20[%d][%d] = %0.2f\n\n",j,i,Yex,j,i,p->Y20[j][i]);
	      return 1;
	    }

	  Yex = p->Ybar2[j][i] * pow( MC[i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	    pow(p->Nbar[i]*om,p->theta/(p->theta-1.0)) * D[j][i] * pow(ev[i][ii].Z,p->theta/(p->theta-1.0));
	  if(fabs(Yex-p->Y20[j][i])>TINYSQ)
	    {
	      fprintf(logfile,"Error! Y2[%d][%d] = %0.2f, Y20[%d][%d] = %0.2f\n\n",j,i,Yex,j,i,p->Y20[j][i]);
	      return 1;
	    }

	  if(nokappa)
	    {
	      p->kappa0[i][ii]=0.0;
	      p->kappa1[i][ii]=0.0;
	    }
	  else
	    {
	      p->kappa0[i][ii] = Dt[j][i] * exp((p->theta-1.0)*ev[i][ii].zp);
	      p->kappa1[i][ii] = p->kappa0[i][ii];

	      //p->kappa0[i][ii] = p->kappa0[i][ii]*1.05;
	      //p->kappa1[i][ii] = p->kappa1[i][ii]/1.05;
	    }
	}
    }

  solver_n = 21;
  alloc_solver_mem();
  uint status = 0;

  solver_x->data[0] = log(p->sigma[0]);
  solver_x->data[1] = log(p->sigma[1]);
  solver_x->data[2] = log(p->sigma[2]);
  solver_x->data[3] = log(p->Ybar2[0][1]);
  solver_x->data[4] = log(p->Ybar2[0][2]);
  solver_x->data[5] = log(p->Ybar2[1][0]);
  solver_x->data[6] = log(p->Ybar2[1][2]);
  solver_x->data[7] = log(p->Ybar2[2][0]);
  solver_x->data[8] = log(p->Ybar2[2][1]);
  solver_x->data[9] = log(p->kappa0[0][0]);
  solver_x->data[10] = log(p->kappa0[0][1]);
  solver_x->data[11] = log(p->kappa0[1][0]);
  solver_x->data[12] = log(p->kappa0[1][1]);
  solver_x->data[13] = log(p->kappa0[2][0]);
  solver_x->data[14] = log(p->kappa0[2][1]);
  solver_x->data[15] = log(p->kappa1[0][0]);
  solver_x->data[16] = log(p->kappa1[0][1]);
  solver_x->data[17] = log(p->kappa1[1][0]);
  solver_x->data[18] = log(p->kappa1[1][1]);
  solver_x->data[19] = log(p->kappa1[2][0]);
  solver_x->data[20] = log(p->kappa1[2][1]);

  if(eval_calfn_once || nokappa)
    {
      status = calfunc_f(solver_x,NULL,f0[0]);
    }
  else
    {
      solver_verbose=0;
      par=0;
      gsl_multiroot_function_fdf f = {&calfunc_f,&calfunc_df,&calfunc_fdf,21,NULL};
      
      if(eqkappa)
	{
	  fprintf(logfile,KBLU "Calibrating firm parameters...\n" RESET);
	  status = find_root_deriv_mkl(&f);
	  if(status)
	    {
	      fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
	      status=1;
	    }
	  status = calfunc_f(solver_x,NULL,f0[0]);
	}
      else
	{
	  eqkappa=1;
	  fprintf(logfile,KBLU "One-shot calibration for kappa0=kappa1...\n" RESET);
	  status = find_root_deriv_mkl(&f);
	  if(status)
	    {
	      fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
	      status=1;
	    }
	  else
	    {
	      eqkappa=0;
	      for(i=0; i<NC; i++)
		{
		  for(ii=0; ii<2; ii++)
		    {
		      exit_rate_target2[i][ii] = 1.0-export_participation_rate_target[i][ii];
		    }
		}

	      fprintf(logfile,KBLU "\nHomotopy process for kappa0!=kappa1 beginning...\n\n" RESET);

	      double grid[NC][NC-1][homotopy_times];
	      for(i=0; i<NC; i++)
		{
		  for(ii=0; ii<2; ii++)
		    {
		      linspace(exit_rate_target[i][ii],
			       exit_rate_target2[i][ii],
			       homotopy_times,
			       grid[i][ii]);
		    }
		}
	      int h;
	      for(h=homotopy_times-1;h>=0;h--)
		{
		  fprintf(logfile,"\tExit rate targets:");
		  for(i=0; i<NC; i++)
		    {
		      fprintf(logfile," (");
		      for(ii=0; ii<2; ii++)
			{
			  exit_rate_target2[i][ii] = grid[i][ii][h];
			  fprintf(logfile," %0.5f", exit_rate_target2[i][ii]);
			}
		      fprintf(logfile," )");
		    }
		  status = find_root_deriv_mkl(&f);
		  if(status)
		    {
		      fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
		      status=1;
		      break;
		    }
		  fprintf(logfile,"\tModel exit rates:");
		  for(i=0; i<NC; i++)
		    {
		      fprintf(logfile,"(");
		      for(ii=0; ii<2; ii++)
			{
			  fprintf(logfile," %0.5f",ev[i][ii].Fzm);
			}
		      fprintf(logfile,")");
		    }
		  fprintf(logfile,"\n\n");
		}
	    }
	}

    }
  solver_verbose=1;
  
  free_solver_mem();
 
  return status;
}

uint calibrate()
{
  if(recal==1 || nokappa || eqkappa || psi_sens || zeta_sens)
    {  
      set_nontargeted_params();

      if(load_iomat())
	{
	  fprintf(logfile, KRED "\nFailed to load iomat!\n" RESET);
	  return 1;
	}
      if(store_base_period_values())
	{
	  fprintf(logfile, KRED "\nFailed to store base period values!\n" RESET);
	  return 1;
	}

      if(calibrate_agg_params())
	{
	  fprintf(logfile, KRED "\nFailed to calibrate agg params!\n" RESET);
	  return 1;
	}

      if(calibrate_firm_params())
	{
	  fprintf(logfile, KRED "\nFailed to calibrate firm params!\n" RESET);
	  return 1;
	}

      if(calibrate_hh_params())
	{
	  fprintf(logfile, KRED "\nFailed to calibrate hh params!\n" RESET);
	  return 1;
	}

      if(nokappa==0 && eqkappa==0 && psi_sens==0 && zeta_sens==0)
	{
	  if(write_params_to_file())
	    {
	      fprintf(logfile, KRED "\nFailed to write parameters to file!\n" RESET);
	      return 1;
	    }
	}
    }
  else
    {
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
      if(copy_params(&(ppp0_ch[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp0_ch!\n" RESET);
	  return 1;
	}

#ifdef CH_SIMPLE
      if(copy_params(&(ppp1_ch[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp1_ch!\n" RESET);
	  return 1;
	}
      if(copy_params(&(ppp2_ch[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp1_ch!\n" RESET);
	  return 1;
	}
#else
      uint ix;
      for(ix=0; ix<NHIST_CH-1; ix++)
	{
	  if(copy_params(&(ppp1_ch[ix][it]),&(ppp0[0])))
	    {
	      fprintf(logfile, KRED "\nFailed to copy ppp1_ch!\n" RESET);
	      return 1;
	    }	  
	}
#endif

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

      stoch_params_ch * sp_ch = &(sppp_ch[it]); 
      for(ih=0; ih<NHIST_CH; ih++)
	{
	  if(copy_params(&((sp_ch->ppp)[ih]),&(ppp0_ch[0])))
	    {
	      fprintf(logfile, KRED "\nFailed to copy stoch_params_ch!\n" RESET);
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

  p->sigma[0] = exp(calpars[0]);
  p->sigma[1] = exp(calpars[1]);
  p->sigma[2] = exp(calpars[2]);
  p->Ybar2[0][1] = exp(calpars[3]);
  p->Ybar2[0][2] = exp(calpars[4]);
  p->Ybar2[1][0] = exp(calpars[5]);
  p->Ybar2[1][2] = exp(calpars[6]);
  p->Ybar2[2][0] = exp(calpars[7]);
  p->Ybar2[2][1] = exp(calpars[8]);  
  p->kappa0[0][0] = exp(calpars[9]);
  p->kappa0[0][1] = exp(calpars[10]);
  p->kappa0[1][0] = exp(calpars[11]);
  p->kappa0[1][1] = exp(calpars[12]);
  p->kappa0[2][0] = exp(calpars[13]);
  p->kappa0[2][1] = exp(calpars[14]);
  p->kappa1[0][0] = exp(calpars[15]);
  p->kappa1[0][1] = exp(calpars[16]);
  p->kappa1[1][0] = exp(calpars[17]);
  p->kappa1[1][1] = exp(calpars[18]);
  p->kappa1[2][0] = exp(calpars[19]);
  p->kappa1[2][1] = exp(calpars[20]);

  uint i, ii, j;
  double MC[NC];
  double D[NC][NC];
  double Db[NC][NC];
  double Dh[NC][NC];
  double Dt[NC][NC];

  // recalculate export demand scale factors
  // be equal
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
	  Dt[j][i] = Db[j][i]/p->theta;
      	}
    }

  for(i=0; i<NC; i++)
    {
      double Zi=0.0;
      if(nokappa)
	{
	  Zi=1.0;
	}
      else
	{
	  Zi = exp(p->sigma[i]*p->sigma[i]*(p->theta-1.0)*(p->theta-1.0)/2.0);
	}

      p->Ybar2[i][i] = pow(p->theta/(p->theta-1.0),p->theta-1.0) * 
	pow(MC[i],p->theta-1.0) * (1.0/p->Nbar[i]) * (1.0/Zi);
      p->Ybar2[i][i] = pow(p->Ybar2[i][i],1.0/(p->theta-1.0));

      D[i][i] = pow(p->Ybar2[i][i], p->theta-1.0) * p->Y20[i][i];
      Db[i][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * D[i][i] * pow(MC[i],1.0-p->theta);
      Dt[i][i] = Db[i][i]/p->theta;
      
      double Phome = (1.0/p->Ybar2[i][i]) * pow(p->Nbar[i],1.0/(1.0-p->theta)) * 
	(p->theta/(p->theta-1.0))*MC[i] * pow(Zi,1.0/(1.0-p->theta));

      if(fabs(Phome-1.0)>TINYSQ)
	{
	  fprintf(logfile,"Error! P[%d][%d] = %0.2f\n\n",i,i,Phome);
	  return 1;
	}
      
      double Yhome = p->Nbar[i] * Db[i][i] * Zi;
      if(fabs(Yhome-p->Y20[i][i])>TINYSQ)
	{
	  fprintf(logfile,"Error! P*Y2[%d][%d] = %0.2f, P*Y20[%d][%d] = %0.2f\n\n",i,i,Yhome,i,i,p->Y20[i][i]);
	  return 1;
	}

      Yhome = p->Ybar2[i][i] * pow( MC[i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	pow(p->Nbar[i],p->theta/(p->theta-1.0)) * D[i][i] * pow(Zi,p->theta/(p->theta-1.0));
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

	  if(ev_steady_state(p->sigma[i],
			     p->theta,
			     p->kappa0[i][ii],
			     p->kappa1[i][ii],
			     1.0,
			     p->beta,
			     Dt[j][i],
			     p->n0[i][ii],
			     &(ev[i][ii])))
	    {
	      fprintf(logfile, KRED "Failed to solve for steady state exporter moments for country/dest %d/%d!\n" RESET,i,j);
	      return 1;
	    }
	}
   
    }

  // factor demand
  double exrate[NC][NC-1];
  double exitrate[NC][NC-1];
  double exportdiff[NC][NC-1];
  double ratio[NC];
  double Kd[NC];
  double Ld[NC];
  double Md[NC];

  for(i=0; i<NC; i++)
    {      
      if(leontief)
	{
	  Kd[i] = (p->gamma[i]/MC[i]) *
	    pow(p->R0[i]*(1.0-p->alpha)/1.0/p->alpha,p->alpha-1.0)*p->Nbar[i]*Dh[i][i]*ev[i][0].Zi;

	  Ld[i] = (p->R0[i]*(1.0-p->alpha)/1.0/p->alpha)*Kd[i];

	  Md[i] = ((1.0-p->gamma[i])/MC[i]) * p->Nbar[i]*Dh[i][i]*ev[i][0].Zi;
	}
      else
	{
	  Kd[i] = (p->gamma[i]*p->alpha/p->R0[i]) * p->Nbar[i] * Dh[i][i] * ev[i][0].Zi;
	  Ld[i] = (p->gamma[i]*(1.0-p->alpha)/1.0) * p->Nbar[i] * Dh[i][i] * ev[i][0].Zi;
	  Md[i] = ((1.0-p->gamma[i])/1.0) * p->Nbar[i] * Dh[i][i] * ev[i][0].Zi;
	}

      double ynon = 0.0;
      double yex = 0.0;
      double lf = 0.0;

      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  exrate[i][ii] = ev[i][ii].n;
	  exitrate[i][ii] = ev[i][ii].Fzm;

	  double om = p->omega[i];
	  if(ii>0)
	    {
	      om = 1.0-p->omega[i];
	    }

	  if(leontief)
	    {
	      Kd[i] += (p->gamma[i]/MC[i]) *
		pow(p->R0[i]*(1.0-p->alpha)/1.0/p->alpha,p->alpha-1.0)*p->Nbar[i]*om*Dh[j][i]*ev[i][ii].Z;
	      
	      Ld[i] += (p->gamma[i]/MC[i]) *
		pow(p->R0[i]*(1.0-p->alpha)/1.0/p->alpha,p->alpha)*p->Nbar[i]*om*Dh[j][i]*ev[i][ii].Z;
	      
	      Md[i] += ((1.0-p->gamma[i])/MC[i]) * p->Nbar[i]*om*Dh[j][i]*ev[i][ii].Z;

	    }
	  else
	    {
	      Kd[i] += (p->gamma[i]*p->alpha/p->R0[i]) * p->Nbar[i] * om * Dh[j][i] * ev[i][ii].Z;
	      Ld[i] += (p->gamma[i]*(1.0-p->alpha)/1.0) * p->Nbar[i] * om * Dh[j][i] * ev[i][ii].Z;
	      Md[i] += ((1.0-p->gamma[i])/1.0) * p->Nbar[i] * Dh[j][i] * om * ev[i][ii].Z;
	}

	  if(!nokappa)
	    {
	      lf += p->Nbar[i] * om * (
				       p->kappa0[i][ii]*(1.0-ev[i][ii].n)*(1.0-ev[i][ii].Fzp) + 
				       p->kappa1[i][ii]*ev[i][ii].n*(1.0-ev[i][ii].Fzm)
				       );
	    }

	  double PYj = Db[j][i] * p->Nbar[i] * om * ev[i][ii].Z;
	  exportdiff[i][ii] = PYj - p->Y20[j][i];

	  double Znon1 = ev[i][ii].Zi * 
	    gsl_cdf_ugaussian_P((ev[i][ii].zm - p->sigma[i]*p->sigma[i]*(p->theta-1.0)) / p->sigma[i]);

	  double Znon2 = ev[i][ii].Zi * 
	    gsl_cdf_ugaussian_P((ev[i][ii].zp - p->sigma[i]*p->sigma[i]*(p->theta-1.0)) / p->sigma[i]);

	  ynon = ynon + om * 
	    Db[i][i]* (ev[i][ii].n*Znon1 + (1.0-ev[i][ii].n)*Znon2) / (1.0-ev[i][ii].n);

	  yex = yex + om * 
	    (Db[i][i] + Db[j][i]) * ev[i][ii].Z / ev[i][ii].n;
	}

      ratio[i] = yex/ynon;

      // recompute factors, investment, and consumption as necessary
      p->K0[i] = Kd[i];
      Ld[i] = Ld[i] + lf;
      p->L0[i] = Ld[i];
      p->Lbar[i] = 3.0*p->L0[i];
      p->X0[i] = p->delta*Kd[i];
      p->C0[i] = p->VA0[i] - sum(p->NX0[i],2) - p->X0[i];
    }

  if(!nokappa)
    {
      uint cnt=0;
      for(i=0; i<NC; i++)
	{
	  for(ii=0; ii<2; ii++)
	    {
	      myf[cnt] = (exrate[i][ii] - export_participation_rate_target[i][ii]);
	      cnt = cnt+1;
	    }

	  for(ii=0; ii<2; ii++)
	    {
	      if(eqkappa)
		{
		  myf[cnt] = p->kappa0[i][ii]-p->kappa1[i][ii];
		  cnt=cnt+1;
		}
	      else
		{
		  myf[cnt] = (exitrate[i][ii] - exit_rate_target2[i][ii]);
		  cnt = cnt+1;
		}
	    }

	  for(ii=0; ii<2; ii++)
	    {
	      myf[cnt] = exportdiff[i][ii];
	      cnt = cnt+1;
	    }

	  myf[cnt] = (ratio[i] - 2.5);
	  cnt = cnt+1;
	}
  
      for(i=0; i<cnt; i++)
	{
	  if(gsl_isnan(myf[i]) || gsl_isinf(myf[i]))
	    {
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
