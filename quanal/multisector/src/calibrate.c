#ifndef __CALIBRATE_C__
#define __CALIBRATE_C__

#include "calibrate.h"

double soft_tfp_loss = 0.025;
double hard_tfp_loss = 0.05;

void set_nontargeted_params(params * p)
{
  // parameters common across countries

#ifdef _QUARTERLY
  p->rbgp = 0.02/4.0; // quarterly!!
  p->delta = 0.06/4.0;
  p->etaK = 0.9;
  if(was_flag)
    {
      p->etaM = 0.0;
      p->etaF = 0.0;
    }
  else
    {
      p->etaM = 50.0;
      p->etaF = 50.0;
    }
#else
  p->rbgp = 0.02;
  p->delta = 0.06;
  p->etaK = 0.8;
  if(was_flag)
    {
      p->etaM = 0.0;
      p->etaF = 0.0;
    }
  else
    {
      p->etaM = 200.0;
      p->etaF = 200.0;
    }
#endif

  p->gbgp = 1.00;

  p->rho = 1.0-1.0/0.65;

  p->psi = -1.0;
  if(sensitivity==3)
    {
      p->psi = -4.0;
    }    

  p->cap_adj_cost = 1;

  // elasticities differ across countries...
  // one option is to use KRS international-macro elasticities
  if(sensitivity==4)
    {
      uint i;
      for(i=0; i<NC; i++)
	{
	  p->zeta[i][0] = 1.0-1.0/1.01;
	  p->zeta[i][1] = 1.0-1.0/1.01;
	  p->sig[i][0] = 1.0-1.0/1.01;
	  p->sig[i][1] = 1.0-1.0/1.01;
	}
    }
  // otherwise calculated using Caliendo and Parro (2014) and Costinaut and Rodriguez-Clare (2013)
  // see ../python/elasticities.py
  else
    {
      // UK
      p->zeta[0][0] = 1.0-1.0/7.57693;
      p->zeta[0][1] = 1.0-1.0/5.0;
      p->sig[0][0] = 1.0-1.0/4.8244;
      p->sig[0][1] = 1.0-1.0/5.0;

      // EU
      p->zeta[1][0] = 1.0-1.0/7.4797;
      p->zeta[1][1] = 1.0-1.0/5.0;
      p->sig[1][0] = 1.0-1.0/4.43168;
      p->sig[1][1] = 1.0-1.0/5.0;

      // RW
      p->zeta[2][0] = 1.0-1.0/6.58628;
      p->zeta[2][1] = 1.0-1.0/5.0;
      p->sig[2][0] = 1.0-1.0/5.261687;
      p->sig[2][1] = 1.0-1.0/5.0;
    } 

  SET_ALL_V(p->r0,NC,p->rbgp);
  SET_ALL_V(p->tauk,NC,0.25);
  SET_ALL_V(p->alpha,NC*NS,0.34);
  SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);

  return;
}

void set_tariffs(stoch_params * sppp)
{
  uint t, ih;

  // optimistic scenario: no goods tariffs, small non-tariff barriers
  // EU-USA NTBs calculated using same apprroach is in Dhingra et al. (2016)
  // "small" defined as 50% of US-EU NTBs
  ih=1;
  params * p = &((sppp->ppp)[ih]);
  double ntb_mult = 0.25;
  
  SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC*NS,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      p->ntb_m_ts[t][0][0][1] = ntb_mult*7.0492317424/100.0;
      p->ntb_m_ts[t][0][1][1] = ntb_mult*3.91656354699/100.0;
      p->ntb_f_ts[t][0][0][1] = ntb_mult*12.3159382657/100.0;
      p->ntb_f_ts[t][0][1][1] = ntb_mult*1.5096120818/100.0;

      p->ntb_m_ts[t][1][0][0] = ntb_mult*5.96268004074/100.0;
      p->ntb_m_ts[t][1][1][0] = ntb_mult*5.7548322111/100.0;
      p->ntb_f_ts[t][1][0][0] = ntb_mult*10.4688170164/100.0;
      p->ntb_f_ts[t][1][1][0] = ntb_mult*4.24505962863/100.0;

      if(tfp_flag==1)
	{
	  uint i,s;
	  for(i=0; i<NC; i++)
	    {
	      for(s=0; s<NS; s++)
		{
		  p->a_ts[t][i][s] = 1.0-soft_tfp_loss;
		}
	    }
	}
    }

  // pessimistic scenario: tariffs calculated from EU MFN tariffs as in Dhingra et al. (2016), plus
  // 100% of EU-USA NTBs from above
  ih=2;
  p = &((sppp->ppp)[ih]);
  ntb_mult = 0.75;

  SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->a_ts,(NT+1)*NC*NS,1.0);

  for(t=TBREXIT; t<(NT+1); t++)
    {
      p->ntb_m_ts[t][0][0][1] = ntb_mult*7.0492317424/100.0;
      p->ntb_m_ts[t][0][1][1] = ntb_mult*3.91656354699/100.0;
      p->ntb_f_ts[t][0][0][1] = ntb_mult*12.3159382657/100.0;
      p->ntb_f_ts[t][0][1][1] = ntb_mult*1.5096120818/100.0;

      p->ntb_m_ts[t][1][0][0] = ntb_mult*5.96268004074/100.0;
      p->ntb_m_ts[t][1][1][0] = ntb_mult*5.7548322111/100.0;
      p->ntb_f_ts[t][1][0][0] = ntb_mult*10.4688170164/100.0;
      p->ntb_f_ts[t][1][1][0] = ntb_mult*4.24505962863/100.0;

      p->tau_m_ts[t][0][0][1] = 4.23013777143/100.0;
      p->tau_f_ts[t][0][0][1] = 4.23013777143/100.0;

      p->tau_m_ts[t][1][0][0] = 3.2935061316/100.0;
      p->tau_f_ts[t][1][0][0] = 3.2935061316/100.0;

      if(tfp_flag==1)
	{
	  uint i,s;
	  for(i=0; i<NC; i++)
	    {
	      for(s=0; s<NS; s++)
		{
		  p->a_ts[t][i][s] = 1.0-hard_tfp_loss;
		}
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
      SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC*NS,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  p->ntb_m_ts[t][0][0][1] = ntb_mult*7.0492317424/100.0;
	  p->ntb_m_ts[t][0][1][1] = ntb_mult*3.91656354699/100.0;
	  p->ntb_f_ts[t][0][0][1] = ntb_mult*12.3159382657/100.0;
	  p->ntb_f_ts[t][0][1][1] = ntb_mult*1.5096120818/100.0;

	  p->ntb_m_ts[t][1][0][0] = ntb_mult*5.96268004074/100.0;
	  p->ntb_m_ts[t][1][1][0] = ntb_mult*5.7548322111/100.0;
	  p->ntb_f_ts[t][1][0][0] = ntb_mult*10.4688170164/100.0;
	  p->ntb_f_ts[t][1][1][0] = ntb_mult*4.24505962863/100.0;

	  if(tfp_flag==1)
	    {
	      uint i,s;
	      for(i=0; i<NC; i++)
		{
		  for(s=0; s<NS; s++)
		    {
		      p->a_ts[t][i][s] = 1.0-soft_tfp_loss;
		    }
		}
	    }
	}
    }
  else if(scenario==3 || scenario==5)
    {
      double ntb_mult=0.75;

      SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC*NS,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  p->ntb_m_ts[t][0][0][1] = ntb_mult*7.0492317424/100.0;
	  p->ntb_m_ts[t][0][1][1] = ntb_mult*3.91656354699/100.0;
	  p->ntb_f_ts[t][0][0][1] = ntb_mult*12.3159382657/100.0;
	  p->ntb_f_ts[t][0][1][1] = ntb_mult*1.5096120818/100.0;

	  p->ntb_m_ts[t][1][0][0] = ntb_mult*5.96268004074/100.0;
	  p->ntb_m_ts[t][1][1][0] = ntb_mult*5.7548322111/100.0;
	  p->ntb_f_ts[t][1][0][0] = ntb_mult*10.4688170164/100.0;
	  p->ntb_f_ts[t][1][1][0] = ntb_mult*4.24505962863/100.0;

	  p->tau_m_ts[t][0][0][1] = 4.23013777143/100.0;
	  p->tau_f_ts[t][0][0][1] = 4.23013777143/100.0;

	  p->tau_m_ts[t][1][0][0] = 3.2935061316/100.0;
	  p->tau_f_ts[t][1][0][0] = 3.2935061316/100.0;

	  if(tfp_flag==1)
	    {
	      uint i,s;
	      for(i=0; i<NC; i++)
		{
		  for(s=0; s<NS; s++)
		    {
		      p->a_ts[t][i][s] = 1.0-hard_tfp_loss;
		    }
		}
	    }
	}
    }
  else if(0)
    {
      double ntb_mult=0.75;

      SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC*NS,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  p->ntb_m_ts[t][0][0][1] = ntb_mult*7.0492317424/100.0;
	  p->ntb_m_ts[t][0][1][1] = ntb_mult*3.91656354699/100.0;
	  p->ntb_f_ts[t][0][0][1] = ntb_mult*12.3159382657/100.0;
	  p->ntb_f_ts[t][0][1][1] = ntb_mult*1.5096120818/100.0;

	  p->ntb_m_ts[t][1][0][0] = ntb_mult*5.96268004074/100.0;
	  p->ntb_m_ts[t][1][1][0] = ntb_mult*5.7548322111/100.0;
	  p->ntb_f_ts[t][1][0][0] = ntb_mult*10.4688170164/100.0;
	  p->ntb_f_ts[t][1][1][0] = ntb_mult*4.24505962863/100.0;

	  p->tau_m_ts[t][0][0][1] = (1.0-pi_brexit)*4.23013777143/100.0;
	  p->tau_f_ts[t][0][0][1] = (1.0-pi_brexit)*4.23013777143/100.0;

	  p->tau_m_ts[t][1][0][0] = (1.0-pi_brexit)*3.2935061316/100.0;
	  p->tau_f_ts[t][1][0][0] = (1.0-pi_brexit)*3.2935061316/100.0;

	  if(tfp_flag==1)
	    {
	      uint i,s;
	      for(i=0; i<NC; i++)
		{
		  for(s=0; s<NS; s++)
		    {
		      p->a_ts[t][i][s] = 1.0 - (pi_brexit*soft_tfp_loss+(1.0-pi_brexit)*hard_tfp_loss);
		    }
		}
	    }
	}
    }
  else if(0)
    {
      double ntb_mult=0.75;

      SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);
      SET_ALL_V(p->a_ts,(NT+1)*NC*NS,1.0);

      for(t=TBREXIT; t<(NT+1); t++)
	{
	  p->ntb_m_ts[t][0][0][1] = (1.0-pi_vote)*ntb_mult*7.0492317424/100.0;
	  p->ntb_m_ts[t][0][1][1] = (1.0-pi_vote)*ntb_mult*3.91656354699/100.0;
	  p->ntb_f_ts[t][0][0][1] = (1.0-pi_vote)*ntb_mult*12.3159382657/100.0;
	  p->ntb_f_ts[t][0][1][1] = (1.0-pi_vote)*ntb_mult*1.5096120818/100.0;

	  p->ntb_m_ts[t][1][0][0] = (1.0-pi_vote)*ntb_mult*5.96268004074/100.0;
	  p->ntb_m_ts[t][1][1][0] = (1.0-pi_vote)*ntb_mult*5.7548322111/100.0;
	  p->ntb_f_ts[t][1][0][0] = (1.0-pi_vote)*ntb_mult*10.4688170164/100.0;
	  p->ntb_f_ts[t][1][1][0] = (1.0-pi_vote)*ntb_mult*4.24505962863/100.0;

	  p->tau_m_ts[t][0][0][1] = (1.0-pi_vote)*(1.0-pi_brexit)*4.23013777143/100.0;
	  p->tau_f_ts[t][0][0][1] = (1.0-pi_vote)*(1.0-pi_brexit)*4.23013777143/100.0;

	  p->tau_m_ts[t][1][0][0] = (1.0-pi_vote)*(1.0-pi_brexit)*3.2935061316/100.0;
	  p->tau_f_ts[t][1][0][0] = (1.0-pi_vote)*(1.0-pi_brexit)*3.2935061316/100.0;

	  if(tfp_flag==1)
	    {
	      uint i,s;
	      for(i=0; i<NC; i++)
		{
		  for(s=0; s<NS; s++)
		    {
		      p->a_ts[t][i][s] = 1.0 - (1.0-pi_vote)*(pi_brexit*soft_tfp_loss+(1.0-pi_brexit)*hard_tfp_loss);
		    }
		}
	    }
	}
    }
}

uint load_iomat(params * p)
{
  uint i, j, got;
  double tmp;
  FILE * file = fopen("input/iomat2011.txt","rb");
  if(file)
    {
      for(i=0; i<(NS*NC+2); i++)
	{
	  for(j=0; j<(NS*NC+NF*NC+1); j++)
	    {
	      got = fscanf(file,"%lf",&tmp);
	      if(got != 1)
		{
		  fprintf(logfile,KRED "Error reading IO matrix!\n" RESET);
		  fclose(file);
		  return 1;
		}
#ifdef _QUARTERLY
	      p->iomat[i][j] = tmp/4.0; // quarterly!!
#else
	      p->iomat[i][j] = tmp;
#endif
	      
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

void load_ts_params(params * p)
{

  SET_ALL_V(p->popt_ts,NC*(NT+1),1.0);
  SET_ALL_V(p->popw_ts,NC*(NT+1),1.0);
  SET_ALL_V(p->pope_ts,NC*(NT+1),1.0);

  uint i, s, t;
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  p->gam_ts[0][i][s] = 1.0;
	  p->a_ts[0][i][s] = 1.0;
	  for(t=0; t<NT; t++)
	    {
	      p->gam_ts[t+1][i][s] = p->gam_ts[t][i][s] * p->gbgp;
	      p->a_ts[t+1][i][s] = 1.0;
	    }
	}
    }
}

uint store_base_period_values(params * p)
{
  double mkt_clear_tol = 1.0e-7;
  uint varow = NC*NS;
  uint gorow = NC*NS+1;


  SET_ALL_V(p->y0,NC*NS,0.0);
  SET_ALL_V(p->va0,NC*NS,0.0);
  SET_ALL_V(p->k0,NC*NS,0.0);
  SET_ALL_V(p->l0,NC*NS,0.0);
  SET_ALL_V(p->md0,NC*NS*NS,0.0);
  SET_ALL_V(p->m0,NC*NS,0.0);
  SET_ALL_V(p->m02,NC*NS*NC,0.0);
  SET_ALL_V(p->q0,NC*NS,0.0);
  SET_ALL_V(p->q02,NC*NS*NC,0.0);
  SET_ALL_V(p->ex0,NC*NC,0.0);
  SET_ALL_V(p->im0,NC*NC,0.0);
  SET_ALL_V(p->nx0,NC*NC,0.0);
  SET_ALL_V(p->c0,NC*NS,0.0);
  SET_ALL_V(p->i0,NC*NS,0.0);
  SET_ALL_V(p->ii0,NC,0.0);

  uint i, s, j, r;

  //double uk_k_gdp = 4.4151793399;
  //double eu_k_gdp = 4.54658836978;
  //double rw_k_gdp = 3.20495385658;

  //p->kk0[0] = 100.0*uk_k_gdp;
  //p->kk0[1] = 100.0*eu_k_gdp * (p->va0[1][0]+p->va0[1][1])/(p->va0[0][0]+p->va0[0][1]);
  //p->kk0[2] = 100.0*rw_k_gdp * (p->va0[2][0]+p->va0[2][1])/(p->va0[0][0]+p->va0[0][1]);

  for(i=0; i<NC; i++)
    {
      uint ccol = NC*NS+i;
      uint icol = NC*NS+NC+i;

      //double rdky = p->alpha[i][0]*(p->va0[i][0]+p->va0[i][1]);
      //double dky = p->delta*p->kk0[i];
      //p->tauk[i] = 1.0 - ( (dky+p->r0[i]*p->kk0[i])/rdky );
      //double rky = rdky - dky;
      //p->r0[i] = ((1.0-p->tauk[i])*rdky-dky)/p->kk0[i];


      for(s=0; s<NS; s++)
	{  
	  // first get value added and factors
	  uint scol = i*NS + s;
	  p->y0[i][s] = p->iomat[gorow][scol];
	  p->va0[i][s] = p->iomat[varow][scol];

	  p->l0[i][s] = (1.0 - p->alpha[i][s]) * p->va0[i][s];
	  p->k0[i][s] = p->alpha[i][s] * p->va0[i][s] / ((p->r0[i] + p->delta) / (1.0 - p->tauk[i]));

	  // now get demand for products from different source countries and sectors 
	  for(j=0; j<NC; j++)
	    {
	      p->c0[i][s] = p->c0[i][s] + p->iomat[j*NS+s][ccol];
	      p->i0[i][s] = p->i0[i][s] + p->iomat[j*NS+s][icol];
	      p->q02[i][s][j] = p->iomat[j*NS+s][ccol] + p->iomat[j*NS+s][icol];

	      for(r=0; r<NS; r++)
		{
		  uint rcol = i*NS + r;
		  p->m02[i][s][j] = p->m02[i][s][j] + p->iomat[j*NS+s][rcol];
		  p->md0[i][r][s] = p->md0[i][r][s] + p->iomat[j*NS+s][rcol];
		}
	    }
	  p->q0[i][s] = sum(p->q02[i][s],NC);
	  p->m0[i][s] = sum(p->m02[i][s],NC);
	}
      p->ll0[i] = sum(p->l0[i],NS);
      p->kk0[i] = sum(p->k0[i],NS);
      p->ii0[i] = sum(p->i0[i],NS);
    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(j != i)
		{
		  p->im02[i][s][j] = p->q02[i][s][j] + p->m02[i][s][j];
		  p->ex02[i][s][j] = p->q02[j][s][i] + p->m02[j][s][i];
		  p->nx02[i][s][j] = p->ex02[i][s][j] - p->im02[i][s][j];
		}
	      else
		{
		  p->im02[i][s][j] = 0.0;
		  p->ex02[i][s][j] = 0.0;
		  p->nx02[i][s][j] = 0.0;
		}
	    }
	  p->im0[i][j] = p->im02[i][0][j] + p->im02[i][1][j];
	  p->ex0[i][j] = p->ex02[i][0][j] + p->ex02[i][1][j];
	  p->nx0[i][j] = p->nx02[i][0][j] + p->nx02[i][1][j];
	}
    }

  double tmp=0.0;
  for(i=0; i<NC; i++)
    {
      tmp = (sum(p->va0[i],NS) - (sum(p->q0[i],NS) + sum(p->ex0[i],NC) - sum(p->im0[i],NC)))/sum(p->va0[i],NS);
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "GDP != C+I+NX for country %d, error = %f\n" RESET,i,tmp);
	  return 1;
	}

      for(s=0; s<NS; s++)
	{
	  tmp = (p->y0[i][s] - SUM_3D_DIM1(p->q02,s,i) - SUM_3D_DIM1(p->m02,s,i))/p->y0[i][s];
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,KRED "supply != demand for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }
	}

      for(s=0; s<NS; s++)
	{
	  tmp = p->y0[i][s] - (p->va0[i][s] + SUM(p->md0[i][s],NS));
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,KRED "go != va + m for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = p->m0[i][s] - SUM(p->m02[i][s],NC);
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,KRED "m != sum(m2) for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }
	}
    }

#ifdef _QUARTERLY
  double uk_nfa_gdp = -17.392652/4.0; // quarterly!!
  double eu_nfa_gdp = -17.809808/4.0; // quartlery!!
#else
  double uk_nfa_gdp = -17.392652;
  double eu_nfa_gdp = -17.809808;
#endif

  p->b0[0] = uk_nfa_gdp;
  p->b0[1] = eu_nfa_gdp * (p->va0[1][0]+p->va0[1][1])/(p->va0[0][0]+p->va0[0][1]);
  p->b0[2] = -sum(p->b0,2);

  return 0;

}

uint calibrate_prod_params(params * p)
{
  uint i,s,r;
  double tmp;
  
  SET_ALL_V(p->lam,NC*NS*NS,0.0);
  SET_ALL_V(p->A,NC*NS,0.0);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  p->A[i][s] = p->va0[i][s] / ( pow(p->k0[i][s],p->alpha[i][s])*pow(p->l0[i][s],1.0-p->alpha[i][s]) );
	  p->lam_va[i][s] = p->va0[i][s] / p->y0[i][s];
	  for(r=0; r<NS; r++)
	    {
	      p->lam[i][s][r] = p->md0[i][s][r] / p->y0[i][s];
	    }
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  tmp = prod_go(p->va0[i][s],p->md0[i][s],p->lam_va[i][s],p->lam[i][s]) - p->y0[i][s];
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "prod_go != y0 for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = prod_va(p->k0[i][s],p->l0[i][s],1.0,p->A[i][s],p->alpha[i][s]) - p->va0[i][s];
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "prod_va != va0 for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = p->y0[i][s] - p->va0[i][s] - sum(p->md0[i][s],NS);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "nonzero profits for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = (1.0-sum(p->lam[i][s],NS)) * (1.0-p->alpha[i][s]) * p->A[i][s]/p->lam_va[i][s] *
	    pow(p->k0[i][s],p->alpha[i][s]) * pow(p->l0[i][s],-(p->alpha[i][s])) - 1.0;
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "labor FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = (1.0-sum(p->lam[i][s],NS)) * (p->alpha[i][s]) * p->A[i][s]/p->lam_va[i][s] *
	    pow(p->k0[i][s],p->alpha[i][s]-1.0) * pow(p->l0[i][s],1.0-p->alpha[i][s]) - 
	    (p->r0[i] + p->delta)/(1.0-p->tauk[i]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "capital FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	}
    }

  return 0;
}

uint calibrate_fin_params(params * p)
{
  uint i,s,j,jj,cnt,idx;
  double tmp;
  double tmp1[NC];

  SET_ALL_V(p->mu,NC*NS*NC,0.0);
  SET_ALL_V(p->M,NC*NS,0.0);
  SET_ALL_V(p->G,NC,0.0);
  SET_ALL_V(p->H,NC*NS,0.0);
  SET_ALL_V(p->eps,NC*NF*NS,0.0);
  SET_ALL_V(p->theta,NC*NS*NC,0.0);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  idx=2;
	  for(j=0; j<NC; j++)
	    {
	      tmp1[j] = pow(p->m02[i][s][j]/p->m02[i][s][idx],1.0 - p->zeta[i][s]);
	    }
	  p->mu[i][s][idx] = 1.0/sum(tmp1,NC);
	  cnt=0;
	  for(j=0; j<NC; j++)
	    {
	      cnt=cnt+1;
	      if(j != idx)
		{
		  if(cnt<NC)
		    {
		      p->mu[i][s][j] = p->mu[i][s][idx]*tmp1[j];
		    }
		  else
		    {
		      p->mu[i][s][j] = 1.0 - sum(p->mu[i][s],NC-1);
		    }
		}
	    
	    }
	  tmp = pow( DOT_PROD_EX(p->m02[i][s],p->mu[i][s],NC,p->zeta[i][s]), 1.0/p->zeta[i][s] );
	  p->M[i][s] = p->m0[i][s]/tmp;
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  idx=2;
	  for(j=0; j<NC; j++)
	    {
	      tmp1[j] = pow(p->q02[i][s][j]/p->q02[i][s][idx],1.0 - p->sig[i][s]);
	    }
	  p->theta[i][s][idx] = 1.0/sum(tmp1,NC);
	  cnt=0;
	  for(j=0; j<NC; j++)
	    {
	      cnt=cnt+1;
	      if(j != idx)
		{
		  if(cnt<NC)
		    {
		      p->theta[i][s][j] = p->theta[i][s][idx]*tmp1[j];
		    }
		  else
		    {
		      p->theta[i][s][j] = 1.0 - sum(p->theta[i][s],NC-1);
		    }
		}
	    
	    }
	  tmp = pow( DOT_PROD_EX(p->q02[i][s],p->theta[i][s],NC,p->sig[i][s]), 1.0/p->sig[i][s] );
	  p->H[i][s] = p->q0[i][s]/tmp;
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  p->eps[i][1][s] = p->i0[i][s] / p->ii0[i];
	}
      p->G[i] = p->ii0[i]/ ( pow(p->i0[i][0],p->eps[i][1][0]) * pow(p->i0[i][1],p->eps[i][1][1]) );
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  tmp = p->m0[i][s] - prod_m(p->m02[i][s],p->M[i][s],p->mu[i][s],p->zeta[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Intermediate Armington production function for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = p->q0[i][s] - prod_q(p->q02[i][s],p->H[i][s],p->theta[i][s],p->sig[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Final Armington production function for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  for(j=0; j<NC; j++)
	    {
	      tmp = 1.0 - p->mu[i][s][j] * pow(p->M[i][s],p->zeta[i][s]) * 
		pow(p->m0[i][s]/p->m02[i][s][j],1.0-p->zeta[i][s]);

	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Intermediate Armington FOC for country/sectorcountry %d/%d/%d, error = %f" RESET, i,s,j,tmp);
		  return 1;
		}

	      tmp = 1.0 - p->theta[i][s][j] * pow(p->H[i][s],p->sig[i][s]) * 
		pow(p->q0[i][s]/p->q02[i][s][j],1.0-p->sig[i][s]);

	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Final Armington FOC for country/sectorcountry %d/%d/%d, error = %f" RESET, i,s,j,tmp);
		  return 1;
		}

	    }

	  if(i==0)
	    {
	      j=1;
	      jj=2;
	    }
	  else if(i==1)
	    {
	      j=0;
	      jj=2;
	    }
	  else
	    {
	      j=0;
	      jj=1;
	    }

	  tmp = 1.0 - (p->mu[i][s][j] / p->mu[i][s][jj]) * 
	    pow(p->m02[i][s][jj]/p->m02[i][s][j],1.0-p->zeta[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Intermediate Armington FOC v2 for country/sector/sector %d/%d/%d, error = %f" RESET,
		      i,s,j,tmp);
	      return 1;
	    }
	      
	  tmp = 1.0 - (p->theta[i][s][j] / p->theta[i][s][jj]) * 
	    pow(p->q02[i][s][jj]/p->q02[i][s][j],1.0-p->sig[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Final Armington FOC v2 for country/sector/sector %d/%d/%d, error = %f" RESET,
		      i,s,j,tmp);
	      return 1;
	    }

	  tmp = p->ii0[i] - prod_inv(p->i0[i],p->eps[i][1],p->G[i]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "ii0 != prod_inv for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }

	  for(s=0; s<NS; s++)
	    {
	      tmp = 1.0 - p->eps[i][1][s]*(p->ii0[i]/p->i0[i][s]);
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Investment FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }
	}
    }

  return 0;
  
}

uint calibrate_hh_params(params * p)
{
  uint i, s, t;
  double tmp;
  double tmp1[NS];
  
  for(i=0; i<NC; i++)
    {
      tmp = pow(p->c0[i][0]/p->c0[i][1],1.0 - p->rho);
      p->eps[i][0][0] = tmp/(1.0+tmp);
      p->eps[i][0][1] = 1.0 - p->eps[i][0][0];

      p->lbar0[i] = 3.0 * p->ll0[i];

      if(fixl==1)
	{
	  p->phi[i]=1.0;
	}
      else
	{
	  tmp = (p->eps[i][0][0] * pow(p->c0[i][0], p->rho) + p->eps[i][0][1]*pow(p->c0[i][1], p->rho)) /
	    (p->lbar0[i] - p->ll0[i]) / p->eps[i][0][0] / pow(p->c0[i][0], p->rho-1.0);
	  p->phi[i] = tmp/(1.0+tmp);
	}

      tmp = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) /
	muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 1) - 1.0;
      if(fabs(tmp)>TINY)
	{
	  fprintf(logfile,KRED "HH intratemp FOC 1 for country %d, error = %f" RESET,i,tmp);
	  return 1;
	}

      if(fixl==0)
	{
	  tmp = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) /
	    mul(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho,p->phi[i],p->psi) - 1.0;
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "HH intratemp FOC 2 for country %d, error = %f" RESET,i,tmp);
	      return 1;
	    }
	}

      for(s=0; s<NS; s++)
	{
	  tmp1[s] = p->c0[i][s] * p->gbgp;
	}
      // note: beta = (1.0+rbgp)/gbgp^(phi*psi-1) > 1
      // but, as long as beta*gbgp^(phi*psi) < 1, we can calculate welfare just fine
      p->beta[i] = muc(p->c0[i],p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) /
	muc(tmp1,p->ll0[i],p->lbar0[i],1.0,1.0,p->eps[i][0],p->rho, p->phi[i], p->psi, 0) / (1.0 + p->rbgp);
    }

  p->wedge_speed = 0.9;

  for(t=0; t<(NT+1); t++)
    {
      for(i=0; i<NC; i++)
	{
	  p->lbar_ts[t][i] = p->lbar0[i] * p->popw_ts[t][i];
	}
    }

  return 0;
}

uint calibrate(params * p)
{
  set_nontargeted_params(p);

  if(load_iomat(p))
    {
      return 1;
    }

  load_ts_params(p);

  if(store_base_period_values(p))
    {
      return 1;
    }

  if(calibrate_prod_params(p))
    {
      return 1;
    }

  if(calibrate_fin_params(p))
    {
      return 1;
    }

  if(calibrate_hh_params(p))
    {
      return 1;
    }

  return 0;
}

uint write_params(const params * p)
{
  FILE * file = fopen("output/params.txt","wb");
  if(file)
    {
      uint i;

      fprintf(file,"--------------------------------------------------\n");
      fprintf(file,"Gross output production\n");
      fprintf(file,"--------------------------------------------------\n");
      
      
      fprintf(file,"\nlam_va\n");
      fprintf(file,"\tUK: %0.4f %0.4f\n",p->lam_va[0][0],p->lam_va[0][1]);
      fprintf(file,"\tEU: %0.4f %0.4f\n",p->lam_va[1][0],p->lam_va[1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f\n",p->lam_va[2][0],p->lam_va[2][1]);
      
      fprintf(file,"\nlam\n");
      fprintf(file,"\tUK: %0.4f %0.4f,\t%0.4f %0.4f\n",
	      p->lam[0][0][0],p->lam[0][0][1],p->lam[0][1][0],p->lam[0][1][1]);
      fprintf(file,"\tEU: %0.4f %0.4f,\t%0.4f %0.4f\n",
	      p->lam[1][0][0],p->lam[1][0][1],p->lam[1][1][0],p->lam[1][1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f,\t%0.4f %0.4f\n",
	      p->lam[2][0][0],p->lam[2][0][1],p->lam[2][1][0],p->lam[2][1][1]);
      
      fprintf(file,"\nalpha\n");
      fprintf(file,"\tUK: %0.4f %0.4f\n",p->alpha[0][0],p->alpha[0][1]);
      fprintf(file,"\tEU: %0.4f %0.4f\n",p->alpha[1][0],p->alpha[1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f\n",p->alpha[2][0],p->alpha[2][1]);

      fprintf(file,"\nA\n");
      fprintf(file,"\tUK: %0.4f %0.4f\n",p->A[0][0],p->A[0][1]);
      fprintf(file,"\tEU: %0.4f %0.4f\n",p->A[1][0],p->A[1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f\n",p->A[2][0],p->A[2][1]);

      fprintf(file,"\n--------------------------------------------------\n");
      fprintf(file,"Armington\n");
      fprintf(file,"--------------------------------------------------\n");

      fprintf(file,"\nzeta\n");
      fprintf(file,"\tUK: %0.4f %0.4f\n",p->zeta[0][0],p->zeta[0][1]);
      fprintf(file,"\tEU: %0.4f %0.4f\n",p->zeta[1][0],p->zeta[1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f\n",p->zeta[2][0],p->zeta[2][1]);

      fprintf(file,"\nmu\n");
      i=0;
      fprintf(file,"\tUK: %0.4f %0.4f %0.4f,\t%0.4f %0.4f %0.4f\n",
	      p->mu[i][0][0],p->mu[i][0][1],p->mu[i][0][2],
	      p->mu[i][1][0],p->mu[i][1][1],p->mu[i][1][2]);
      i=1;
      fprintf(file,"\tEU: %0.4f %0.4f %0.4f,\t%0.4f %0.4f %0.4f\n",
	      p->mu[i][0][0],p->mu[i][0][1],p->mu[i][0][2],
	      p->mu[i][1][0],p->mu[i][1][1],p->mu[i][1][2]);
      i=2;
      fprintf(file,"\tRW: %0.4f %0.4f %0.4f,\t%0.4f %0.4f %0.4f\n",
	      p->mu[i][0][0],p->mu[i][0][1],p->mu[i][0][2],
	      p->mu[i][1][0],p->mu[i][1][1],p->mu[i][1][2]);

      fprintf(file,"\nM\n");
      fprintf(file,"\tUK: %0.4f %0.4f\n",p->M[0][0],p->M[0][1]);
      fprintf(file,"\tEU: %0.4f %0.4f\n",p->M[1][0],p->M[1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f\n",p->M[2][0],p->M[2][1]);

      fprintf(file,"\nsig\n");
      fprintf(file,"\tUK: %0.4f %0.4f\n",p->sig[0][0],p->sig[0][1]);
      fprintf(file,"\tEU: %0.4f %0.4f\n",p->sig[1][0],p->sig[1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f\n",p->sig[2][0],p->sig[2][1]);

      fprintf(file,"\ntheta\n");
      i=0;
      fprintf(file,"\tUK: %0.4f %0.4f %0.4f,\t%0.4f %0.4f %0.4f\n",
	      p->theta[i][0][0],p->theta[i][0][1],p->theta[i][0][2],
	      p->theta[i][1][0],p->theta[i][1][1],p->theta[i][1][2]);
      i=1;
      fprintf(file,"\tEU: %0.4f %0.4f %0.4f,\t%0.4f %0.4f %0.4f\n",
	      p->theta[i][0][0],p->theta[i][0][1],p->theta[i][0][2],
	      p->theta[i][1][0],p->theta[i][1][1],p->theta[i][1][2]);
      i=2;
      fprintf(file,"\tRW: %0.4f %0.4f %0.4f,\t%0.4f %0.4f %0.4f\n",
	      p->theta[i][0][0],p->theta[i][0][1],p->theta[i][0][2],
	      p->theta[i][1][0],p->theta[i][1][1],p->theta[i][1][2]);

      fprintf(file,"\nH\n");
      fprintf(file,"\tUK: %0.4f %0.4f\n",p->H[0][0],p->H[0][1]);
      fprintf(file,"\tEU: %0.4f %0.4f\n",p->H[1][0],p->H[1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f\n",p->H[2][0],p->H[2][1]);

      fprintf(file,"\n--------------------------------------------------\n");
      fprintf(file,"Final demand\n");
      fprintf(file,"--------------------------------------------------\n");

      fprintf(file,"\nrho\t%0.4f\n",p->rho);

      fprintf(file,"\neps\n");
      fprintf(file,"\tUK: %0.4f %0.4f,\t%0.4f %0.4f\n",
	      p->eps[0][0][0],p->eps[0][0][1],p->eps[0][1][0],p->eps[0][1][1]);
      fprintf(file,"\tEU: %0.4f %0.4f,\t%0.4f %0.4f\n",
	      p->eps[1][0][0],p->eps[1][0][1],p->eps[1][1][0],p->eps[1][1][1]);
      fprintf(file,"\tRW: %0.4f %0.4f,\t%0.4f %0.4f\n",
	      p->eps[2][0][0],p->eps[2][0][1],p->eps[2][1][0],p->eps[2][1][1]);

      fprintf(file,"\nG\n");
      fprintf(file,"\tUK: %0.4f\n",p->G[0]);
      fprintf(file,"\tEU: %0.4f\n",p->G[1]);
      fprintf(file,"\tRW: %0.4f\n",p->G[2]);


      fprintf(file,"\n--------------------------------------------------\n");
      fprintf(file,"Households\n");
      fprintf(file,"--------------------------------------------------\n");

      fprintf(file,"\ndelta\t%0.4f\n",p->delta);
      fprintf(file,"rbgp\t%0.4f\n",p->rbgp);
      fprintf(file,"etaK\t%0.4f\n",p->etaK);
      fprintf(file,"psi\t%0.4f\n",p->psi);

      fprintf(file,"\ntauk\n");
      fprintf(file,"\tUK: %0.4f\n",p->tauk[0]);
      fprintf(file,"\tEU: %0.4f\n",p->tauk[1]);
      fprintf(file,"\tRW: %0.4f\n",p->tauk[2]);

      fprintf(file,"\nbeta\n");
      fprintf(file,"\tUK: %0.4f\n",p->beta[0]);
      fprintf(file,"\tEU: %0.4f\n",p->beta[1]);
      fprintf(file,"\tRW: %0.4f\n",p->beta[2]);

      fprintf(file,"\nphi\n");
      fprintf(file,"\tUK: %0.4f\n",p->phi[0]);
      fprintf(file,"\tEU: %0.4f\n",p->phi[1]);
      fprintf(file,"\tRW: %0.4f\n",p->phi[2]);

      fprintf(file,"\nr0\n");
      fprintf(file,"\tUK: %0.4f\n",p->r0[0]);
      fprintf(file,"\tEU: %0.4f\n",p->r0[1]);
      fprintf(file,"\tRW: %0.4f\n",p->r0[2]);

      fprintf(file,"\ntauk\n");
      fprintf(file,"\tUK: %0.4f\n",p->tauk[0]);
      fprintf(file,"\tEU: %0.4f\n",p->tauk[1]);
      fprintf(file,"\tRW: %0.4f\n",p->tauk[2]);

      fprintf(file,"\nk0\n");
      fprintf(file,"\tUK: %0.4f\n",p->kk0[0]);
      fprintf(file,"\tEU: %0.4f\n",p->kk0[1]);
      fprintf(file,"\tRW: %0.4f\n",p->kk0[2]);

      fprintf(file,"\nb0\n");
      fprintf(file,"\tUK: %0.4f\n",p->b0[0]);
      fprintf(file,"\tEU: %0.4f\n",p->b0[1]);
      fprintf(file,"\tRW: %0.4f\n",p->b0[2]);

      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write parameters!\n" RESET);
      return 1;
    }
}

#endif
