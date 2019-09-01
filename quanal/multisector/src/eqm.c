#ifndef __EQM_C__
#define __EQM_C__

#include "eqm.h"

// wage (NC), sectoral capital (NC*NS), sectoral labor (NC*NS), gross output prices (NC*NS),
// sectoral consumption (NC*NS)
const uint nbgp = NC+4*NC*NS;

void set_neqm()
{
  // key eqm vars used in solver: wage (NC), sectoral investment (NC*NS, but not in last period),
  // sectoral labor (NC*NS), gross output prices (NC*NS), sectoral consumption (NC*NS), 
  // bonds (NC-1, but not in first period), bond price (1, but not in last period)

  // deterministic no-Brexit counterfactual
  if(scenario==0)
    {
      neqm = (NT+1)*NC + 3*(NT+1)*NC*NS + NT*NC*NS + NT*(NC-1) + NT;
      if(f_adj_cost)
	{
	  neqm = neqm + NT*NC*NS*NC + NT*NC*NS;
	}
      if(m_adj_cost)
	{
	  neqm = neqm + NT*NC*NS*NC + NT*NC*NS;
	}
    }
  // stochastic equilibrium
  else if(scenario == 1)
    {
      // we have TVOTE-TREF periods with only one branch: the period in which the referendum is announced
      uint nn = TVOTE-TREF;
      neqm = nn*(NC + 4*NC*NS + (NC-1) + 1);
      if(f_adj_cost)
	{
	  neqm = neqm + nn*NC*NS*NC + nn*NC*NS;
	}
      if(m_adj_cost)
	{
	  neqm = neqm + nn*NC*NS*NC + nn*NC*NS;
	}

      // in the period of the vote, the tree splits...

      // ...along branch 0, the vote fails and Brexit never occurs... this is a deterministic branch that
      // goes all the way to the end
      nn = (NT+1-TVOTE);
      neqm = neqm + (nn)*NC + 3*(nn)*NC*NS + (nn-1)*NC*NS + (nn-1)*(NC-1) + (nn-1);
      if(f_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
      if(m_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
      
      // ... along branch 1, the vote happens and we have to wait until the Brexit period
      // to find out what happens... this section is TBREXIT-TVOTE periods long
      nn = TBREXIT-TVOTE;
      neqm = neqm + nn*(NC + 4*NC*NS + (NC-1) + 1);
      if(f_adj_cost)
	{
	  neqm = neqm + nn*NC*NS*NC + nn*NC*NS;
	}
      if(m_adj_cost)
	{
	  neqm = neqm + nn*NC*NS*NC + nn*NC*NS;
	}
      
      // after Brexit, we have 2 possible histories that could arise
      nn = (NT+1-TBREXIT);
      neqm = neqm + (nn)*NC + 3*(nn)*NC*NS + (nn-1)*NC*NS + (nn-1)*(NC-1) + (nn-1);
      neqm = neqm + (nn)*NC + 3*(nn)*NC*NS + (nn-1)*NC*NS + (nn-1)*(NC-1) + (nn-1);
      if(f_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
      if(m_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
    }
  // perfect foresight equilibrium with Brexit... from the period of the vote onwards
  else if(scenario == 4 || scenario == 5)
    {
      uint nn = (NT+1-TVOTE);
      neqm = (nn)*NC + 3*(nn)*NC*NS + (nn-1)*NC*NS + (nn-1)*(NC-1) + (nn-1);
      if(f_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
      if(m_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
    }
  // perfect foresight equilibrium with Brexit... from the period of the referendum onwards
  else if(scenario == 2 || scenario == 3)
    {
      uint nn = (NT+1-TREF);
      neqm = (nn)*NC + 3*(nn)*NC*NS + (nn-1)*NC*NS + (nn-1)*(NC-1) + (nn-1);
      if(f_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
      if(m_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*NS*NC + (nn-1)*NC*NS;
	}
    }

}

void init_vars(eqm * e)
{
  SET_ALL_V(e->b_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->cc_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ii_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ll_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->kk_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->cpi_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->pi_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->w_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->rk_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ngdp_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->rgdp_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->iy_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->lp_agg_t,(NT+1)*NC,0.0);

  SET_ALL_V(e->ex_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->im_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->nx_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->exf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->imf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->nxf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->exm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->imm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->nxm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rer_t,(NT+1)*NC*NC,0.0);

  SET_ALL_V(e->rex_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rim_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rexf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rimf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rexm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rimm_t,(NT+1)*NC*NC,0.0);

  SET_ALL_V(e->y_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->va_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->md_t,(NT+1)*NC*NS*NS,0.0);
  SET_ALL_V(e->k_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->l_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->is_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->lp_t,(NT+1)*NC*NS,0.0);

  SET_ALL_V(e->exs_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->ims_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->nxs_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->rexs_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->rims_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->pm_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->m_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->m2_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->p_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->q_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->q2_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->c_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->i_t,(NT+1)*NC*NS,0.0);
}

void copy_vars(eqm * e1, const eqm * e0)
{
  memcpy((double *)(e1->pb_t),(const double *)(e0->pb_t),sizeof(double)*(NT+1));

  memcpy((double *)(e1->b_t),(const double *)(e0->b_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->cc_t),(const double *)(e0->cc_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ii_t),(const double *)(e0->ii_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ll_t),(const double *)(e0->ll_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->kk_t),(const double *)(e0->kk_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->cpi_t),
	 (const double *)(e0->cpi_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->pi_t),(const double *)(e0->pi_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->w_t),(const double *)(e0->w_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->rk_t),(const double *)(e0->rk_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ngdp_t),
	 (const double *)(e0->ngdp_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->rgdp_t),
	 (const double *)(e0->rgdp_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->iy_t),(const double *)(e0->iy_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->lp_agg_t),
	 (const double *)(e0->lp_agg_t),
	 sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->ex_t),
	 (const double *)(e0->ex_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->im_t),
	 (const double *)(e0->im_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->nx_t),
	 (const double *)(e0->nx_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->exf_t),
	 (const double *)(e0->exf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->imf_t),
	 (const double *)(e0->imf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->nxf_t),
	 (const double *)(e0->nxf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->exm_t),
	 (const double *)(e0->exm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->imm_t),
	 (const double *)(e0->imm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->nxm_t),
	 (const double *)(e0->nxm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rer_t),
	 (const double *)(e0->rer_t),
	 sizeof(double)*(NT+1)*NC*NC);

  memcpy((double *)(e1->rex_t),
	 (const double *)(e0->rex_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rim_t),
	 (const double *)(e0->rim_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rexf_t),
	 (const double *)(e0->rexf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rimf_t),
	 (const double *)(e0->rimf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rexm_t),
	 (const double *)(e0->rexm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rimm_t),
	 (const double *)(e0->rimm_t),
	 sizeof(double)*(NT+1)*NC*NC);


  memcpy((double *)(e1->y_t),
	 (const double *)(e0->y_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->py_t),
	 (const double *)(e0->py_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->va_t),
	 (const double *)(e0->va_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->k_t),
	 (const double *)(e0->k_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->l_t),
	 (const double *)(e0->l_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->is_t),
	 (const double *)(e0->is_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->md_t),
	 (const double *)(e0->md_t),
	 sizeof(double)*(NT+1)*NC*NS*NS);
  memcpy((double *)(e1->lp_t),
	 (const double *)(e0->lp_t),
	 sizeof(double)*(NT+1)*NC*NS);

  memcpy((double *)(e1->exs_t),
	 (const double *)(e0->exs_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->ims_t),
	 (const double *)(e0->ims_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->nxs_t),
	 (const double *)(e0->nxs_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->rexs_t),
	 (const double *)(e0->rexs_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->rims_t),
	 (const double *)(e0->rims_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->pm_t),
	 (const double *)(e0->pm_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->m_t),
	 (const double *)(e0->m_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->m2_t),
	 (const double *)(e0->m2_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->p_t),
	 (const double *)(e0->p_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->q_t),
	 (const double *)(e0->q_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->q2_t),
	 (const double *)(e0->q2_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->c_t),
	 (const double *)(e0->c_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->i_t),
	 (const double *)(e0->i_t),
	 sizeof(double)*(NT+1)*NC*NS);

  memcpy((double *)(e1->welfare_t),
	 (const double *)(e0->welfare_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare2_t),
	 (const double *)(e0->welfare2_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare_cost_t),
	 (const double *)(e0->welfare_cost_t),
	 sizeof(double)*(NT+1)*NC);
}

uint stack_bgp_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t = NT;
  
  COPY_SUBVECTOR_LOG(myx+nx,e->w_t[t],NC);
  nx=nx+NC;

  COPY_SUBVECTOR_LOG(myx+nx,e->k_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->l_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->py_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->c_t[t],NC*NS);
  nx=nx+NC*NS;

  if(nx != nbgp)
    {
      fprintf(logfile,KRED "Error stacking bgp vars! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }

    return 0;
}

uint unstack_bgp_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t = NT;

  copy_subvector_exp( (double *)(e->w_t[t]), myx+nx, NC);
  nx=nx+NC;

  COPY_SUBVECTOR_EXP(e->k_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->l_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->py_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->c_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  if(nx != nbgp)
    {
      fprintf(logfile,KRED "Error stacking bgp vars! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }

    return 0;
}

uint stack_eqm_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;
  if(scenario == 2 || scenario == 3)
    {
      nn = NT+1-TREF;
      t0 = TREF;
    }
  else if(scenario == 4 || scenario == 5)
    {
      nn = NT+1-TVOTE;
      t0 = TVOTE;
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->is_t[t0]),(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->l_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->py_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->c_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->b_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->pb_t[t0]),(nn-1));
  nx = nx + (nn-1);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error stacking eqm vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

    return 0;
}

uint unstack_eqm_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;

  if(scenario == 2 || scenario == 3)
    {
      nn = NT+1-TREF;
      t0 = TREF;
    }
  else if(scenario == 4 || scenario == 5)
    {
      nn = NT+1-TVOTE;
      t0 = TVOTE;
    }

  COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;
  
  COPY_SUBVECTOR_EXP(&(e->is_t[t0]),myx+nx,(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->l_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->py_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->c_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->b_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->pb_t[t0]),myx+nx,(nn-1));
  nx = nx + (nn-1);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;

      COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;

      COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking eqm vars! nx = %d, neqm = %d\n" RESET ,nx,neqm);
      return 1;
    }

  return 0;
}

uint stack_stoch_eqm_vars(double * myx, const stoch_eqm * se)
{
  uint nx = 0;

  // first stack the deterministic part for the periods between the referendum announcement and the vote
  uint t0 = TREF;
  uint nn = TVOTE-TREF;
  const eqm * e = &((se->eee)[0]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->is_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->l_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->py_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->c_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<=TVOTE; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->b_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->pb_t[t0]),(nn));
  nx = nx + (nn);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;
    }

  // now stack the deterministic part associated with the "no vote" branch
  t0 = TVOTE;
  nn = NT+1-TVOTE;  
  e = &((se->eee)[0]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->is_t[t0]),(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->l_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->py_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->c_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->b_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->pb_t[t0]),(nn-1));
  nx = nx + (nn-1);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;
      
      COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;
      
      COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  // now stack the TBREXIT-TVOTE periods associated with the "yes vote" branch
  t0 = TVOTE;
  nn = TBREXIT-TVOTE;
  e = &((se->eee)[1]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->is_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->l_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->py_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->c_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  for(t=t0+1; t<=TBREXIT; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->b_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->pb_t[t0]),(nn));
  nx = nx + (nn);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;
    }

  // now stack the post-Brexit part for the hard and soft branches
  t0 = TBREXIT;
  nn = NT+1-TBREXIT;

  uint ih;
  for(ih=1; ih<NHIST; ih++)
    {
      e = &((se->eee)[ih]);

      COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->is_t[t0]),(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->l_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->py_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->c_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      *(myx+nx) = e->b_t[t][i];
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_LOG(myx+nx,&(e->pb_t[t0]),(nn-1));
      nx = nx + (nn-1);

      if(f_adj_cost)
	{
	  COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn-1)*NC*NS*NC);
	  nx = nx + (nn-1)*NC*NS*NC;

	  COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn-1)*NC*NS);
	  nx = nx + (nn-1)*NC*NS;
	}

      if(m_adj_cost)
	{
	  COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn-1)*NC*NS*NC);
	  nx = nx + (nn-1)*NC*NS*NC;

	  COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn-1)*NC*NS);
	  nx = nx + (nn-1)*NC*NS;
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error stacking stoch_eqm vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

    return 0;
}

uint unstack_stoch_eqm_vars(stoch_eqm * se, const double * myx)
{
  uint nx = 0;

  // first unstack the deterministic part for the periods between the referendum annoumcement
  // and the vot
  uint i = 0;
  uint t = 0;
  uint t0 = TREF;
  uint nn = TVOTE-TREF;
  eqm * e = &((se->eee)[0]);

  // when we unstack we want to copy over to all histories for the deterministic part,
  // since they have to be identical for this period
  int ih;
  for(ih=NHIST-1; ih>=0; ih--)
    {
      e = &((se->eee)[ih]);
      nx = 0;

      COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;
  
      COPY_SUBVECTOR_EXP(&(e->is_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->l_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->py_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->c_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      for(t=t0+1; t<=TVOTE; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->b_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->pb_t[t0]),myx+nx,(nn));
      nx = nx + (nn);

      if(f_adj_cost)
	{
	  COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn)*NC*NS*NC);
	  nx = nx + nn*NC*NS*NC;

	  COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn)*NC*NS);
	  nx = nx + nn*NC*NS;
	}

      if(m_adj_cost)
	{
	  COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn)*NC*NS*NC);
	  nx = nx + nn*NC*NS*NC;

	  COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn)*NC*NS);
	  nx = nx + nn*NC*NS;
	}
    }

  // now do the deterministic part for the "stay" branch
  nn = NT+1-TVOTE;
  t0 = TVOTE;

  e = &((se->eee)[0]);
  COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;
  
  COPY_SUBVECTOR_EXP(&(e->is_t[t0]),myx+nx,(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->l_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->py_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->c_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->b_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->pb_t[t0]),myx+nx,(nn-1));
  nx = nx + (nn-1);

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;

      COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn-1)*NC*NS*NC);
      nx = nx + (nn-1)*NC*NS*NC;

      COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  // now do the TBREXIT-TVOTE period... here we need to copy to both branches 1 and 2 since they
  // must be identical here
  nn = TBREXIT-TVOTE;
  t0 = TVOTE;
  uint nx2=nx;
  for(ih=(NHIST-1); ih>=1; ih--)
    {
      e = &((se->eee)[ih]);
      nx = nx2;

      COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;
  
      COPY_SUBVECTOR_EXP(&(e->is_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->l_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->py_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->c_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      for(t=t0+1; t<=TBREXIT; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->b_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->pb_t[t0]),myx+nx,(nn));
      nx = nx + (nn);

      if(f_adj_cost)
	{
	  COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn)*NC*NS*NC);
	  nx = nx + nn*NC*NS*NC;

	  COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn)*NC*NS);
	  nx = nx + nn*NC*NS;
	}

      if(m_adj_cost)
	{
	  COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn)*NC*NS*NC);
	  nx = nx + nn*NC*NS*NC;

	  COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn)*NC*NS);
	  nx = nx + nn*NC*NS;
	}
    }

  // last, do the 2 post-Brexit branches
  nn = NT+1-TBREXIT;
  t0 = TBREXIT;
  
  for(ih=1; ih<NHIST; ih++)
    {
      e = &((se->eee)[ih]);

      COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;
  
      COPY_SUBVECTOR_EXP(&(e->is_t[t0]),myx+nx,(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->l_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->py_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->c_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->b_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->pb_t[t0]),myx+nx,(nn-1));
      nx = nx + (nn-1);

      if(f_adj_cost)
	{
	  COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn-1)*NC*NS*NC);
	  nx = nx + (nn-1)*NC*NS*NC;

	  COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn-1)*NC*NS);
	  nx = nx + (nn-1)*NC*NS;
	}

      if(m_adj_cost)
	{
	  COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn-1)*NC*NS*NC);
	  nx = nx + (nn-1)*NC*NS*NC;

	  COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn-1)*NC*NS);
	  nx = nx + (nn-1)*NC*NS;
	}

    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking stoch_eqm vars! nx = %d, neqm = %d\n" RESET ,nx,neqm);
      return 1;
    }

  return 0;
}

uint set_initial_bgp_guess()
{
  uint i, t, s;
  eqm * e = &(eee0[0]);
  params * p = &(ppp0[0]);
  t=NT;
  
  for(i=0; i<NC; i++)
    {
      e->w_t[t][i] = 1.0;
      for(s=0; s<NS; s++)
	{
	  e->c_t[t][i][s] = p->c0[i][s] * p->gam_ts[t][i][s] * (p->lbar_ts[t][i]/p->lbar0[i]);
	  e->l_t[t][i][s] = p->l0[i][s] * (p->lbar_ts[t][i]/p->lbar0[i]);
	  e->k_t[t][i][s] = p->k0[i][s] * p->gam_ts[t][i][s] * (p->lbar_ts[t][i]/p->lbar0[i]);
	  e->py_t[t][i][s] = 1.0;
	}
    }

  if(stack_bgp_vars(solver_x->data,e))
    {
      fprintf(logfile,KRED "Failed to create guess for balanced growth path!\n" RESET);
      return 1;
    }
  else
    {
      return 0;
    }
}

uint set_initial_eqm_guess()
{
  uint i,s,t,j;
  double bb[NC];

  eqm * e = &(eee0[0]);
  params * p = &(ppp0[0]);

  bb[0] = p->b0[0] * pow(p->gbgp,NT);
  bb[1] = p->b0[1] * pow(p->gbgp,NT);
  bb[2] = p->b0[2] * pow(p->gbgp,NT);

  free_solver_mem();
  solver_n = nbgp;
  alloc_solver_mem();
  if(solve_bgp(bb))
    {
      fprintf(logfile, KRED "Error solving for balanced growth path!\n");
      return 1;
    }
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

  // first construct bond guess... a little awkward to logspace this because we have to deal with
  // absolute values
  double tmpb0[NC] = {fabs(p->b0[0]),fabs(p->b0[1]),fabs(p->b0[2])};
  double tmpb1[NC] = {fabs(bb[0]),fabs(bb[1]),fabs(bb[2])};
  LOGSPACE_2D(tmpb0,tmpb1,NT+1,NC,e->b_t);
  for(i=0; i<(NC-1); i++)
    {
      if(fabs(p->b0[i])<1.0e-6)
	{
	  for(t=0; t<(NT+1); t++)
	    {
	      e->b_t[t][i] = 0.0;
	    }
	}
      else
	{
	  if(p->b0[i] < -TINY)
	    {
	      for(t=0; t<(NT+1); t++)
		{
		  e->b_t[t][i] = -e->b_t[t][i];
		}
	    }
	}
    }
  set_all_v(e->pb_t,NT+1,e->pb_t[NT]);

  // now construct guesses for prices real variables
  double tmpp[NC] = {1.0,1.0,1.0};
  double tmpp2[NC][NS] = {{1.0,1.0},{1.0,1.0},{1.0,1.0}};
  LOGSPACE_2D(p->k0,e->k_t[NT],NT+1,NC*NS,e->k_t);
  LOGSPACE_2D(p->l0,e->l_t[NT],NT+1,NC*NS,e->l_t);
  LOGSPACE_2D(p->c0,e->c_t[NT],NT+1,NC*NS,e->c_t);
  LOGSPACE_2D(p->q0,e->q_t[NT],NT+1,NC*NS,e->q_t);
  LOGSPACE_2D(p->m0,e->m_t[NT],NT+1,NC*NS,e->m_t);
  LINSPACE_2D(tmpp,e->w_t[NT],NT+1,NC,e->w_t);
  LINSPACE_2D(tmpp2,e->py_t[NT],NT+1,NC*NS,e->py_t);
  LINSPACE_2D(tmpp2,e->p_t[NT],NT+1,NC*NS,e->p_t);
  LINSPACE_2D(tmpp2,e->pm_t[NT],NT+1,NC*NS,e->pm_t);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(t=0; t<NT; t++)
	    {
	      e->is_t[t][i][s] = e->k_t[t+1][i][s] - (1.0-p->delta)*e->k_t[t][i][s];
	    }
	}
    }

  for(t=0; t<(NT+1); t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  e->q2_t[t][i][s][j] = (p->q02[i][s][j]/p->q0[i][s]) * e->q_t[t][i][s];
		  e->m2_t[t][i][s][j] = (p->m02[i][s][j]/p->m0[i][s]) * e->m_t[t][i][s];
		}
	    }
      
	}
    }

  if(stack_eqm_vars(solver_x->data,e))
    {
      fprintf(logfile,KRED "Failed to create guess for balanced growth path!\n" RESET);
      return 1;
    }
  else
    {      
      return 0;
    }
}

uint write_bgp_vars(const eqm * e, const params * p)
{
  uint t = NT;

#ifdef _QUARTERLY
  FILE * file = fopen("output/bgp.txt","wb");
#else
  FILE * file = fopen("output/bgp_quarterly.txt","wb");
#endif

  if(file)
    {

      fprintf(file,"private bonds (exo. state): %0.12f\t%0.12f\t%0.12f\n",e->b_t[t][0],e->b_t[t][1],e->b_t[t][2]);
      fprintf(file,"bond price: %0.12f\n",e->pb_t[t]);
      fprintf(file,"wages: %0.12f\t%0.12f\t%0.12f\n",e->w_t[t][0],e->w_t[t][1],e->w_t[t][2]);

      uint i=0;
      fprintf(file,"UK capital: %0.12f\t%0.12f\n",e->k_t[t][i][0],e->k_t[t][i][1]);      
      i=1;
      fprintf(file,"EU capital: %0.12f\t%0.12f\n",e->k_t[t][i][0],e->k_t[t][i][1]);      
      i=2;
      fprintf(file,"RW capital: %0.12f\t%0.12f\n",e->k_t[t][i][0],e->k_t[t][i][1]);      

      i=0;
      fprintf(file,"UK labor: %0.12f\t%0.12f\n",e->l_t[t][i][0],e->l_t[t][i][1]);
      i=1;
      fprintf(file,"EU labor: %0.12f\t%0.12f\n",e->l_t[t][i][0],e->l_t[t][i][1]);
      i=2;
      fprintf(file,"RW labor: %0.12f\t%0.12f\n",e->l_t[t][i][0],e->l_t[t][i][1]);

      i=0;
      fprintf(file,"UK go: %0.12f\t%0.12f\n",e->y_t[t][i][0],e->y_t[t][i][1]);
      i=1;
      fprintf(file,"EU go: %0.12f\t%0.12f\n",e->y_t[t][i][0],e->y_t[t][i][1]);
      i=2;
      fprintf(file,"RW go: %0.12f\t%0.12f\n",e->y_t[t][i][0],e->y_t[t][i][1]);

      i=0;
      fprintf(file,"UK go prices: %0.12f\t%0.12f\n",e->py_t[t][i][0],e->py_t[t][i][1]);
      i=1;
      fprintf(file,"EU go prices: %0.12f\t%0.12f\n",e->py_t[t][i][0],e->py_t[t][i][1]);
      i=2;
      fprintf(file,"RW go prices: %0.12f\t%0.12f\n",e->py_t[t][i][0],e->py_t[t][i][1]);
      
      i=0;
      fprintf(file,"UK NX: %0.12f\t%0.12f\t%0.12f\n",e->nx_t[t][i][0],e->nx_t[t][i][1],e->nx_t[t][i][2]);
      i=1;
      fprintf(file,"EU NX: %0.12f\t%0.12f\t%0.12f\n",e->nx_t[t][i][0],e->nx_t[t][i][1],e->nx_t[t][i][2]);
      i=2;
      fprintf(file,"RW NX: %0.12f\t%0.12f\t%0.12f\n",e->nx_t[t][i][0],e->nx_t[t][i][1],e->nx_t[t][i][2]);

      i=0;
      fprintf(file,"UK GDP: %0.12f\n",e->rgdp_t[t][i]);
      i=1;
      fprintf(file,"EU GDP: %0.12f\n",e->rgdp_t[t][i]);
      i=2;
      fprintf(file,"RW GDP: %0.12f\n",e->rgdp_t[t][i]);

      i=0;
      fprintf(file,"UK K/Y: %0.12f\n",100.0*e->kk_t[t][i]/e->rgdp_t[t][i]);
      i=1;
      fprintf(file,"EU K/Y: %0.12f\n",100.0*e->kk_t[t][i]/e->rgdp_t[t][i]);
      i=2;
      fprintf(file,"RW K/Y: %0.12f\n",100.0*e->kk_t[t][i]/e->rgdp_t[t][i]);

      fclose(file);

      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write balanced growth path!\n" RESET);
      return 1;
    }
}

uint write_eqm_vars(const eqm * e, char * fname, uint i)
{
  char * fname2 = concat("output/",fname);

  char * fname3;
#ifdef _QUARTERLY
  if(sensitivity==0 && m_adj_cost==0 && f_adj_cost==0 && fixl<2 && tfp_flag==0)
    {
      fname3 = concat(fname2,"_quarterly.csv");
    }
  else if(m_adj_cost==1 || f_adj_cost==1)
    {
      fname3 = concat(fname2,"_frictions_quarterly.csv");
    }
  else if(sensitivity==1)
    {
      fname3 = concat(fname2,"_pi0.25_quarterly.csv");
    }
  else if(sensitivity==2)
    {
      fname3 = concat(fname2,"_pi0.75_quarterly.csv");
    }
  else if(sensitivity==3)
    {
      fname3 = concat(fname2,"_rra5_quarterly.csv");
    }
  else if(sensitivity==4)
    {
      fname3 = concat(fname2,"_cobb_quarterly.csv");
    }
  else if(fixl==2)
    {
      fname3 = concat(fname2,"_stickyw_quarterly.csv");
    }
  else if(tfp_flag==1)
    {
      fname3 = concat(fname2,"_tfp_quarterly.csv");
    }
#else
  if(sensitivity==0 && m_adj_cost==0 && f_adj_cost==0 && fixl<2 && tfp_flag==0)
    {
      fname3 = concat(fname2,".csv");
    }
  else if(m_adj_cost==1 || f_adj_cost==1)
    {
      fname3 = concat(fname2,"_frictions.csv");
    }
  else if(sensitivity==1)
    {
      fname3 = concat(fname2,"_pi0.25.csv");
    }
  else if(sensitivity==2)
    {
      fname3 = concat(fname2,"_pi0.75.csv");
    }
  else if(sensitivity==3)
    {
      fname3 = concat(fname2,"_rra5.csv");
    }
  else if(sensitivity==4)
    {
      fname3 = concat(fname2,"_cobb.csv");
    }
  else if(fixl==2)
    {
      fname3 = concat(fname2,"_stickyw.csv");
    }
  else if(tfp_flag==1)
    {
      fname3 = concat(fname2,"_tfp.csv");
    }
#endif

  FILE * file = fopen(fname3,"wb");
  free(fname2);

  if(file)
    {
      uint s,j,t;
      fprintf(file,"period,rgdp,ngdp,iy,c,W,sW,cW,ceW1,ceW2");
      for(s=0;s<NS;s++)
	{
	  fprintf(file,",y%d,va%d,i%d,c%d,l%d",s,s,s,s,s);
	}
      for(j=0; j<NC; j++)
	{
	  if(j!=i)
	    {
	      fprintf(file,",rer%d,nx%d,rex%d,rim%d,rexf%d,rimf%d",
		      j,j,j,j,j,j);
	      for(s=0;s<NS;s++)
		{
		  fprintf(file,",rexs%d-%d,rims%d-%d,rexsf%d-%d,rimsf%d-%d",s,j,s,j,s,j,s,j);
		}
	    }
	}
      fprintf(file,"\n");

      for(t=0;t<(NT+1);t++)
	{
	  fprintf(file,"%d,",t);
	  fprintf(file,"%0.16f,",e->rgdp_t[t][i]);
	  fprintf(file,"%0.16f,",e->ngdp_t[t][i]);
	  fprintf(file,"%0.16f,",e->iy_t[t][i]);
	  fprintf(file,"%0.16f,",e->cc_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare2_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare_cost_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare3_t[t][i]);
	  fprintf(file,"%0.16f",e->welfare4_t[t][i]);
	  for(s=0; s<NS; s++)
	    {
	      fprintf(file,",%0.16f,",e->y_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->va_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->is_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->c_t[t][i][s]);
	      fprintf(file,"%0.16f",e->l_t[t][i][s]);
	    }
	  for(j=0; j<NC; j++)
	    {
	      if(j!=i)
		{
		  fprintf(file,",%0.16f,",e->rer_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->nx_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->rex_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->rim_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->rexf_t[t][i][j]);
		  fprintf(file,"%0.16f",e->rimf_t[t][i][j]);
		  for(s=0; s<NS; s++)
		    {
		      fprintf(file,",%0.16f,",e->rexs_t[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->rims_t[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->q2_t[t][i][s][j]);
		      fprintf(file,"%0.16f",e->q2_t[t][j][s][i]);
		    }
		}
	    }
	  fprintf(file,"\n");
	}
      
      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write equilibrium vars!\n" RESET);
      return 1;
    }
}

void set_vars(eqm * e, const params * p, uint t, uint bgp)
{
  uint i,s,r,j;

  e->b_t[t][2] = -(e->b_t[t][0]+e->b_t[t][1]);

  SET_ALL_V(e->ngdp_t[t],NC,0.0);
  SET_ALL_V(e->rgdp_t[t],NC,0.0);
  SET_ALL_V(e->ex_t[t],NC*NC,0.0);
  SET_ALL_V(e->im_t[t],NC*NC,0.0);
  SET_ALL_V(e->nx_t[t],NC*NC,0.0);
  SET_ALL_V(e->exf_t[t],NC*NC,0.0);
  SET_ALL_V(e->imf_t[t],NC*NC,0.0);
  SET_ALL_V(e->nxf_t[t],NC*NC,0.0);
  SET_ALL_V(e->exm_t[t],NC*NC,0.0);
  SET_ALL_V(e->imm_t[t],NC*NC,0.0);
  SET_ALL_V(e->nxm_t[t],NC*NC,0.0);
  SET_ALL_V(e->exs_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->ims_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->nxs_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->rex_t[t],NC*NC,0.0);
  SET_ALL_V(e->rim_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexf_t[t],NC*NC,0.0);
  SET_ALL_V(e->rimf_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexm_t[t],NC*NC,0.0);
  SET_ALL_V(e->rimm_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexs_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->rims_t[t],NC*NS*NC,0.0);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
  
	  e->va_t[t][i][s] = (p->a_ts[t][i][s]) * prod_va(e->k_t[t][i][s],e->l_t[t][i][s],p->gam_ts[t][i][s],p->A[i][s],p->alpha[i][s]);
	  e->y_t[t][i][s] = e->va_t[t][i][s]/p->lam_va[i][s];
	  e->ngdp_t[t][i] = e->ngdp_t[t][i] + e->py_t[t][i][s]*e->y_t[t][i][s];
	  e->rgdp_t[t][i] = e->rgdp_t[t][i] + e->y_t[t][i][s];
	  e->lp_t[t][i][s] = e->va_t[t][i][s]/e->l_t[t][i][s];

	  for(r=0; r<NS; r++)
	    {
	      e->md_t[t][i][s][r] = e->y_t[t][i][s]*p->lam[i][s][r];	      
	    }
	}

      e->ll_t[t][i] = e->l_t[t][i][0] + e->l_t[t][i][1];

      // adjustment costs are in units of labor, so we need to incremement labor demand by them
      // ??? do we need to do something to GDP as well with this?
      if(f_adj_cost && t<NT)
	{
	  for(s=0; s<NS; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  if(j!=i)
		    {
		      double tmp = 0.0;
		      if(t>0)
			{
			  tmp = e->q2_t[t][i][s][j]/e->q2_t[t-1][i][s][j] - 1.0;
			}
		      else
			{
			  tmp = e->q2_t[t][i][s][j]/p->q02[i][s][j] - 1.0;
			}
		      e->ll_t[t][i] = e->ll_t[t][i] + 
			(p->etaF/2.0) * tmp * tmp;
		    }
		}
	    }
	}

      if(m_adj_cost && t<NT)
	{
	  for(s=0; s<NS; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  if(j!=i)
		    {
		      double tmp = 0.0;
		      if(t>0)
			{
			  tmp = e->m2_t[t][i][s][j]/e->m2_t[t-1][i][s][j] - 1.0;
			}
		      else
			{
			  tmp = e->m2_t[t][i][s][j]/p->m02[i][s][j] - 1.0;
			}
		      e->ll_t[t][i] = e->ll_t[t][i] + 
			(p->etaM/2.0) * tmp * tmp;
		    }
		}
	    }
	}

      if(t<NT)
	{
	  if(t==(NT-1) || p->cap_adj_cost==0)
	    {
	      for(s=0; s<NS; s++)
		{
		  e->k_t[t+1][i][s] = (1.0-p->delta) * e->k_t[t][i][s] + e->is_t[t][i][s];
		}
	    }
	  else
	    {
	      for(s=0; s<NS; s++)
		{
		  e->k_t[t+1][i][s] = (1.0-p->delta) * e->k_t[t][i][s] + 
		    phiK(e->is_t[t][i][s]/e->k_t[t][i][s],p->delta,p->gbgp,p->etaK) * e->k_t[t][i][s];
		}
	    }
	}
      else
	{
	  for(s=0; s<NS; s++)
	    {
	      e->is_t[t][i][s] = (p->gbgp-1.0+p->delta) * e->k_t[t][i][s];
	    }
	}
      e->kk_t[t][i] = e->k_t[t][i][0] + e->k_t[t][i][1];
      e->ii_t[t][i] = e->is_t[t][i][0] + e->is_t[t][i][1];

      if(!m_adj_cost || t==NT)
	{
	  for(s=0; s<NS; s++)
	    {
	      e->pm_t[t][i][s] = 0.0;
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_m_ts[t][i][s][j]+p->ntb_m_ts[t][i][s][j];
		  e->pm_t[t][i][s] = e->pm_t[t][i][s] + 
		    pow(p->mu[i][s][j],1.0/(1.0-p->zeta[i][s])) * 
		    pow(tc*e->py_t[t][j][s],p->zeta[i][s]/(p->zeta[i][s]-1.0));
		}
	      e->pm_t[t][i][s] = (1.0/p->M[i][s]) * 
		pow(e->pm_t[t][i][s],(p->zeta[i][s]-1.0)/p->zeta[i][s]);
	    }
	}

      if(!f_adj_cost || t==NT)
	{
	  for(s=0; s<NS; s++)
	    {
	      e->p_t[t][i][s] = 0.0;
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_f_ts[t][i][s][j]+p->ntb_f_ts[t][i][s][j];
		  e->p_t[t][i][s] = e->p_t[t][i][s] + 
		    pow(p->theta[i][s][j],1.0/(1.0-p->sig[i][s])) * 
		    pow(tc*e->py_t[t][j][s],p->sig[i][s]/(p->sig[i][s]-1.0));
		}
	      e->p_t[t][i][s] = (1.0/p->H[i][s]) * 
		pow(e->p_t[t][i][s],(p->sig[i][s]-1.0)/p->sig[i][s]);
	    }
	}

      // households' stochastic discount factor for dynamic firm's problem with adjustment costs
      if(t>0)
	{
	  double mutp = muc(
		e->c_t[t][i],
		e->ll_t[t][i],
		p->lbar_ts[t][i],
		p->pope_ts[t][i],
		p->popw_ts[t][i],
		p->eps[i][0],
		p->rho,
		p->phi[i],
		p->psi,
		0);
	  double mut = muc(
		  e->c_t[t-1][i],
		  e->ll_t[t-1][i],
		  p->lbar_ts[t-1][i],
		  p->pope_ts[t-1][i],
		  p->popw_ts[t-1][i],
		  p->eps[i][0],
		  p->rho,
		  p->phi[i],
		  p->psi,
		  0);

	  e->Q_t[t-1][i] = p->beta[i] * (mutp / e->p_t[t][i][0]) / (mut / e->p_t[t-1][i][0]);
	}
      if(t==NT)
	{
	  e->Q_t[t][i] = p->beta[i];
	}

      e->pi_t[t][i] = pow(e->p_t[t][i][0]/p->eps[i][1][0],p->eps[i][1][0]) * 
	pow(e->p_t[t][i][1]/p->eps[i][1][1],p->eps[i][1][1]) / p->G[i]; 
      e->iy_t[t][i] = e->pi_t[t][i]*e->ii_t[t][i]/e->ngdp_t[t][i];

      for(s=0; s<NS; s++)
	{
	  e->i_t[t][i][s] = e->pi_t[t][i] * p->eps[i][1][s] * e->ii_t[t][i]/e->p_t[t][i][s];
	  e->q_t[t][i][s] = e->c_t[t][i][s] + e->i_t[t][i][s];
	  e->m_t[t][i][s] = e->md_t[t][i][0][s] + e->md_t[t][i][1][s];

	  if(!m_adj_cost || t==NT)
	    {
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_m_ts[t][i][s][j]+p->ntb_m_ts[t][i][s][j];
		  e->m2_t[t][i][s][j] = e->m_t[t][i][s] * 
		    pow(tc*e->py_t[t][j][s],1.0/(p->zeta[i][s]-1.0)) * 
		    pow(e->pm_t[t][i][s]*p->mu[i][s][j]*pow(p->M[i][s],p->zeta[i][s]),1.0/(1.0-p->zeta[i][s]));
		}
	    }

	  if(!f_adj_cost || t==NT)
	    {
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_f_ts[t][i][s][j]+p->ntb_f_ts[t][i][s][j];
		  e->q2_t[t][i][s][j] = e->q_t[t][i][s] * 
		    pow(tc*e->py_t[t][j][s],1.0/(p->sig[i][s]-1.0)) * 
		    pow(e->p_t[t][i][s]*p->theta[i][s][j]*pow(p->H[i][s],p->sig[i][s]),1.0/(1.0-p->sig[i][s]));
		}
	    }
	}

      e->cpi_t[t][i] = (e->p_t[t][i][0]*p->c0[i][0] + e->p_t[t][i][1]*p->c0[i][1])/
	(p->c0[i][0]+p->c0[i][1]);
      e->cc_t[t][i] = e->c_t[t][i][0] + e->c_t[t][i][1];
      
      if(t == 0)
	{
	  e->rk_t[t][i] = p->r0[i] + p->delta;
	}
      else if(t==NT)
	{
	  if(i==0)
	    {
	      e->pb_t[t] = e->cpi_t[t][i]/(1.0+p->rbgp);
	    }
	  if(bgp)
	    {
	      e->rk_t[t][i] = e->pi_t[t][i]*e->cpi_t[t][0]/e->pb_t[t] - (1.0-p->delta)*e->pi_t[t][i];
	    }
	  else
	    {
	      e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->pb_t[t-1] - (1.0-p->delta)*e->pi_t[t][i];
	    }
	}
      else
	{
	  e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->pb_t[t-1] - (1.0-p->delta)*e->pi_t[t][i];
	}
    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  e->ngdp_t[t][i] = e->ngdp_t[t][i] - 
	    e->py_t[t][j][0]*e->m2_t[t][i][0][j] - 
	    e->py_t[t][j][1]*e->m2_t[t][i][1][j];

	  e->rgdp_t[t][i] = e->rgdp_t[t][i] - e->m2_t[t][i][0][j] - e->m2_t[t][i][1][j];

	  if(j!=i)
	    {
	      e->rer_t[t][i][j] = e->cpi_t[t][j]/e->cpi_t[t][i];
	      for(s=0; s<NS; s++)
		{
		  double m = e->py_t[t][i][s]*e->m2_t[t][j][s][i];
		  double f = e->py_t[t][i][s]*e->q2_t[t][j][s][i];
	
		  e->exs_t[t][i][s][j] = m+f;
		  e->exm_t[t][i][j] = e->exm_t[t][i][j] + m;
		  e->exf_t[t][i][j] = e->exf_t[t][i][j] + f;
		  e->ex_t[t][i][j] = e->ex_t[t][i][j] + m+f;

		  m = e->m2_t[t][j][s][i];
		  f = e->q2_t[t][j][s][i];
		  e->rexs_t[t][i][s][j] = m+f;
		  e->rexm_t[t][i][j] = e->rexm_t[t][i][j] + m;
		  e->rexf_t[t][i][j] = e->rexf_t[t][i][j] + f;
		  e->rex_t[t][i][j] = e->rex_t[t][i][j] + m+f;

		  m = e->py_t[t][j][s]*e->m2_t[t][i][s][j];
		  f = e->py_t[t][j][s]*e->q2_t[t][i][s][j];
		  e->ims_t[t][i][s][j] = m+f;
		  e->imm_t[t][i][j] = e->imm_t[t][i][j] + m;
		  e->imf_t[t][i][j] = e->imf_t[t][i][j] + f;
		  e->im_t[t][i][j] = e->im_t[t][i][j] + m+f;

		  m = e->m2_t[t][i][s][j];
		  f = e->q2_t[t][i][s][j];
		  e->rims_t[t][i][s][j] = m+f;
		  e->rimm_t[t][i][j] = e->rimm_t[t][i][j] + m;
		  e->rimf_t[t][i][j] = e->rimf_t[t][i][j] + f;
		  e->rim_t[t][i][j] = e->rim_t[t][i][j] + m+f;

		  e->nxs_t[t][i][s][j] = e->exs_t[t][i][s][j] - e->ims_t[t][i][s][j];
		}
	      e->nxm_t[t][i][j] = e->exm_t[t][i][j] - e->imm_t[t][i][j];
	      e->nxf_t[t][i][j] = e->exf_t[t][i][j] - e->imf_t[t][i][j];
	      e->nx_t[t][i][j] = e->ex_t[t][i][j] - e->im_t[t][i][j];
	    }
	  else
	    {
	      e->rer_t[t][i][j] = 1.0;
	    }
	}

      e->lp_agg_t[t][i] = e->rgdp_t[t][i] / e->ll_t[t][i];
    }  
}

uint eval_bgp_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,s=0,t=NT,nx=0;
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);

  e->b_t[t][0] = bbgp[0];
  e->b_t[t][1] = bbgp[1];
  unstack_bgp_vars(e,myx);
  set_vars(e,p,t,1);
  nx=0;

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  for(i=0; i<(NC-1); i++)
    {
      myf[nx] = bop(p,e,t,i);
      nx = nx+1;
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  myf[nx] = mpk_rk(p,e,t,i,s);
	  nx=nx+1;

	  myf[nx] = mpl_w(p,e,t,i,s);
	  nx=nx+1;

	  myf[nx] = mkt_clear_y(p,e,t,i,s);
	  nx=nx+1;
	}

      myf[nx] = mucg_mucs(p,e,t,i);
      nx=nx+1;

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;
    }

  if(nx != nbgp)
    {
      fprintf(logfile,KRED "Error evaluating bgp eqns! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }
  
  return 0;
}

uint solve_bgp(double bb[NC-1])
{
  bbgp[0] = bb[0];
  bbgp[1] = bb[1];
  solver_n = nbgp;
  alloc_solver_mem();
  set_initial_bgp_guess();
  gsl_multiroot_function_fdf f = {&bgp_func_f,&bgp_func_df,&bgp_func_fdf,nbgp,NULL};
  uint status = find_root_deriv_mkl(&f);

  free_solver_mem();
  write_bgp_vars(&(eee0[0]),&(ppp0[0]));
  return status;
}

uint eval_eqm_conds(const double * myx, double * myf, uint tn)
{
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);
  uint i=0,s=0,t=NT,nx=0;
  uint t0 = 0;

  if(scenario == 2)
    {
      t0 = TREF;
      e = &(eee1[tn]);
      p = &(ppp1[tn]);
    }
  else if(scenario == 3)
    {
      t0 = TREF;
      e = &(eee2[tn]);
      p = &(ppp2[tn]);
    }
  else if(scenario == 4)
    {
      t0 = TVOTE;
      e = &(eee1[tn]);
      p = &(ppp1[tn]);
    }
  else if(scenario == 5)
    {
      t0 = TVOTE;
      e = &(eee2[tn]);
      p = &(ppp2[tn]);
    }  

  unstack_eqm_vars(e,myx);

  e->b_t[0][0] = p->b0[0];
  e->b_t[0][1] = p->b0[1];
  e->b_t[0][2] = p->b0[2];

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  e->k_t[0][i][s] = p->k0[i][s];
	}
    }

  for(t=t0; t<(NT+1); t++)
    {
      set_vars(e,p,t,0);
    }

  nx=0;
  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<(NC-1); i++)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(t<NT)
		{
		  myf[nx] = mpk_rk(p,e,t+1,i,s);
		  nx=nx+1;
		}

	      myf[nx] = mpl_w(p,e,t,i,s);
	      nx=nx+1;

	      myf[nx] = mkt_clear_y(p,e,t,i,s);
	      nx=nx+1;
	    }

	  myf[nx] = mucg_mucs(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;
	    }

	  if(m_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_m_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_m2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }

	  if(f_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_q_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_q2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error evaluating eqm eqns! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }
  
  return 0;
}

uint solve_eqm()
{
  char * sname;

#ifdef _QUARTERLY
  if(f_adj_cost==0 && m_adj_cost==0 && tfp_flag==0 && fixl!=2)
    {
      sname = "output/seed_quarterly.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==0)
    {
      sname = "output/seed_ffrictions_quarterly.bin";
    }
  else if(f_adj_cost==0 && m_adj_cost==1)
    {
      sname = "output/seed_mfrictions_quarterly.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==1)
    {
      sname = "output/seed_frictions_quarterly.bin";
    }
  else if(fixl==2)
    {
      sname = "output/seed_stickyw_quarterly.bin";
    }
  else if(tfp_flag==1)
    {
      sname = "output/seed_tfp_quarterly.bin";
    }
#else
  if(f_adj_cost==0 && m_adj_cost==0 && tfp_flag==0 && fixl!=2)
    {
      sname = "output/seed.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==0)
    {
      sname = "output/seed_ffrictions.bin";
    }
  else if(f_adj_cost==0 && m_adj_cost==1)
    {
      sname = "output/seed_mfrictions.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==1)
    {
      sname = "output/seed_frictions.bin";
    }
  else if(fixl==2)
    {
      sname = "output/seed_stickyw.bin";
    }
  else if(tfp_flag==1)
    {
      sname = "output/seed_tfp.bin";
    }
#endif

  // if we are solving for the no-Brexit counterfactual we must construct an initial guess for the solver
  if(scenario==0)
    {
      if(read_seed==1)
	{
	  free_solver_mem();
	  solver_n = neqm;
	  alloc_solver_mem();

	  if(read_vec_bin(solver_x->data, neqm, sname))
	    {
	      fprintf(logfile,KRED "Error loading equilibrium guess from seed file!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
      else
	{
	  if(set_initial_eqm_guess())
	    {
	      fprintf(logfile,KRED "Error constructing equilibrium guess!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
    }
  // otherwise we should use the solution from the previous exercise as the initial guess
  else
    {      
      free_solver_mem();
      solver_n = neqm;
      alloc_solver_mem();

      stoch_eqm * se = &(sss[0]);
      eqm * e;

      if(scenario == 2)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1[it]), &(eee0[0])  );
	    }
	}
      else if(scenario == 3)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee2[it]), &(eee0[0])  );
	    }
	}
      else if(scenario == 4)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1[it]), &((se->eee)[0]) );
	    }
	}
      else if(scenario == 5)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee2[it]), &((se->eee)[0]) );
	    }
	}

      // now determine which eqm struct to use as the initial guess...
      // ...if we are in the optimistic scenario, use the first path...
      if(scenario == 2 || scenario == 4)
	{
	  e = &((se->eee)[1]);
	}
      // ...otherwise use the second path
      else if(scenario == 3 || scenario == 5)
	{
	  e = &((se->eee)[2]);
	}
      if(stack_eqm_vars(solver_x->data,e))
	{
	  fprintf(logfile,KRED "Failed to stack variables from previous exercise!\n" RESET);
	  free_solver_mem();
	  return 1;
	}
    }

  uint status = 0;
  if(eval_eqm_once_flag)
    {
      status = eqm_func_f(solver_x,NULL,f0[0]);
#ifdef _QUARTERLY
      write_vec_txt(f0[0]->data,solver_n,"output/F_quarterly.txt");
#else
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
#endif
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      gsl_multiroot_function_fdf f = {&eqm_func_f,&eqm_func_df,&eqm_func_fdf,neqm,NULL};
      status = find_root_deriv_mkl(&f);
      if(status)
	fprintf(logfile,KRED "Error solving for equilibrium!\n" RESET);

      if(scenario==0 && write_seed==1 && !status)
	{
	  write_vec_bin(solver_x->data, neqm, sname);
	}

    }

  free_solver_mem();

  return status;
}

uint eval_stoch_eqm_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,s=0,t=NT,nx=0,ih;

  stoch_eqm * se = &(sss[tn]);
  stoch_params * sp = &(sppp[tn]);
  eqm *e, *e2;
  params *p, *p2;

  unstack_stoch_eqm_vars(se,myx);

  uint t0 = TREF;

  for(ih=0; ih<NHIST; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      for(t=t0; t<(NT+1); t++)
	{
	  set_vars(e,p,t,0);
	}
    }

  nx=0;
  
  // first we do the equations for the periods after referendum announcement but before the vote
  // during these periods, we only have to evaluate equilibrium conditions for one history
  // (all histories are the same up until the vote)
  //
  // there is no uncertainty except in the period immediately before the vote, so we do all the usual equations
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);
  for(t=t0; t<(TVOTE-1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<(NC-1); i++)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(t<NT)
		{
		  myf[nx] = mpk_rk(p,e,t+1,i,s);
		  nx=nx+1;
		}

	      myf[nx] = mpl_w(p,e,t,i,s);
	      nx=nx+1;

	      myf[nx] = mkt_clear_y(p,e,t,i,s);
	      nx=nx+1;
	    }

	  myf[nx] = mucg_mucs(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;
	    }

	  if(m_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_m_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_m2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }

	  if(f_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_q_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_q2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }
	}
    }

  // in the period before the vote, we have to deal with the uncertainty about whether we will get
  // a "stay" vote (with probability pi_vote) or a "leave" vote (with probability 1-pi_vote)
  t = TVOTE-1;
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);
  e2 = &(se->eee[1]);
  p2 = &(sp->ppp[1]);

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  // no need to worry about uncertainty with the balance of payments... b_t[i][t+1] is determined in period t
  // so it is the same in all histories
  for(i=0; i<(NC-1); i++)
    {
      myf[nx] = bop(p,e,t,i);
      nx = nx+1;
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  // since the MPK equation is forward looking, we have uncertainty!
	  // fortunately with the code structure it is very easy:
	  myf[nx] = pi_vote * mpk_rk(p,e,t+1,i,s) + (1.0-pi_vote) * mpk_rk(p2,e2,t+1,i,s);
	  nx=nx+1;

	  myf[nx] = mpl_w(p,e,t,i,s);
	  nx=nx+1;

	  myf[nx] = mkt_clear_y(p,e,t,i,s);
	  nx=nx+1;
	}

      myf[nx] = mucg_mucs(p,e,t,i);
      nx=nx+1;

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      // more uncertainty here:
      myf[nx] = pi_vote * euler(p,e,t,i) + (1.0-pi_vote) * euler(p2,e2,t,i);
      nx = nx+1;

      if(m_adj_cost)
	{
	  for(s=0; s<NS; s++)
	    {
	      myf[nx] = prod_m_chk(p,e,t,i,s);
	      nx = nx+1;
	      
	      // even more here!
	      uint j;
	      for(j=0; j<NC; j++)
		{
		  myf[nx] = pi_vote * foc_m2(p,e,t,i,s,j) + (1.0-pi_vote) * foc_m2(p2,e2,t,i,s,j);
		  nx = nx+1;
		}
	    }
	}

      if(f_adj_cost)
	{
	  for(s=0; s<NS; s++)
	    {
	      myf[nx] = prod_q_chk(p,e,t,i,s);
	      nx = nx+1;
	      
	      // even more here!
	      uint j;
	      for(j=0; j<NC; j++)
		{
		  myf[nx] = pi_vote * foc_q2(p,e,t,i,s,j) + (1.0-pi_vote) * foc_q2(p2,e2,t,i,s,j);
		  nx = nx+1;
		}
	    }
	}
    }

  // next we do the equations for the "stay" path... no uncertainty here and it goes all the way to the steady state
  t0 = TVOTE;
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);

  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<(NC-1); i++)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(t<NT)
		{
		  myf[nx] = mpk_rk(p,e,t+1,i,s);
		  nx=nx+1;
		}

	      myf[nx] = mpl_w(p,e,t,i,s);
	      nx=nx+1;

	      myf[nx] = mkt_clear_y(p,e,t,i,s);
	      nx=nx+1;
	    }

	  myf[nx] = mucg_mucs(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;
	    }

	  if(m_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_m_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_m2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }

	  if(f_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_q_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_q2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }
	}
    }
    
  // now we do the TBREXIT-TVOTE periods between vote and Brexit... here there is uncertainty again in the
  // period before Brexit only
  t0 = TVOTE;
  e = &(se->eee[1]);
  p = &(sp->ppp[1]);
  for(t=t0; t<(TBREXIT-1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<(NC-1); i++)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(t<NT)
		{
		  myf[nx] = mpk_rk(p,e,t+1,i,s);
		  nx=nx+1;
		}

	      myf[nx] = mpl_w(p,e,t,i,s);
	      nx=nx+1;

	      myf[nx] = mkt_clear_y(p,e,t,i,s);
	      nx=nx+1;
	    }

	  myf[nx] = mucg_mucs(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;
	    }

	  if(m_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_m_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_m2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }

	  if(f_adj_cost && t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  myf[nx] = prod_q_chk(p,e,t,i,s);
		  nx = nx+1;

		  uint j;
		  for(j=0; j<NC; j++)
		    {
		      myf[nx] = foc_q2(p,e,t,i,s,j);
		      nx = nx+1;
		    }
		}
	    }
	}
    }

  // in the period before the Brexit, we have to deal with the uncertainty about whether we will get
  // hard or soft Brexit
  t = TBREXIT-1;
  e = &(se->eee[1]);
  p = &(sp->ppp[1]);
  e2 = &(se->eee[2]);
  p2 = &(sp->ppp[2]);

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  for(i=0; i<(NC-1); i++)
    {
      myf[nx] = bop(p,e,t,i);
      nx = nx+1;
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  // since the MPK equation is forward looking, we have uncertainty!
	  myf[nx] = pi_brexit * mpk_rk(p,e,t+1,i,s) + (1.0-pi_brexit) * mpk_rk(p2,e2,t+1,i,s);
	  nx=nx+1;

	  myf[nx] = mpl_w(p,e,t,i,s);
	  nx=nx+1;

	  myf[nx] = mkt_clear_y(p,e,t,i,s);
	  nx=nx+1;
	}

      myf[nx] = mucg_mucs(p,e,t,i);
      nx=nx+1;

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      // more uncertainty here:
      myf[nx] = pi_brexit * euler(p,e,t,i) + (1.0-pi_brexit) * euler(p2,e2,t,i);
      nx = nx+1;

      if(m_adj_cost)
	{
	  for(s=0; s<NS; s++)
	    {
	      myf[nx] = prod_m_chk(p,e,t,i,s);
	      nx = nx+1;
	      
	      // even more here!
	      uint j;
	      for(j=0; j<NC; j++)
		{
		  myf[nx] = pi_brexit * foc_m2(p,e,t,i,s,j) + (1.0-pi_brexit) * foc_m2(p2,e2,t,i,s,j);
		  nx = nx+1;
		}
	    }
	}

      if(f_adj_cost)
	{
	  for(s=0; s<NS; s++)
	    {
	      myf[nx] = prod_q_chk(p,e,t,i,s);
	      nx = nx+1;
	      
	      // even more here!
	      uint j;
	      for(j=0; j<NC; j++)
		{
		  myf[nx] = pi_brexit * foc_q2(p,e,t,i,s,j) + (1.0-pi_brexit) * foc_q2(p2,e2,t,i,s,j);
		  nx = nx+1;
		}
	    }
	}
    }

  // after Brexit, we have to worry about the separate possible histories
  // but there is no uncertainty to worry about... everything is deterministic
  // within a given history from Brexit onwards
  t0 = TBREXIT;
  for(ih=1; ih<NHIST; ih++)
    {
      e = &(se->eee[ih]);
      p = &(sp->ppp[ih]);

      for(t=t0; t<(NT+1); t++)
	{
	  myf[nx] = price_norm(e,t);
	  nx=nx+1;

	  for(i=0; i<(NC-1); i++)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  for(i=0; i<NC; i++)
	    {
	      for(s=0; s<NS; s++)
		{
		  if(t<NT)
		    {
		      myf[nx] = mpk_rk(p,e,t+1,i,s);
		      nx=nx+1;
		    }

		  myf[nx] = mpl_w(p,e,t,i,s);
		  nx=nx+1;

		  myf[nx] = mkt_clear_y(p,e,t,i,s);
		  nx=nx+1;
		}

	      myf[nx] = mucg_mucs(p,e,t,i);
	      nx=nx+1;

	      myf[nx] = muc_mul(p,e,t,i);
	      nx=nx+1;

	      if(t<NT)
		{
		  myf[nx] = euler(p,e,t,i);
		  nx = nx+1;
		}

	      if(m_adj_cost && t<NT)
		{
		  for(s=0; s<NS; s++)
		    {
		      myf[nx] = prod_m_chk(p,e,t,i,s);
		      nx = nx+1;

		      uint j;
		      for(j=0; j<NC; j++)
			{
			  myf[nx] = foc_m2(p,e,t,i,s,j);
			  nx = nx+1;
			}
		    }
		}

	      if(f_adj_cost && t<NT)
		{
		  for(s=0; s<NS; s++)
		    {
		      myf[nx] = prod_q_chk(p,e,t,i,s);
		      nx = nx+1;

		      uint j;
		      for(j=0; j<NC; j++)
			{
			  myf[nx] = foc_q2(p,e,t,i,s,j);
			  nx = nx+1;
			}
		    }
		}
	    }
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error evaluating stochastic eqm eqns! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }
  
  return 0;
}

uint solve_stoch_eqm()
{  
  char * sname;
#ifdef _QUARTERLY
  if(f_adj_cost==0 && m_adj_cost==0 && tfp_flag==0 && fixl!=2)
    {
      sname = "output/stoch_seed_quarterly.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==0)
    {
      sname = "output/stoch_seed_ffrictions_quarterly.bin";
    }
  else if(f_adj_cost==0 && m_adj_cost==1)
    {
      sname = "output/stoch_seed_mfrictions_quarterly.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==1)
    {
      sname = "output/stoch_seed_frictions_quarterly.bin";
    }
  else if(fixl==2)
    {
      sname = "output/stoch_seed_stickyw_quarterly.bin";
    }
  else if(tfp_flag==1)
    {
      sname = "output/stoch_seed_tfp_quarterly.bin";
    }
#else
  if(f_adj_cost==0 && m_adj_cost==0 && tfp_flag==0 && fixl!=2)
    {
      sname = "output/stoch_seed.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==0)
    {
      sname = "output/stoch_seed_ffrictions.bin";
    }
  else if(f_adj_cost==0 && m_adj_cost==1)
    {
      sname = "output/stoch_seed_mfrictions.bin";
    }
  else if(f_adj_cost==1 && m_adj_cost==1)
    {
      sname = "output/stoch_seed_frictions.bin";
    }
  else if(fixl==2)
    {
      sname = "output/stoch_seed_stickyw.bin";
    }
  else if(tfp_flag==1)
    {
      sname = "output/stoch_seed_tfp.bin";
    }
#endif

  // otherwise we should use the solution from the no-Brexit counterfactual as the initial guess
  // first we need to copy everything from eee[0] to the stochastic equilibrium sss[it] for each
  // thread it (we have to do it for all thread because otherwise those threads' variables for
  // t = 0 through t = TREF will never get set
  uint it,ih;
  for(it=0; it<NTH; it++)
    {
      for(ih=0; ih<NHIST; ih++)
	{
	  copy_vars( &( (&(sss[it]))->eee[ih] ) , &(eee0[0]) );
	}
    }
  
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

  if(read_stoch_seed)
    {
      if(read_vec_bin(solver_x->data, neqm, sname))
	{
	  fprintf(logfile,KRED "Error loading stochastic equilibrium guess from seed file!\n" RESET);
	  free_solver_mem();
	  return 1;
	}
    }
  else
    {
      if(stack_stoch_eqm_vars(solver_x->data,&(sss[0])))
	{
	  fprintf(logfile,KRED "Failed to stack variables from previous exercise!\n" RESET);
	  free_solver_mem();
	  return 1;
	}
    }

  // now we solve the stochastic model!!
  uint status = 0;
  if(eval_eqm_once_flag)
    {
      status = stoch_eqm_func_f(solver_x,NULL,f0[0]);
#ifdef _QUARTERLY
      write_vec_txt(f0[0]->data,solver_n,"output_quarterly/F.txt");
#else
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
#endif
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      gsl_multiroot_function_fdf f = {&stoch_eqm_func_f,&stoch_eqm_func_df,&stoch_eqm_func_fdf,neqm,NULL};
      status = find_root_deriv_mkl(&f);
      if(status)
	fprintf(logfile,KRED "Error solving for stochastic equilibrium!\n" RESET);

      if(!status && write_stoch_seed==1)
	{
	  write_vec_bin(solver_x->data, neqm, sname);
	}
    }

  free_solver_mem();

  return status;
}

void calc_welfare(eqm * e, const params * p)
{
  int t, i;
  for(i=0; i<NC; i++)
    {
      t=NT;
      e->welfare_t[t][i] = (1.0/(1.0-p->beta[i]*pow(p->gbgp,p->phi[i]*p->psi))) * 
	(pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
	     p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
	     p->phi[i]*p->psi/p->rho) * 
	 pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));	
      
      for(t=(NT-1); t>=0; t--)
	{
	  e->welfare_t[t][i] = p->beta[i] * e->welfare_t[t+1][i] + 
	    (pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
		 p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
		 p->phi[i]*p->psi/p->rho) * 
	     pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
	  
	  e->welfare_t[t+1][i] = pow(e->welfare_t[t+1][i],1.0/p->psi);
	}
      e->welfare_t[0][i] = pow(e->welfare_t[0][i],1.0/p->psi);
    }

  if(scenario == 0)
    {
      for(i=0; i<NC; i++)
	{
	  t=NT;
	  e->welfare_cost_t[t][i] = (1.0/(1.0-e->pb_t[t]*p->gbgp)) * e->cpi_t[t][i] * e->cc_t[t][i];
	  
	  for(t=(NT-1); t>=0; t--)
	    {
	      e->welfare_cost_t[t][i] = e->cpi_t[t][i]*e->cc_t[t][i] + e->pb_t[t] * e->welfare_cost_t[t+1][i];
	    }
	  for(t=0; t<NT; t++)
	    {
	      e->welfare_cost_t[t][i] = e->pb_t[t] * e->welfare_cost_t[t+1][i];
	    }

 	}
    }
}
 
void calc_stoch_welfare()
{
  int t, i, ih;
  eqm *e, *e0, *e2;
  const params *p;

  stoch_eqm * se = &(sss[0]);
  const stoch_params * sp = &(sppp[0]);

  // backwords from steady state to period of Brexit
  for(ih=0; ih<NHIST; ih++)
    { 
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(i=0; i<NC; i++)
	{
	  t=NT;
	  e->welfare2_t[t][i] = (1.0/(1.0-p->beta[i]*pow(p->gbgp,p->phi[i]*p->psi))) * 
	    (pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
		 p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
		 p->phi[i]*p->psi/p->rho) * 
	     pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));	
      
	  for(t=(NT-1); t>=TBREXIT; t--)
	    {
	      e->welfare2_t[t][i] = p->beta[i] * e->welfare2_t[t+1][i] + 
		(pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
		     p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
		     p->phi[i]*p->psi/p->rho) * 
		 pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
	  
	      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	}
    }

  // last pre-Brexit period for "stay" branch
  t=TBREXIT-1;
  e = &( (se->eee)[0] );
  p = &( (sp->ppp)[0] );
  for(i=0; i<NC; i++)
    {      
      e->welfare2_t[t][i] = p->beta[i] * e->welfare2_t[t+1][i] + 
	(pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
	     p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
	     p->phi[i]*p->psi/p->rho) * 
	 pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
      
      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
    }

  // last pre-Brexit period for "leave" branch
  t = TBREXIT-1;
  e = &( (se->eee)[1] );
  p = &( (sp->ppp)[1] );
  e2 = &( (se->eee)[2] );
  for(i=0; i<NC; i++)
    {
      e->welfare2_t[t][i] = (pi_brexit * p->beta[i] * e->welfare2_t[t+1][i]) + 
	(1.0-pi_brexit) * (p->beta[i]*e2->welfare2_t[t+1][i]) + 
	(pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
	     p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
	     p->phi[i]*p->psi/p->rho) * 
	 pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));

      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
      e2->welfare2_t[t+1][i] = pow(e2->welfare2_t[t+1][i],1.0/p->psi);
    }

  // backwords from TBREXIT-2 to vote period
  for(ih=0; ih<NHIST; ih++)
    { 
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      
      for(i=0; i<NC; i++)
	{
	  for(t=TBREXIT-2; t>=TVOTE; t--)
	    {
	      e->welfare2_t[t][i] = p->beta[i] * e->welfare2_t[t+1][i] + 
		(pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
		     p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
		     p->phi[i]*p->psi/p->rho) * 
		 pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
	      
	      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	}
    }

  // vote period (all branches)
  // Brexit period for "leave" branch
  t = TVOTE-1;
  e0 = &( (se->eee)[0] );
  e2 = &( (se->eee)[1] );
  for(ih=0; ih<NC; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      
      for(i=0; i<NC; i++)
	{
	  e->welfare2_t[t][i] = (pi_vote * p->beta[i] * e0->welfare2_t[t+1][i]) + 
	    (1.0-pi_vote) * (p->beta[i]*e2->welfare2_t[t+1][i]) + 
	    (pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
		 p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
		 p->phi[i]*p->psi/p->rho) * 
	     pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
	  
	  //e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	  //e2->welfare2_t[t+1][i] = pow(e2->welfare2_t[t+1][i],1.0/p2->psi);
	}
    }
  for(ih=0; ih<NC; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      for(i=0; i<NC; i++)
	{
	  e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	}
    }

  // backwards from TVOTE-1 to 0
  for(ih=0; ih<NHIST; ih++)
    { 
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(i=0; i<NC; i++)
	{      
	  for(t=(TVOTE-2); t>=0; t--)
	    {
	      e->welfare2_t[t][i] = p->beta[i] * e->welfare2_t[t+1][i] + 
		(pow(p->eps[i][0][0] * pow(e->c_t[t][i][0]/p->pope_ts[t][i],p->rho) + 
		     p->eps[i][0][1] * pow(e->c_t[t][i][1]/p->pope_ts[t][i],p->rho),
		     p->phi[i]*p->psi/p->rho) * 
		 pow((p->lbar_ts[t][i]-e->ll_t[t][i])/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
	  
	      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	  e->welfare2_t[0][i] = pow(e->welfare2_t[0][i],1.0/p->psi);
	}
    }
}

void calc_ce_welfare1()
{
  int t, i;

  stoch_eqm *se = &(sss[0]);
  const stoch_params *sp = &(sppp[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[1] );
  const eqm *e_stay = &( (se->eee)[0] );
  const eqm *e_soft = &(eee1[0]);
  const eqm * e_hard = &(eee2[0]);

  double c0_stay, c0_soft, c0_hard, c0;
  double c1_stay, c1_soft, c1_hard, c1;
  double ll_stay, ll_soft, ll_hard, ll;

  for(i=0; i<NC; i++)
    {
      t=NT;

      c0_stay = e_stay->c_t[t][i][0];
      c1_stay = e_stay->c_t[t][i][1];
      ll_stay = e_stay->ll_t[t][i];

      c0_soft = e_soft->c_t[t][i][0];
      c1_soft = e_soft->c_t[t][i][1];
      ll_soft = e_soft->ll_t[t][i];

      c0_hard = e_hard->c_t[t][i][0];
      c1_hard = e_hard->c_t[t][i][1];
      ll_hard = e_hard->ll_t[t][i];

      c0 = pi_vote*c0_stay +
	(1.0-pi_vote)*pi_brexit*c0_soft + 
	(1.0-pi_vote)*(1.0-pi_brexit)*c0_hard;

      c1 = pi_vote*c1_stay +
	(1.0-pi_vote)*pi_brexit*c1_soft + 
	(1.0-pi_vote)*(1.0-pi_brexit)*c1_hard;

      ll = pi_vote*ll_stay +
	(1.0-pi_vote)*pi_brexit*ll_soft + 
	(1.0-pi_vote)*(1.0-pi_brexit)*ll_hard;

      e->welfare3_t[t][i] = (1.0/(1.0-p->beta[i]*pow(p->gbgp,p->phi[i]*p->psi))) * 
	(pow(p->eps[i][0][0] * pow(c0/p->pope_ts[t][i],p->rho) + 
	     p->eps[i][0][1] * pow(c1/p->pope_ts[t][i],p->rho),
	     p->phi[i]*p->psi/p->rho) * 
	 pow((p->lbar_ts[t][i]-ll)/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));	
      
      for(t=(NT-1); t>=TREF; t--)
	{
	  c0_stay = e_stay->c_t[t][i][0];
	  c1_stay = e_stay->c_t[t][i][1];
	  ll_stay = e_stay->ll_t[t][i];

	  c0_soft = e_soft->c_t[t][i][0];
	  c1_soft = e_soft->c_t[t][i][1];
	  ll_soft = e_soft->ll_t[t][i];

	  c0_hard = e_hard->c_t[t][i][0];
	  c1_hard = e_hard->c_t[t][i][1];
	  ll_hard = e_hard->ll_t[t][i];

	  c0 = pi_vote*c0_stay +
	    (1.0-pi_vote)*pi_brexit*c0_soft + 
	    (1.0-pi_vote)*(1.0-pi_brexit)*c0_hard;

	  c1 = pi_vote*c1_stay +
	    (1.0-pi_vote)*pi_brexit*c1_soft + 
	    (1.0-pi_vote)*(1.0-pi_brexit)*c1_hard;

	  ll = pi_vote*ll_stay +
	    (1.0-pi_vote)*pi_brexit*ll_soft + 
	    (1.0-pi_vote)*(1.0-pi_brexit)*ll_hard;

	  e->welfare3_t[t][i] = p->beta[i] * e->welfare3_t[t+1][i] + 
	    (pow(p->eps[i][0][0] * pow(c0/p->pope_ts[t][i],p->rho) + 
		 p->eps[i][0][1] * pow(c1/p->pope_ts[t][i],p->rho),
		 p->phi[i]*p->psi/p->rho) * 
	     pow((p->lbar_ts[t][i]-ll)/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
	  
	  e->welfare3_t[t+1][i] = pow(e->welfare3_t[t+1][i],1.0/p->psi);
	}
      e->welfare3_t[TREF][i] = pow(e->welfare3_t[TREF][i],1.0/p->psi);
    }
}

void calc_ce_welfare2()
{
  int t, i;

  stoch_eqm *se = &(sss[0]);
  const stoch_params *sp = &(sppp[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[1] );
  const eqm *e_soft = &(eee1[0]);
  const eqm * e_hard = &(eee2[0]);

  double c0_soft, c0_hard, c0;
  double c1_soft, c1_hard, c1;
  double ll_soft, ll_hard, ll;

  for(i=0; i<NC; i++)
    {
      t=NT;

      c0_soft = e_soft->c_t[t][i][0];
      c1_soft = e_soft->c_t[t][i][1];
      ll_soft = e_soft->ll_t[t][i];

      c0_hard = e_hard->c_t[t][i][0];
      c1_hard = e_hard->c_t[t][i][1];
      ll_hard = e_hard->ll_t[t][i];

      c0 = pi_brexit*c0_soft + (1.0-pi_brexit)*c0_hard;
      c1 = pi_brexit*c1_soft + (1.0-pi_brexit)*c1_hard;
      ll = pi_brexit*ll_soft + (1.0-pi_brexit)*ll_hard;

      e->welfare4_t[t][i] = (1.0/(1.0-p->beta[i]*pow(p->gbgp,p->phi[i]*p->psi))) * 
	(pow(p->eps[i][0][0] * pow(c0/p->pope_ts[t][i],p->rho) + 
	     p->eps[i][0][1] * pow(c1/p->pope_ts[t][i],p->rho),
	     p->phi[i]*p->psi/p->rho) * 
	 pow((p->lbar_ts[t][i]-ll)/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));	
      
      for(t=(NT-1); t>=TVOTE; t--)
	{
	  c0_soft = e_soft->c_t[t][i][0];
	  c1_soft = e_soft->c_t[t][i][1];
	  ll_soft = e_soft->ll_t[t][i];

	  c0_hard = e_hard->c_t[t][i][0];
	  c1_hard = e_hard->c_t[t][i][1];
	  ll_hard = e_hard->ll_t[t][i];

	  c0 = pi_brexit*c0_soft + (1.0-pi_brexit)*c0_hard;
	  c1 = pi_brexit*c1_soft + (1.0-pi_brexit)*c1_hard;
	  ll = pi_brexit*ll_soft + (1.0-pi_brexit)*ll_hard;

	  e->welfare4_t[t][i] = p->beta[i] * e->welfare4_t[t+1][i] + 
	    (pow(p->eps[i][0][0] * pow(c0/p->pope_ts[t][i],p->rho) + 
		 p->eps[i][0][1] * pow(c1/p->pope_ts[t][i],p->rho),
		 p->phi[i]*p->psi/p->rho) * 
	     pow((p->lbar_ts[t][i]-ll)/p->popw_ts[t][i],(1.0-p->phi[i])*p->psi));
	  
	  e->welfare4_t[t+1][i] = pow(e->welfare4_t[t+1][i],1.0/p->psi);
	}
      e->welfare4_t[TVOTE][i] = pow(e->welfare4_t[TVOTE][i],1.0/p->psi);
    }
}


///////////////////////////////

int bgp_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
   //fcnt = fcnt + 1;
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_bgp_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int bgp_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&bgp_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int bgp_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(bgp_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&bgp_func_f, x, J, 0))
	{
	  return GSL_EFAILED;
	}
      else
	{
	  return GSL_SUCCESS;
	}
    }
}

int eqm_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
  //fcnt = fcnt + 1;
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_eqm_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int eqm_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&eqm_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int eqm_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(eqm_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&eqm_func_f, x, J, 0))
	{
	  return GSL_EFAILED;
	}
      else
	{
	  return GSL_SUCCESS;
	}
     }
}

int stoch_eqm_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
  //fcnt = fcnt + 1;
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_stoch_eqm_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int stoch_eqm_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&stoch_eqm_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int stoch_eqm_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(stoch_eqm_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&stoch_eqm_func_f, x, J, 0))
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
