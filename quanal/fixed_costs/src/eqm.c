#ifndef __EQM_C__
#define __EQM_C__

#include "eqm.h"

// wage (NC), capital (NC), gross output (NC), bilateral prices (NC*NC)
const uint nbgp = NC*3 + NC*NC;

// key eqm vars used in solver: wage (NC), gross output (NC), investment (NC, but not in last period),
// biateral prices (NC*NC), bonds (NC-1), but not in first period, rental rates on capital (NC, not in last period),
// bond prices (NC, not in last period)
void set_neqm()
{
  uint nn = 0;
  if(scenario != 1 && scenario != 9)
    {
      if(scenario==0 || scenario==6)
	{
	  nn = NT+1;
	}
      else if(scenario==2 || scenario==3)
	{
	  nn = NT+1 - TREF;
	}
      else if(scenario==4 || scenario==5)
	{
	  nn = NT+1 - TVOTE;
	}
      else if(scenario==7 || scenario==8 || (scenario>=1000 && scenario<2000))
	{
	  nn = NT+1 - TCH0;
	}
      else if(scenario==10 || scenario==11 || scenario>=2000)
	{
	  nn = NT+1 - TCH1;
	}

      neqm = nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;
    }
  else if(scenario==1)
    {
      // we have TVOTE-TREF periods with only one branch between period in which referendum is announced
      // and period in which vote takes place
      nn = TVOTE-TREF;
      neqm = nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;

      // in the period of the vote, the tree splits...
      
      // along branch 0, the vote fails and Brexit never occurs... this is a deterministic branch
      // that goes all the way to the steady state
      nn = (NT+1) - TVOTE;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;

      // along branch 1, the vote goes through and we have to wait until TBREXIT to find out whether
      // we have hard or soft Brexit... this section is TBREXIT-TVOTE periods long
      nn = TBREXIT-TVOTE;
      neqm += nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;

      // after Brexit, we have 2 possible histories that can arise, represented by two distinct
      // deterministic sections
      nn = (NT+1) - TBREXIT;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;  
    }
  else if(scenario==9)
    {
#ifdef CH_SIMPLE
      nn = TCH1-TCH0;
      neqm = nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;      

      nn = (NT+1)-TCH1;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;

      nn = TCH2-TCH1;
      neqm += nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;

      nn = (NT+1) - TCH2;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;        
#else
      // periods after announcement but before implementation
      nn = TCH1-TCH0;
      neqm = nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;      

      // "no implementation" branch after separation
      nn = (NT+1)-TCH1;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;

      // all other branches with reversion (or not, for ih==NHIST_CH-1) after separation
      uint ih;
      for(ih=1; ih<NHIST_CH; ih++)
	{
	  nn = NT+1-(TCH1+ih);
	  if(ih==NHIST_CH-1)
	    {
	      nn = NT+1-(TCH1+ih-1);
	    }
	  neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;
	}

      // periods with uncertainty (parts where 1 or more branches overlaps with last branch
      nn = TCH2-TCH1;
      neqm += nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;
#endif
    }
}

void init_vars(eqm * e)
{
  SET_ALL_V(e->Q_t,(NT+1)*NC,0.0);

  SET_ALL_V(e->B_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->C_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->MUC_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->X_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->L_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->K_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->Y_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->M_t,(NT+1)*NC,0.0);

  SET_ALL_V(e->P_t,(NT+1)*NC,1.0);
  SET_ALL_V(e->W_t,(NT+1)*NC,1.0);
  SET_ALL_V(e->R_t,(NT+1)*NC,1.0);
  SET_ALL_V(e->Lambda_t,(NT+1)*NC,1.0);
  SET_ALL_V(e->MC_t,(NT+1)*NC,1.0);

  SET_ALL_V(e->NGDP_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->RGDP_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->XY_t,(NT+1)*NC,0.0);

  SET_ALL_V(e->EX_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->IM_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->NX_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->RER_t,(NT+1)*NC*NC,1.0);
  SET_ALL_V(e->TOT_t,(NT+1)*NC*NC,1.0);

  SET_ALL_V(e->P2_t,(NT+1)*NC*NC,1.0);
  SET_ALL_V(e->Y2_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->Y2s_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->Kd_t,(NT+1)*NC,0.0);

  SET_ALL_V(e->D_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->Db_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->Dh_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->Dt_t,(NT+1)*NC*NC,0.0);

  uint t;
  for(t=0; t<NT+1; t++)
    {
      uint i;
      for(i=0; i<NC; i++)
	{
	  uint ii;
	  for(ii=0; ii<2; ii++)
	    {
	      e->ev_t[t][i][ii].zp=0.0;
	      e->ev_t[t][i][ii].zm=0.0;
	      e->ev_t[t][i][ii].Zi=0.0;
	      e->ev_t[t][i][ii].Zp=0.0;
	      e->ev_t[t][i][ii].Zm=0.0;
	      e->ev_t[t][i][ii].Z=0.0;
	      e->ev_t[t][i][ii].Fzp=0.0;
	      e->ev_t[t][i][ii].Fzm=0.0;
	      e->ev_t[t][i][ii].dV=0.0;
	      e->ev_t[t][i][ii].n=0.0;
	    }
	}
    }
  
  SET_ALL_V(e->welfare_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->welfare2_t,(NT+1)*NC,0.0); 
  SET_ALL_V(e->welfare3_t,(NT+1)*NC,0.0); 
  SET_ALL_V(e->welfare4_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->welfare_cost_t,(NT+1)*NC,0.0);
}

void copy_vars(eqm * e1, const eqm * e0)
{
  memcpy((double *)(e1->Q_t),(const double *)(e0->Q_t),sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->B_t),(const double *)(e0->B_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->C_t),(const double *)(e0->C_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->MUC_t),(const double *)(e0->C_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->X_t),(const double *)(e0->X_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->Lf_t),(const double *)(e0->Lf_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->L_t),(const double *)(e0->L_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->K_t),(const double *)(e0->K_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->Y_t),(const double *)(e0->Y_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->M_t),(const double *)(e0->M_t),sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->P_t),(const double *)(e0->P_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->W_t),(const double *)(e0->W_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->R_t),(const double *)(e0->R_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->Lambda_t),(const double *)(e0->Lambda_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->MC_t),(const double *)(e0->MC_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->NGDP_t),(const double *)(e0->NGDP_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->RGDP_t),(const double *)(e0->RGDP_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->XY_t),(const double *)(e0->XY_t),sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->EX_t),(const double *)(e0->EX_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->IM_t),(const double *)(e0->IM_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->EXR_t),(const double *)(e0->EX_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->IMR_t),(const double *)(e0->IM_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->NX_t),(const double *)(e0->NX_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->RER_t),(const double *)(e0->RER_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->TOT_t),(const double *)(e0->TOT_t),sizeof(double)*(NT+1)*NC*NC);

  memcpy((double *)(e1->P2_t),(const double *)(e0->P2_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->Y2_t),(const double *)(e0->Y2_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->Y2s_t),(const double *)(e0->Y2s_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->Kd_t),(const double *)(e0->Kd_t),sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->D_t),(const double *)(e0->D_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->Db_t),(const double *)(e0->Db_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->Dh_t),(const double *)(e0->Dh_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->Dt_t),(const double *)(e0->Dt_t),sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->exrate_t),(const double *)(e0->exrate_t),sizeof(double)*(NT+1)*NC*(NC-1));
  memcpy((double *)(e1->exitrate_t),(const double *)(e0->exitrate_t),sizeof(double)*(NT+1)*NC*(NC-1));


  uint t, i, ii;
  for(t=0; t<(NT+1); t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(ii=0; ii<2; ii++)
	    {
	      copy_exporter_vars(&(e1->ev_t[t][i][ii]),(const exporter_vars *)(&(e0->ev_t[t][i][ii])));
	    }
	}
    }

  memcpy((double *)(e1->welfare_t),(const double *)(e0->welfare_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare2_t),(const double *)(e0->welfare2_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare3_t),(const double *)(e0->welfare3_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare4_t),(const double *)(e0->welfare4_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare_cost_t),(const double *)(e0->welfare_cost_t),sizeof(double)*(NT+1)*NC);


}

uint stack_bgp_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t = NT;
  
  COPY_SUBVECTOR_LOG(myx+nx,e->W_t[t],NC);
  nx=nx+NC;

  COPY_SUBVECTOR_LOG(myx+nx,e->Y_t[t],NC);
  nx=nx+NC;

  COPY_SUBVECTOR_LOG(myx+nx,e->K_t[t],NC);
  nx=nx+NC;

  COPY_SUBVECTOR_LOG(myx+nx,e->P2_t[t],NC*NC);
  nx=nx+NC*NC;

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

  copy_subvector_exp( (double *)(e->W_t[t]), myx+nx, NC);
  nx=nx+NC;

  COPY_SUBVECTOR_EXP(e->Y_t[t],myx+nx,NC);
  nx=nx+NC;

  COPY_SUBVECTOR_EXP(e->K_t[t],myx+nx,NC);
  nx=nx+NC;

  COPY_SUBVECTOR_EXP(e->P2_t[t],myx+nx,NC*NC);
  nx=nx+NC*NC;

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
  else if(scenario==7 || scenario==8 || (scenario>=1000 && scenario<2000))
    {
      nn = NT+1-TCH0;
      t0 = TCH0;
    }
  else if(scenario==10 || scenario==11 || scenario >= 2000)
    {
      nn = NT+1-TCH1;
      t0 = TCH1;
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),NC*(nn-1));
  nx = nx + NC*(nn-1);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),NC*(nn-1));
  nx = nx + NC*(nn-1);

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
  uint t0 =0;
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
  else if(scenario==7 || scenario==8 || (scenario>=1000 && scenario<2000))
    {
      nn = NT+1-TCH0;
      t0 = TCH0;
    }
  else if(scenario==10 || scenario==11 || scenario>=2000)
    {
      nn = NT+1-TCH1;
      t0 = TCH1;
    }

  COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->B_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,NC*(nn-1));
  nx = nx + NC*(nn-1);

  COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,NC*(nn-1));
  nx = nx + NC*(nn-1);

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking eqm vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}

uint stack_stoch_eqm_vars(double * myx, const stoch_eqm * se)
{
  uint i = 0;
  uint t = 0;
  uint nx = 0;
  uint t0 = TREF;
  uint nn = TVOTE-TREF;

  // first stack the deterministic part for the periods between the referendum announcement and the vote
  const eqm * e = &((se->eee)[0]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TVOTE; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),NC*nn);
  nx = nx + NC*nn;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),NC*nn);
  nx = nx + NC*nn;

  // now stack the deterministic part associated with the "no vote" branch
  t0 = TVOTE;
  nn = NT+1-TVOTE;  
  e = &((se->eee)[0]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),NC*(nn-1));
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  // now stack the TBREXIT-TVOTE periods associated with the "yes vote" branch
  t0 = TVOTE;
  nn = TBREXIT-TVOTE;
  e = &((se->eee)[1]);
  
  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TBREXIT; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  // now stack the post-Brexit part for the hard and soft branches
  t0 = TBREXIT;
  nn = NT+1-TBREXIT;

  uint ih;
  for(ih=1; ih<NHIST; ih++)
    {
      e = &((se->eee)[ih]);

      COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      *(myx+nx) = e->B_t[t][i];
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error stacking stoch_eqm vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}

#ifdef CH_SIMPLE
uint stack_stoch_eqm_ch_vars(double * myx, const stoch_eqm_ch * se)
{
  uint i = 0;
  uint t = 0;
  uint nx = 0;
  uint t0 = TCH0;
  uint nn = TCH1-TCH0;

  // first stack the deterministic part for the periods between announcement and implementation
  const eqm * e = &((se->eee)[0]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TCH1; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),NC*nn);
  nx = nx + NC*nn;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),NC*nn);
  nx = nx + NC*nn;

  // now stack the deterministic part associated with the "no implementation" branch
  t0 = TCH1;
  nn = NT+1-TCH1;  
  e = &((se->eee)[0]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),NC*(nn-1));
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  // now stack the TCH2-TCH1 periods associated with the "implementation" branch
  t0 = TCH1;
  nn = TCH2-TCH1;
  e = &((se->eee)[1]);
  
  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TCH2; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  // now stack the TCH2-onward branches for both possibilities
  t0 = TCH2;
  nn = NT+1-TCH2;

  uint ih;
  for(ih=1; ih<NHIST_CH; ih++)
    {
      e = &((se->eee)[ih]);

      COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      *(myx+nx) = e->B_t[t][i];
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error stacking stoch_eqm_ch vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}
#else
uint stack_stoch_eqm_ch_vars(double * myx, const stoch_eqm_ch * se)
{
  int ih = 0;
  uint i = 0;
  uint t = 0;
  uint nx = 0;
  uint t0 = TCH0;
  uint nn = TCH1-TCH0;

  const eqm * e = &((se->eee)[0]);

  // first stack the deterministic part for the periods between announcement and implementation
  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TCH1; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),NC*nn);
  nx = nx + NC*nn;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),NC*nn);
  nx = nx + NC*nn;

  // now stack the deterministic part associated with the "no implementation" branch
  t0 = TCH1;
  nn = NT+1-TCH1;  
  e = &((se->eee)[0]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),NC*(nn-1));
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  // now stack TCH2-TCH1 part where multiple branches coincide...
  // need only one copy of this, which we store in the last branch
  t0 = TCH1;
  nn = TCH2-TCH1;
  e = &((se->eee)[NHIST_CH-1]);

  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TCH2; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->B_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  // now the parts where the other branches are separate
  for(ih=1; ih<NHIST_CH; ih++)
    {
      e = &((se->eee)[ih]);
      t0 = TCH1+ih;
      nn = NT+1-(TCH1+ih);

      if(ih==NHIST_CH-1)
	{
	  t0 = TCH1+ih-1;
	  nn = NT+1-(TCH1+ih-1);
	}
      
      COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      *(myx+nx) = e->B_t[t][i];
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_LOG(myx+nx,&(e->Q_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->R_t[t0]),(nn-1)*NC);
      nx = nx + (nn-1)*NC;
    }
  
  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error stacking stoch_eqm_ch vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}
#endif

uint unstack_stoch_eqm_vars(stoch_eqm * se, const double * myx)
{
  uint i = 0;
  uint t = 0;
  uint nx = 0;
  uint t0 = TREF;
  uint nn = TVOTE-TREF;

  // first stack the deterministic part for the periods between the referendum announcement and the vote
  eqm * e = &((se->eee)[0]);

  // when we unstack we want to copy over to all histories for the deterministic part,
  // since they have to be identical for this period
  int ih;
  for(ih=NHIST-1; ih>=0; ih--)
    {
      e = &((se->eee)[ih]);
      nx = 0;

      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<=TVOTE; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,nn*NC);
      nx = nx + nn*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,nn*NC);
      nx = nx + nn*NC;
    }

  // now stack the deterministic part associated with the "no vote" branch
  t0 = TVOTE;
  nn = NT+1-TVOTE;  
  e = &((se->eee)[0]);

  COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->B_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  // now do the TBREXIT-TVOTE period... here we need to copy to both branches 1 and 2 since they
  // must be identical here
  t0 = TVOTE;
  nn = TBREXIT-TVOTE;
  uint nx2=nx;
  for(ih=(NHIST-1); ih>=1; ih--)
    {
      nx=nx2;
      e = &((se->eee)[ih]);
  
      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<=TBREXIT; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;
    }

  // now stack the post-Brexit part for the hard and soft branches
  t0 = TBREXIT;
  nn = NT+1-TBREXIT;

  for(ih=1; ih<NHIST; ih++)
    {
      e = &((se->eee)[ih]);

      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking stoch_eqm vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}

#ifdef CH_SIMPLE
uint unstack_stoch_eqm_ch_vars(stoch_eqm_ch * se, const double * myx)
{
  uint i = 0;
  uint t = 0;
  uint nx = 0;
  uint t0 = TCH0;
  uint nn = TCH1-TCH0;

  // first stack the deterministic part for the periods between the referendum announcement and the vote
  eqm * e = &((se->eee)[0]);

  // when we unstack we want to copy over to all histories for the deterministic part,
  // since they have to be identical for this period
  int ih;
  for(ih=NHIST_CH-1; ih>=0; ih--)
    {
      e = &((se->eee)[ih]);
      nx = 0;

      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<=TCH1; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,nn*NC);
      nx = nx + nn*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,nn*NC);
      nx = nx + nn*NC;
    }

  // now stack the deterministic part associated with the "no implementation" branch
  t0 = TCH1;
  nn = NT+1-TCH1;  
  e = &((se->eee)[0]);

  COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->B_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  // now do the TCH2-TCH1 period... here we need to copy to both branches 1 and 2 since they
  // must be identical here
  t0 = TCH1;
  nn = TCH2-TCH1;
  uint nx2=nx;
  for(ih=(NHIST_CH-1); ih>=1; ih--)
    {
      nx=nx2;
      e = &((se->eee)[ih]);
  
      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<=TCH2; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;
    }

  // now stack the post-implementation branches
  t0 = TCH2;
  nn = NT+1-TCH2;

  for(ih=1; ih<NHIST_CH; ih++)
    {
      e = &((se->eee)[ih]);

      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking stoch_eqm_ch vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}
#else
uint unstack_stoch_eqm_ch_vars(stoch_eqm_ch * se, const double * myx)
{
  uint i = 0;
  uint t = 0;
  uint nx = 0;
  uint t0 = TCH0;
  uint nn = TCH1-TCH0;
  eqm * e = &((se->eee)[0]);

  /*fprintf(logfile,"\n\1");
  for(t=0; t<TCH2+4; t++)
    {
      fprintf(logfile,KRED "%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n" RESET,t,
	      (se->eee)[0].B_t[t][0],
	      (se->eee)[1].B_t[t][0],
	      (se->eee)[2].B_t[t][0],
	      (se->eee)[3].B_t[t][0],
	      (se->eee)[4].B_t[t][0]);
    }
    fprintf(logfile,"\n");*/

  // first stack the deterministic part for the periods between the referendum announcement and the vote
  //  

  // when we unstack we want to copy over to all histories for the deterministic part,
  // since they have to be identical for this period
  int ih;
  for(ih=NHIST_CH-1; ih>=0; ih--)
    {
      e = &((se->eee)[ih]);
      nx = 0;

      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<=TCH1; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,nn*NC);
      nx = nx + nn*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,nn*NC);
      nx = nx + nn*NC;
    }

  // now stack the deterministic part associated with the "no reform" branch
  t0 = TCH1;
  nn = NT+1-TCH1;  
  e = &((se->eee)[0]);

  COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->B_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn-1)*NC);
  nx = nx + (nn-1)*NC;

  /*fprintf(logfile,"\n\2");
  for(t=0; t<TCH2+4; t++)
    {
      fprintf(logfile,KRED "%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n" RESET,t,
	      (se->eee)[0].B_t[t][0],
	      (se->eee)[1].B_t[t][0],
	      (se->eee)[2].B_t[t][0],
	      (se->eee)[3].B_t[t][0],
	      (se->eee)[4].B_t[t][0]);
    }
    fprintf(logfile,"\n");*/

  // now do the branches with different implementation lengths...
  // first we copy the stuff for the TCH1-TCH2 period stored in the last branch...
  t0 = TCH1;
  nn = TCH2-TCH1;  
  e = &((se->eee)[NHIST_CH-1]);
  
  COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TCH2; t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->B_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  /*  fprintf(logfile,"\n\3");
  for(t=0; t<TCH2+4; t++)
    {
      fprintf(logfile,KRED "%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n" RESET,t,
	      (se->eee)[0].B_t[t][0],
	      (se->eee)[1].B_t[t][0],
	      (se->eee)[2].B_t[t][0],
	      (se->eee)[3].B_t[t][0],
	      (se->eee)[4].B_t[t][0]);
    }
    fprintf(logfile,"\n");*/

  // ...then we copy it to all the other branches
  for(ih=1; ih<NHIST_CH-1; ih++)
    {
      //nn = ih;
      eqm * e2 = &((se->eee)[ih]);
      COPY_SUBVECTOR(e2->W_t[t0],e->W_t[t0],nn*NC);
      COPY_SUBVECTOR(e2->Y_t[t0],e->Y_t[t0],nn*NC);
      COPY_SUBVECTOR(e2->X_t[t0],e->X_t[t0],nn*NC);
      COPY_SUBVECTOR(e2->P2_t[t0],e->P2_t[t0],nn*NC*NC);
      for(t=t0+1; t<=TCH2; t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e2->B_t[t][i] = e->B_t[t][i];
	    }
	}
      COPY_SUBVECTOR(e2->Q_t[t0],e->Q_t[t0],nn*NC);
      COPY_SUBVECTOR(e2->R_t[t0],e->R_t[t0],nn*NC);
    }

  // now stack the post-reversion part for each branch... this overwrites what we did above, but
  // only partially
  for(ih=1; ih<NHIST_CH; ih++)
    {
      t0 = TCH1+ih;
      nn = NT+1-t0;
      if(ih==NHIST_CH-1)
	{
	  t0 = TCH1+ih-1;
	  nn = NT+1-t0;
	}
      e = &((se->eee)[ih]);

      COPY_SUBVECTOR_EXP(&(e->W_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->Y_t[t0]),myx+nx,(nn)*NC);
      nx = nx + (nn)*NC;

      COPY_SUBVECTOR_EXP(&(e->X_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->P2_t[t0]),myx+nx,(nn)*NC*NC);
      nx = nx + (nn)*NC*NC;

      for(t=t0+1; t<(NT+1); t++)
	{
	  for(i=0; i<(NC-1); i++)
	    {
	      e->B_t[t][i] = *(myx+nx);
	      nx = nx+1;
	    }
	}

      COPY_SUBVECTOR_EXP(&(e->Q_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->R_t[t0]),myx+nx,(nn-1)*NC);
      nx = nx + (nn-1)*NC;
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking stoch_eqm_ch vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  /*fprintf(logfile,"\n\4");
  for(t=0; t<TCH2+4; t++)
    {
      fprintf(logfile,KRED "%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n" RESET,t,
	      (se->eee)[0].B_t[t][0],
	      (se->eee)[1].B_t[t][0],
	      (se->eee)[2].B_t[t][0],
	      (se->eee)[3].B_t[t][0],
	      (se->eee)[4].B_t[t][0]);
    }
    fprintf(logfile,"\n");*/

  return 0;
} 
#endif

uint set_initial_bgp_guess()
{
  uint i, t, j;
  eqm * e;
  params * p;

  if(scenario==0)
    {
      e = &(eee0[0]);
      p = &(ppp0[0]);
    }
  else if(scenario==6)
    {
      e = &(eee0_ch[0]);
      p = &(ppp0_ch[0]);
    }
  t=NT;
  
  for(i=0; i<NC; i++)
    {
      e->W_t[t][i] = 1.0;
      e->K_t[t][i] = p->K0[i];
      e->Y_t[t][i] = p->Y0[i];
      
      for(j=0; j<NC; j++)
	{
	  e->P2_t[t][i][j] = 1.0;
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
  uint i,t;
  double bb[NC];
  eqm * e;
  params * p;

  if(scenario==0)
    {
      e = &(eee0[0]);
      p = &(ppp0[0]);
    }
  else if(scenario==6)
    {
      e = &(eee0_ch[0]);
      p = &(ppp0_ch[0]);
    }

  bb[0] = p->B0[0];
  bb[1] = p->B0[1];
  bb[2] = p->B0[2];

  par = 0;
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

#ifdef CH_SIMPLE
  if(scenario==6)
    {
      uint it;
      for(it=0; it<NTH; it++)
	{
	  for(i=0; i<NC; i++)
	    {
	      ppp0_ch[it].K0[i] = e->K_t[NT][i];
	      ppp0_ch[it].n0[i][0] = e->ev_t[NT][i][0].n;
	      ppp0_ch[it].n0[i][1] = e->ev_t[NT][i][1].n;

	      ppp1_ch[it].K0[i] = e->K_t[NT][i];
	      ppp1_ch[it].n0[i][0] = e->ev_t[NT][i][0].n;
	      ppp1_ch[it].n0[i][1] = e->ev_t[NT][i][1].n;

	      ppp2_ch[it].K0[i] = e->K_t[NT][i];
	      ppp2_ch[it].n0[i][0] = e->ev_t[NT][i][0].n;
	      ppp2_ch[it].n0[i][1] = e->ev_t[NT][i][1].n;

	      stoch_params_ch * sp_ch = &(sppp_ch[it]); 
	      uint ih;
	      for(ih=0; ih<NHIST_CH; ih++)
		{
		  (sp_ch->ppp)[ih].K0[i] = e->K_t[NT][i];
		  (sp_ch->ppp)[ih].n0[i][0] = e->ev_t[NT][i][0].n;
		  (sp_ch->ppp)[ih].n0[i][1] = e->ev_t[NT][i][1].n;
		}
	    }
	}
    }
#else
  if(scenario==6)
    {
      uint it;
      for(it=0; it<NTH; it++)
	{
	  for(i=0; i<NC; i++)
	    {
	      ppp0_ch[it].K0[i] = e->K_t[NT][i];
	      ppp0_ch[it].n0[i][0] = e->ev_t[NT][i][0].n;
	      ppp0_ch[it].n0[i][1] = e->ev_t[NT][i][1].n;

	      uint is;
	      for(is=0; is<NHIST_CH-1; is++)
		{
		  ppp1_ch[is][it].K0[i] = e->K_t[NT][i];
		  ppp1_ch[is][it].n0[i][0] = e->ev_t[NT][i][0].n;
		  ppp1_ch[is][it].n0[i][1] = e->ev_t[NT][i][1].n;
		}

	      stoch_params_ch * sp_ch = &(sppp_ch[it]); 
	      uint ih;
	      for(ih=0; ih<NHIST_CH; ih++)
		{
		  (sp_ch->ppp)[ih].K0[i] = e->K_t[NT][i];
		  (sp_ch->ppp)[ih].n0[i][0] = e->ev_t[NT][i][0].n;
		  (sp_ch->ppp)[ih].n0[i][1] = e->ev_t[NT][i][1].n;
		}
	    }
	}
    }
#endif

  // first construct bond guess... a little awkward to logspace this because we have to deal with
  // absolute values
  double tmpb0[NC] = {fabs(p->B0[0]),fabs(p->B0[1]),fabs(p->B0[2])};
  double tmpb1[NC] = {fabs(bb[0]),fabs(bb[1]),fabs(bb[2])};
  LOGSPACE_2D(tmpb0,tmpb1,NT+1,NC,e->B_t);
  for(i=0; i<(NC-1); i++)
    {
      if(fabs(p->B0[i])<1.0e-6)
	{
	  for(t=0; t<(NT+1); t++)
	    {
	      e->B_t[t][i] = 0.0;
	    }
	}
      else
	{
	  if(p->B0[i] < -TINY)
	    {
	      for(t=0; t<(NT+1); t++)
		{
		  e->B_t[t][i] = -e->B_t[t][i];
		}
	    }
	}
    }

  // now construct guesses for prices real variables
  if(scenario==0)
    {
      double tmpp[NC] = {1.0,1.0,1.0};
      double tmpp2[NC][NC] = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
      LOGSPACE_2D(p->K0,e->K_t[NT],NT+1,NC,e->K_t);
      LOGSPACE_2D(p->Y0,e->Y_t[NT],NT+1,NC,e->Y_t);
      LINSPACE_2D(tmpp,e->W_t[NT],NT+1,NC,e->W_t);
      LINSPACE_2D(p->R0,e->R_t[NT],NT+1,NC,e->R_t);
      LINSPACE_2D(tmpp2,e->P2_t[NT],NT+1,NC*NC,e->P2_t);
      SET_ALL_V(e->Q_t,(NT+1)*NC,p->beta);
    }
  else if(scenario==6)
    {
      LOGSPACE_2D(e->K_t[NT],e->K_t[NT],NT+1,NC,e->K_t);
      LOGSPACE_2D(e->Y_t[NT],e->Y_t[NT],NT+1,NC,e->Y_t);
      LINSPACE_2D(e->W_t[NT],e->W_t[NT],NT+1,NC,e->W_t);
      LINSPACE_2D(e->R_t[NT],e->R_t[NT],NT+1,NC,e->R_t);
      LINSPACE_2D(e->P2_t[NT],e->P2_t[NT],NT+1,NC*NC,e->P2_t);
      SET_ALL_V(e->Q_t,(NT+1)*NC,p->beta);
    }

  for(i=0; i<NC; i++)
    {
      for(t=0; t<NT; t++)
	{
	  e->X_t[t][i] = e->K_t[t+1][i] - (1.0-p->delta)*e->K_t[t][i];
	}
    }

  if(stack_eqm_vars(solver_x->data,e))
    {
      fprintf(logfile,KRED "Failed to create guess for equilibrium!\n" RESET);
      return 1;
    }
  else
    {      
      return 0;
    }
}

uint write_bgp_vars(const eqm * e, const char * fname)
{
  uint t = NT;

  char * fname2 = concat("output/",fname);
  char * fname3;

  if(nokappa)
    {
      fname3 = concat(fname2,"_nokappa.txt");
    }
  else if(eqkappa)
    {
      fname3 = concat(fname2,"_eqkappa.txt");
    }
  else if(low_pi)
    {
      fname3 = concat(fname2,"_lowpi.txt");
    }
  else if(high_pi)
    {
      fname3 = concat(fname2,"_highpi.txt");
    }
  else if(fin_aut)
    {
      fname3 = concat(fname2,"_finaut.txt");
    }
  else if(comp_mkts)
    {
      fname3 = concat(fname2,"_compmkts.txt");
    }
  else if(zeta_sens)
    {
      fname3 = concat(fname2,"_szeta.txt");
    }
  else if(psi_sens)
    {
      fname3 = concat(fname2,"_spsi.txt");
    }
  else
    {
      fname3 = concat(fname2,"_baseline.txt");
    }

  FILE * file = fopen(fname3,"wb");
  
  free(fname2);
  free(fname3);
  
  if(file)
    {
      fprintf(file,"bonds (exog) : %0.3f\t%0.3f\t%0.3f\n",e->B_t[t][0],e->B_t[t][1],e->B_t[t][2]);
  
      fprintf(file,"\nAggregates   : UK\t\tEU\t\tRW\n");
      fprintf(file,"capital      : %0.3f\t\t%0.3f\t%0.3f\n",e->K_t[t][0],e->K_t[t][1],e->K_t[t][2]);
      fprintf(file,"labor        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->L_t[t][0],e->L_t[t][1],e->L_t[t][2]);
      fprintf(file,"intermediates: %0.3f\t\t%0.3f\t\t%0.3f\n",e->M_t[t][0],e->M_t[t][1],e->M_t[t][2]);
      fprintf(file,"gross output : %0.3f\t\t%0.3f\t\t%0.3f\n",e->Y_t[t][0],e->Y_t[t][1],e->Y_t[t][2]);
      fprintf(file,"gdp          : %0.3f\t\t%0.3f\t\t%0.3f\n",e->RGDP_t[t][0],e->RGDP_t[t][1],e->RGDP_t[t][2]);
      fprintf(file,"investment   : %0.3f\t\t%0.3f\t\t%0.3f\n",e->X_t[t][0],e->X_t[t][1],e->X_t[t][2]);
      fprintf(file,"trd lab      : %0.3f\t\t%0.3f\t\t%0.3f\n",e->Lf_t[t][0],e->Lf_t[t][1],e->Lf_t[t][2]);
      fprintf(file,"consumption  : %0.3f\t\t%0.3f\t\t%0.3f\n",e->C_t[t][0],e->C_t[t][1],e->C_t[t][2]);
      fprintf(file,"wages        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->W_t[t][0],e->W_t[t][1],e->W_t[t][2]);
      fprintf(file,"cpi          : %0.3f\t\t%0.3f\t\t%0.3f\n",e->P_t[t][0],e->P_t[t][1],e->P_t[t][2]);

      fprintf(file,"\nBilateral flows\n");
      uint i=0;
      fprintf(file,"UK Y2        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->Y2_t[t][i][0],e->Y2_t[t][i][1],e->Y2_t[t][i][2]);
      i=1;
      fprintf(file,"EU Y2        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->Y2_t[t][i][0],e->Y2_t[t][i][1],e->Y2_t[t][i][2]);
      i=2;
      fprintf(file,"RW Y2        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->Y2_t[t][i][0],e->Y2_t[t][i][1],e->Y2_t[t][i][2]);

      i=0;
      fprintf(file,"UK P2        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->P2_t[t][i][0],e->P2_t[t][i][1],e->P2_t[t][i][2]);
      i=1;
      fprintf(file,"EU P2        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->P2_t[t][i][0],e->P2_t[t][i][1],e->P2_t[t][i][2]);
      i=2;
      fprintf(file,"RW P2        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->P2_t[t][i][0],e->P2_t[t][i][1],e->P2_t[t][i][2]);

      i=0;
      fprintf(file,"UK NX        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->NX_t[t][i][0],e->NX_t[t][i][1],e->NX_t[t][i][2]);
      i=1;
      fprintf(file,"EU NX        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->NX_t[t][i][0],e->NX_t[t][i][1],e->NX_t[t][i][2]);
      i=2;
      fprintf(file,"RW NX        : %0.3f\t\t%0.3f\t\t%0.3f\n",e->NX_t[t][i][0],e->NX_t[t][i][1],e->NX_t[t][i][2]);

      fprintf(file,"\nExporter dynamics\n");
      i=0;
      fprintf(file,"UK exp. rate : %0.3f\t\t%0.3f\n",100*e->exrate_t[t][i][0],100*e->exrate_t[t][i][1]);
      i=1;
      fprintf(file,"EU exp. rate : %0.3f\t\t%0.3f\n",100*e->exrate_t[t][i][0],100*e->exrate_t[t][i][1]);
      i=2;
      fprintf(file,"RW exp. rate : %0.3f\t\t%0.3f\n",100*e->exrate_t[t][i][0],100*e->exrate_t[t][i][1]);

      i=0;
      fprintf(file,"UK exit rate : %0.3f\t\t%0.3f\n",100*e->exitrate_t[t][i][0],100*e->exitrate_t[t][i][1]);
      i=1;
      fprintf(file,"EU exit rate : %0.3f\t\t%0.3f\n",100*e->exitrate_t[t][i][0],100*e->exitrate_t[t][i][1]);
      i=2;
      fprintf(file,"RW exit rate : %0.3f\t\t%0.3f\n",100*e->exitrate_t[t][i][0],100*e->exitrate_t[t][i][1]);

      fclose(file);

      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write balanced growth path!\n" RESET);
      return 1;
    }
}

uint write_eqm_vars(const params * p, const eqm * e, const char * fname, uint i)
{
  char * fname2 = concat("output/",fname);
  char * fname3;

  if(nokappa)
    {
      fname3 = concat(fname2,"_nokappa.csv");
    }
  else if(eqkappa)
    {
      fname3 = concat(fname2,"_eqkappa.csv");
    }
  else if(low_pi)
    {
      fname3 = concat(fname2,"_lowpi.csv");
    }
  else if(high_pi)
    {
      fname3 = concat(fname2,"_highpi.csv");
    }
  else if(fin_aut)
    {
      fname3 = concat(fname2,"_finaut.csv");
    }
  else if(comp_mkts)
    {
      fname3 = concat(fname2,"_compmkts.csv");
    }
  else if(zeta_sens)
    {
      fname3 = concat(fname2,"_szeta.csv");
    }
  else if(psi_sens)
    {
      fname3 = concat(fname2,"_spsi.csv");
    }
  else if(scenario>=6 && (int)TCH2>(int)TBREXIT)
    {
      fname3 = concat(fname2,"_long.csv");
    }
  else
    {
      fname3 = concat(fname2,"_baseline.csv");
    }

  FILE * file = fopen(fname3,"wb");
  
  free(fname2);
  free(fname3);
  
  uint ii, j, t;
  if(file)
    {
      fprintf(file,"period,rgdp,ngdp,y,x,lf,c,W,sW,cW,ceW1,ceW2");
      for(ii=0; ii<(NC-1); ii++)
	{
	  if(i==0)
	    {
	      if(ii==0)
		{
		  j=1;
		}
	      else
		{
		  j=2;
		}
	    }
	  else if(i==1)
	    {
	      if(ii==0)
		{
		  j=0;
		}
	      else
		{
		  j=2;
		}
	    }
	  else if(i==2)
	    {
	      if(ii==0)
		{
		  j=0;
		}
	      else
		{
		  j=1;
		}
	    }
	  fprintf(file,",rer%d,tot%d,nx%d,rex%d,rim%d,exrate%d,exitrate%d,tau%d,ntb%d",j,j,j,j,j,j,j,j,j);
	}
      fprintf(file,"\n");

      for(t=0;t<(NT+1);t++)
	{
	  fprintf(file,"%d,",t);
	  fprintf(file,"%0.16f,",e->RGDP_t[t][i]);
	  fprintf(file,"%0.16f,",e->NGDP_t[t][i]);
	  fprintf(file,"%0.16f,",e->Y_t[t][i]);
	  fprintf(file,"%0.16f,",e->X_t[t][i]);
	  fprintf(file,"%0.16f,",e->Lf_t[t][i]);
	  fprintf(file,"%0.16f,",e->C_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare2_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare_cost_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare3_t[t][i]);
	  fprintf(file,"%0.16f",e->welfare4_t[t][i]);
	  for(ii=0; ii<(NC-1); ii++)
	    {
	      if(i==0)
		{
		  if(ii==0)
		    {
		      j=1;
		    }
		  else
		    {
		      j=2;
		    }
		}
	      else if(i==1)
		{
		  if(ii==0)
		    {
		      j=0;
		    }
		  else
		    {
		      j=2;
		    }
		}
	      else if(i==2)
		{
		  if(ii==0)
		    {
		      j=0;
		    }
		  else
		    {
		      j=1;
		    }
		}
	      fprintf(file,",%0.16f,",e->RER_t[t][i][j]);
	      fprintf(file,"%0.16f,",e->TOT_t[t][i][j]);
	      fprintf(file,"%0.16f,",e->NX_t[t][i][j]);
	      fprintf(file,"%0.16f,",e->EXR_t[t][i][j]);
	      fprintf(file,"%0.16f,",e->IMR_t[t][i][j]);
	      fprintf(file,"%0.16f,",e->exrate_t[t][i][ii]);
	      fprintf(file,"%0.16f,",e->exitrate_t[t][i][ii]);
	      fprintf(file,"%0.16f,",p->tau_ts[t][i][j]);
	      fprintf(file,"%0.16f",p->ntb_ts[t][i][j]);
	    }
	  fprintf(file,"\n");
	}

      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write equilibrium vars!\n" RESET);
      return 1;
    }
}

uint set_vars1(eqm * e, const params * p, uint t, uint bgp)
{
  uint i,j;

  // initialize some things to zero
  SET_ALL_V(e->EX_t[t],NC*NC,0.0);
  SET_ALL_V(e->IM_t[t],NC*NC,0.0);
  SET_ALL_V(e->NX_t[t],NC*NC,0.0);
  SET_ALL_V(e->EXR_t[t],NC*NC,0.0);
  SET_ALL_V(e->IMR_t[t],NC*NC,0.0);

  // bond market clearing implies the third bond
  e->B_t[t][2] = -(e->B_t[t][0]+e->B_t[t][1]);

  // first set the aggregate price level based on the bilateral prices
  for(i=0; i<NC; i++)
    {
      e->P_t[t][i] = 0.0;
      for(j=0; j<NC; j++)
	{
	  e->P_t[t][i] += p->mu[i][j] * pow((1.0+p->tau_ts[t][i][j])*e->P2_t[t][i][j],1.0-p->zeta);
	}
      e->P_t[t][i] = pow(e->P_t[t][i],1.0/(1.0-p->zeta)) / p->Ybar[i];
    }

  // now compute bilateral aggregate quantities using aggregator's FOC
  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  e->Y2_t[t][i][j] = pow(
				 (1.0+p->tau_ts[t][i][j])*e->P2_t[t][i][j] 
				 / (p->mue[i][j] * p->Ybare[i] * e->P_t[t][i]),
				 -p->zeta) * e->Y_t[t][i];

	  if(i!=j)
	    {
	      e->IM_t[t][i][j] = e->P2_t[t][i][j]*e->Y2_t[t][i][j];
	      e->IMR_t[t][i][j] = e->Y2_t[t][i][j];
	    }
	  else
	    {
	      e->IM_t[t][i][j] = 0.0;
	    }
	}

      double test = arm_agg(e->Y2_t[t][i],p->Ybar[i],p->mue[i],p->fzeta,p->ifzeta);
      if(fabs(test-e->Y_t[t][i])>1.0e-6)
	{
	  fprintf(logfile,KRED "Error! Y != ArmAgg(Y2), country=%d!\n" RESET,i);
	  return 1;
	}

      
    }
  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  //e->EX_t[t][i][j] = (1.0+p->ntb_ts[t][j][i])*e->IM_t[t][j][i];
	  e->EX_t[t][i][j] = e->IM_t[t][j][i];
	  e->EXR_t[t][i][j] = e->EX_t[t][i][j]/e->P2_t[t][j][i];
	  e->NX_t[t][i][j] = e->EX_t[t][i][j] - e->IM_t[t][i][j];
	  e->RER_t[t][i][j] = e->P_t[t][j] / e->P_t[t][i];
	  e->TOT_t[t][i][j] = (1.0+p->tau_ts[t][j][i])*e->P2_t[t][j][i] / e->P2_t[t][i][j] / (1.0+p->tau_ts[t][i][j]);
	}
    }

  // compute the investment stuff...
  for(i=0; i<NC; i++)
    {
      if(t<NT)
	{
	  if(t==(NT-1) || p->cap_adj_cost==0)
	    {
	      e->K_t[t+1][i] = (1.0-p->delta) * e->K_t[t][i] + e->X_t[t][i];
	    }
	  else
	    {
	      e->K_t[t+1][i] = (1.0-p->delta) * e->K_t[t][i] + 
		phiK(e->X_t[t][i]/e->K_t[t][i],p->delta,p->etaK) * e->K_t[t][i];
	    }
	}
      else
	{
	  e->X_t[t][i] = p->delta * e->K_t[t][i];
	}
	
      if(t==NT)
	{
	  e->Q_t[t][i] = e->P_t[t][0]*p->beta;

	  if(bgp)
	    {
	      e->R_t[t][i] = e->P_t[t][i] * e->P_t[t][0]/e->Q_t[t][i] 
		- e->P_t[t][i]*(1.0-p->delta);
	    }
	  else
	    {
	      e->R_t[t][i] = e->P_t[t-1][i] * e->P_t[t][0]/e->Q_t[t-1][i] 
		- e->P_t[t][i]*(1.0-p->delta);	      
	    }
	  e->Lambda_t[t][i] = p->beta;
	}
      else
	{
	  // note that we can do this despite the presence of uncertainty, because the firm value function
	  // updating uses separate branches already
	  e->Lambda_t[t][i] = e->Q_t[t][i];
	}
    }  

  // now compute the constants we need to compute firm-level stuff
  for(i=0; i<NC; i++)
    {
      if(leontief)
	{
	  e->MC_t[t][i] = p->gamma[i] *
	    (pow(e->R_t[t][i]/p->alpha,p->alpha) *
	     pow(e->W_t[t][i]/(1.0-p->alpha),1.0-p->alpha)
	     )
	    + (1.0-p->gamma[i])*e->P_t[t][i];
	}
      else
	{
	  e->MC_t[t][i] = 
	    pow(e->R_t[t][i]/(p->alpha*p->gamma[i]),p->gamma[i]*p->alpha) * 
	    pow(e->W_t[t][i]/(p->gamma[i]*(1.0-p->alpha)),p->gamma[i]*(1.0-p->alpha)) * 
	    pow(e->P_t[t][i]/(1.0-p->gamma[i]),1.0-p->gamma[i]);
	}

      for(j=0; j<NC; j++)
	{
	  e->D_t[t][j][i] = pow(p->Ybar2[j][i], p->theta-1.0) * pow(e->P2_t[t][j][i],p->theta) * e->Y2_t[t][j][i];

	  e->Db_t[t][j][i] = pow(p->theta/(p->theta-1.0),1.0-p->theta) * 
	    pow(1.0+p->ntb_ts[t][j][i],1.0-p->theta) * e->D_t[t][j][i] * pow(e->MC_t[t][i],1.0-p->theta);

	  e->Dh_t[t][j][i] = pow(p->theta/(p->theta-1.0),-p->theta) * 
	  pow(1.0+p->ntb_ts[t][j][i],1.0-p->theta) * e->D_t[t][j][i] * pow(e->MC_t[t][i],1.0-p->theta);

	  e->Dt_t[t][j][i] = e->Db_t[t][j][i]/p->theta;
	} 
    }

  return 0;
}

uint set_vars2(eqm * e, const params * p, uint t, uint bgp)
{
  uint i, ii, j;

  // steady state exporter dynamics
  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  if(t==NT)
	    {
	      if(ev_steady_state(p->sigma[i],
				 p->theta,
				 p->kappa0[i][ii],
				 p->kappa1[i][ii],
				 e->W_t[t][i],
				 e->Lambda_t[t][i],
				 e->Dt_t[t][j][i],
				 p->n0[i][ii],
				 &(e->ev_t[t][i][ii])))
		{
		  fprintf(logfile,
			  KRED "Failed to solve for steady state exporter moments for country/dest %d/%d!\n" RESET,
			  i,j);
		  return 1;
		}
	    }
	  else
	    {
	      calc_z_Z_Fz(p->sigma[i],p->theta,p->kappa0[i][ii],p->kappa1[i][ii],
			  e->W_t[t][i],e->Dt_t[t][j][i],e->ev_t[t+1][i][ii].dV,
			  &(e->ev_t[t][i][ii]));
	      
	      double lam=0.0;
	      if(t==0)
		{
		  lam = p->beta;
		}
	      else
		{
		  lam = e->Lambda_t[t-1][i];
		}
	      update_dV(p->kappa0[i][ii],p->kappa1[i][ii],e->W_t[t][i],lam,
			e->Dt_t[t][j][i],e->ev_t[t+1][i][ii].dV,&(e->ev_t[t][i][ii]));
	    }
	}
   
    }

  return 0;
}

uint set_stoch_vars2(eqm * e, const eqm * e_good, const eqm * e_bad, const params * p, uint t, double pi)
{
  uint i, ii, j;
  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  double W_good=0.0;
	  double W_bad=0.0;
	  double Lambda_good = 0.0;
	  double Lambda_bad = 0.0;
	  double Dtj_good=0.0;
	  double Dtj_bad=0.0;
	  double dVp_good = 0.0;
	  double dVp_bad = 0.0;

	  W_good = e_good->W_t[t][i];
	  W_bad = e_bad->W_t[t][i];

	  Lambda_good = e_good->Lambda_t[t-1][i];
	  Lambda_bad = e_bad->Lambda_t[t-1][i];

	  Dtj_good = e_good->Dt_t[t][j][i];
	  Dtj_bad = e_bad->Dt_t[t][j][i];

	  dVp_good = e_good->ev_t[t+1][i][ii].dV;
	  dVp_bad = e_bad->ev_t[t+1][i][ii].dV;

	  calc_stoch_z_Z_Fz(p->sigma[i],p->theta,p->kappa0[i][ii],p->kappa1[i][ii],
			    W_good,W_bad,Dtj_good,Dtj_bad,dVp_good,dVp_bad,pi,
			    &(e->ev_t[t][i][ii]));

	  update_stoch_dV(p->kappa0[i][ii],p->kappa1[i][ii],
			  W_good,W_bad,Lambda_good,Lambda_bad,
			  Dtj_good,Dtj_bad,dVp_good,dVp_bad,pi,
			  &(e->ev_t[t][i][ii]));
	}
    }

  return 0;

}


void set_vars3(eqm * e, const params * p, uint t, uint bgp)
{

  uint i,j,ii;
  for(i=0; i<NC; i++)
    {      
      if(leontief)
	{
	  e->Kd_t[t][i] = (p->gamma[i]/e->MC_t[t][i]) *
	    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha-1.0) *
	    p->Nbar[i] * e->Dh_t[t][i][i] * e->ev_t[t][i][0].Zi;
	      
	  e->L_t[t][i] = (p->gamma[i]/e->MC_t[t][i]) *
	    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha) *
	    p->Nbar[i] * e->Dh_t[t][i][i] * e->ev_t[t][i][0].Zi;
	      
	  e->M_t[t][i] = ((1.0-p->gamma[i])/e->MC_t[t][i]) * 
	    p->Nbar[i] * e->Dh_t[t][i][i] * e->ev_t[t][i][0].Zi;
	}
      else
	{
	  e->Kd_t[t][i] = (p->gamma[i]*p->alpha/e->R_t[t][i]) * p->Nbar[i] * e->Dh_t[t][i][i] * e->ev_t[t][i][0].Zi;
	  e->L_t[t][i] = (p->gamma[i]*(1.0-p->alpha)/e->W_t[t][i]) * p->Nbar[i] * e->Dh_t[t][i][i] * e->ev_t[t][i][0].Zi;
	  e->M_t[t][i] = ((1.0-p->gamma[i])/e->P_t[t][i]) * p->Nbar[i] * e->Dh_t[t][i][i] * e->ev_t[t][i][0].Zi;
	}

      e->Y2s_t[t][i][i] =  p->Ybar2[i][i] * pow(e->MC_t[t][i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	pow(p->Nbar[i],p->theta/(p->theta-1.0)) * e->D_t[t][i][i] * pow(e->ev_t[t][i][0].Zi,p->theta/(p->theta-1.0));

      e->VA_t[t][i] = e->Db_t[t][i][i] * p->Nbar[i] * e->ev_t[t][i][0].Zi;
      e->Pi_t[t][i] = e->Dt_t[t][i][i] * p->Nbar[i] * e->ev_t[t][i][0].Zi;
     
      e->Lf_t[t][i] = 0.0;

      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  double n0;
	  if(bgp)
	    {
	      n0 = e->ev_t[t][i][ii].n;
	    }
	  else if(t==0)
	    {
	      n0 = p->n0[i][ii];
	    }
	  else
	    {
	      n0 = e->ev_t[t-1][i][ii].n;
	    }

	  update_n(n0,&(e->ev_t[t][i][ii]));
	  e->exrate_t[t][i][ii] = e->ev_t[t][i][ii].n;
	  e->exitrate_t[t][i][ii] = e->ev_t[t][i][ii].Fzm;

	  double om = p->omega[i];
	  if(ii>0)
	    {
	      om = 1.0-p->omega[i];
	    }

	  if(leontief)
	    {
	      e->Kd_t[t][i] += (p->gamma[i]/e->MC_t[t][i]) *
		pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha-1.0) *
		p->Nbar[i] * om * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;
	      
	      e->L_t[t][i] += (p->gamma[i]/e->MC_t[t][i]) *
		pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha) *
		p->Nbar[i] * om * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;
	      
	      e->M_t[t][i] += ((1.0-p->gamma[i])/e->MC_t[t][i]) * 
		p->Nbar[i] * om * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;
	    }
	  else
	    {
	      e->Kd_t[t][i] += (p->gamma[i]*p->alpha/e->R_t[t][i]) * 
		p->Nbar[i] * om * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;

	      e->L_t[t][i] += (p->gamma[i]*(1.0-p->alpha)/e->W_t[t][i]) * 
		p->Nbar[i] * om * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;

	      e->M_t[t][i] += ((1.0-p->gamma[i])/e->P_t[t][i]) * p->Nbar[i] * 
		e->Dh_t[t][j][i] * om * e->ev_t[t][i][ii].Z;
	    }

	  if(!nokappa)
	    {
	      e->Lf_t[t][i] += p->Nbar[i] * om * (
						  p->kappa0[i][ii] * (1.0-n0) * (1.0-e->ev_t[t][i][ii].Fzp) + 
						  p->kappa1[i][ii] * n0 * (1.0-e->ev_t[t][i][ii].Fzm)
						  );
	    }

	  e->Y2s_t[t][j][i] =  p->Ybar2[j][i] * pow(e->MC_t[t][i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	    pow(p->Nbar[i]*om,p->theta/(p->theta-1.0)) * pow(1.0+p->ntb_ts[t][j][i],-p->theta) * e->D_t[t][j][i] * 
	      pow(e->ev_t[t][i][ii].Z,p->theta/(p->theta-1.0));

	  e->Pi_t[t][i] += e->Dt_t[t][j][i] * p->Nbar[i] * om * e->ev_t[t][i][0].Z;
	  e->T_t[t][i] = p->tau_ts[t][i][j] * e->IM_t[t][i][j];
	  e->VA_t[t][i] += e->EX_t[t][i][j];
	}

      e->L_t[t][i] += e->Lf_t[t][i];
      e->NGDP_t[t][i] = e->P_t[t][i]*(e->Y_t[t][i] - e->M_t[t][i]);
      e->RGDP_t[t][i] = e->Y_t[t][i] - e->M_t[t][i];
      e->XY_t[t][i] = 100.0*e->P_t[t][i]*e->X_t[t][i]/e->NGDP_t[t][i];
      e->C_t[t][i] = e->Y_t[t][i] - e->X_t[t][i] - e->M_t[t][i];

      e->MUC_t[t][i] = muc(e->C_t[t][i],
			   e->L_t[t][i],
			   p->Lbar[i],
			   p->phi[i],
			   p->psi);

      e->VA_t[t][i] -= e->P_t[t][i]*e->M_t[t][i];

      e->Inc_t[t][i] = e->W_t[t][i]*e->L_t[t][i] + e->R_t[t][i]*e->K_t[t][i]
	+ e->Pi_t[t][i] + e->T_t[t][i];     

      e->Exp_t[t][i] = e->P_t[t][i]*(e->C_t[t][i] + e->X_t[t][i])
	+ sum(e->NX_t[t][i],NC);
    }
}

uint eval_bgp_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,j=0,t=NT,nx=0;
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);
  if(scenario==6)
    {
      e = &(eee0_ch[tn]);
      p = &(ppp0_ch[tn]);
    }

  e->B_t[t][0] = bbgp[0];
  e->B_t[t][1] = bbgp[1];
  unstack_bgp_vars(e,myx);
  if(set_vars1(e,p,t,1))
    {
      fprintf(logfile,KRED "Error calling set_vars1!\n" RESET);
      return 1;
    }
  if(set_vars2(e,p,t,1))
    {
      fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
      return 1;
    }
  set_vars3(e,p,t,1);
  nx=0;

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  for(i=0; i<NC; i++)
    {
      if(i!=NC-1)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      myf[nx] = mkt_clear_K(e,t,i);
      nx=nx+1;

      for(j=0; j<NC; j++)
	{	  
	  myf[nx] = mkt_clear_Y2(e,t,i,j);
	  nx=nx+1;
	}
    }

  if(nx != nbgp)
    {
      fprintf(logfile,KRED "Error evaluating bgp eqns! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }
  
  return 0;
}

uint eval_eqm_conds(const double * myx, double * myf, uint tn)
{
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);
  uint i=0,j=0,nx=0;
  int t = 0;
  int t0 = 0;

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
  else if(scenario == 6)
    {
      t0 = 0;
      e = &(eee0_ch[tn]);
      p = &(ppp0_ch[tn]);
    }  
#ifdef CH_SIMPLE
  else if(scenario == 7)
    {
      t0 = TCH0;
      e = &(eee1_ch[tn]);
      p = &(ppp1_ch[tn]);
    }  
  else if(scenario == 8)
    {
      t0 = TCH0;
      e = &(eee2_ch[tn]);
      p = &(ppp2_ch[tn]);
    }  
  else if(scenario == 10)
    {
      t0 = TCH1;
      e = &(eee1_ch[tn]);
      p = &(ppp1_ch[tn]);
    }  
  else if(scenario == 11)
    {
      t0 = TCH1;
      e = &(eee2_ch[tn]);
      p = &(ppp2_ch[tn]);
    }  
#else
  else if(scenario>=1000 && scenario <2000)
    {
      t0 = TCH0;
      e = &(eee1_ch[scenario-1000][tn]);
      p = &(ppp1_ch[scenario-1000][tn]);
    }
  else if(scenario>=2000)
    {
      t0 = TCH1;
      e = &(eee1_ch[scenario-2000][tn]);
      p = &(ppp1_ch[scenario-2000][tn]);
    }
#endif

  unstack_eqm_vars(e,myx);

  e->B_t[0][0] = p->B0[0];
  e->B_t[0][1] = p->B0[1];
  e->B_t[0][2] = p->B0[2];

  if(scenario==0 || scenario==6)
    {
      for(i=0; i<NC; i++)
	{
	  e->K_t[0][i] = p->K0[i];
	}
    }

  for(t=t0; t<(NT+1); t++)
    {
      if(set_vars1(e,p,t,0))
	{
	  fprintf(logfile,KRED "Error calling set_vars1!\n" RESET);
	  return 1;
	}
    }

  for(t=NT; t>=t0; t--)
    {
      if(set_vars2(e,p,t,0))
	{
	  fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
	  return 1;
	}
    }
  for(t=t0; t<(NT+1); t++)
    {
      set_vars3(e,p,t,0);
    }

  nx=0;
  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      if(t<(NT-1))
		{
		  myf[nx] = noarb(p,e,t,i);
		  nx=nx+1;
		}

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
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

uint eval_stoch_eqm_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,j=0,t=NT,nx=0,ih=0;

  stoch_eqm * se = &(sss[tn]);
  stoch_params * sp = &(sppp[tn]);
  eqm *e, *e0, *e2;
  params *p, *p2;

  unstack_stoch_eqm_vars(se,myx);

  uint t0 = TREF;
  for(ih=0; ih<NHIST; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(t=t0; t<(NT+1); t++)
	{
	  if(set_vars1(e,p,t,0))
	    {
	      fprintf(logfile,KRED "Error calling set_vars1!\n" RESET);
	      return 1;
	    }
	}
    }

  for(t=NT; t>=t0; t--)
    {
      for(ih=0; ih<NHIST; ih++)
	{
	  e = &( (se->eee)[ih] );
	  p = &( (sp->ppp)[ih] );

	  if(t==TVOTE)
	    {
	      e0 = &(se->eee[0]);
	      e2 = &(se->eee[1]);
	      if(set_stoch_vars2(e,e0,e2,p,t,pi_vote))
		{
		  return 1;
		}
	    }
	  else if(e != &(se->eee[0]) && t==TBREXIT)
	    {
	      e0 = &(se->eee[1]);
	      e2 = &(se->eee[2]);
	      if(set_stoch_vars2(e,e0,e2,p,t,pi_brexit))
		{
		  return 1;
		}
	    }
	  else
	    {
	      if(set_vars2(e,p,t,0))
		{
		  fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
		  return 1;
		}
	    }
	}
    }

  for(ih=0; ih<NHIST; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(t=t0; t<(NT+1); t++)
	{
	  set_vars3(e,p,t,0);
	}
    }

  nx = 0;
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

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      myf[nx] = noarb(p,e,t,i);
	      nx=nx+1;

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
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

  for(i=0; i<NC; i++)
    {
      if(i!=NC-1)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      myf[nx] = mkt_clear_K(e,t,i);
      nx=nx+1;

      myf[nx] = pi_vote * euler(p,e,t,i) + (1.0-pi_vote) * euler(p2,e2,t,i);
      nx = nx+1;

      myf[nx] = pi_vote * noarb(p,e,t,i) + (1.0-pi_vote) * noarb(p2,e2,t,i);
      nx = nx+1;

      if(i>0)
	{
	  myf[nx] = fin_mkt(e,t,i);
	  nx=nx+1;
	}

      for(j=0; j<NC; j++)
	{	  
	  myf[nx] = mkt_clear_Y2(e,t,i,j);
	  nx=nx+1;
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

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      if(t<(NT-1))
		{
		  myf[nx] = noarb(p,e,t,i);
		  nx=nx+1;
		}

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
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

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      myf[nx] = noarb(p,e,t,i);
	      nx=nx+1;

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
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

  for(i=0; i<NC; i++)
    {
      if(i!=NC-1)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      myf[nx] = mkt_clear_K(e,t,i);
      nx=nx+1;

      myf[nx] = pi_brexit * euler(p,e,t,i) + (1.0-pi_brexit) * euler(p2,e2,t,i);
      nx = nx+1;

      myf[nx] = pi_brexit * noarb(p,e,t,i) + (1.0-pi_brexit) * noarb(p2,e2,t,i);
      nx = nx+1;

      if(i>0)
	{
	  myf[nx] = fin_mkt(e,t,i);
	  nx=nx+1;
	}

      for(j=0; j<NC; j++)
	{	  
	  myf[nx] = mkt_clear_Y2(e,t,i,j);
	  nx=nx+1;
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

	  for(i=0; i<NC; i++)
	    {
	      if(i!=NC-1)
		{
		  myf[nx] = bop(p,e,t,i);
		  nx = nx+1;
		}

	      myf[nx] = muc_mul(p,e,t,i);
	      nx=nx+1;

	      myf[nx] = mkt_clear_K(e,t,i);
	      nx=nx+1;

	      if(t<NT)
		{
		  myf[nx] = euler(p,e,t,i);
		  nx = nx+1;

		  if(t<(NT-1))
		    {
		      myf[nx] = noarb(p,e,t,i);
		      nx=nx+1;
		    }

		  if(i>0)
		    {
		      myf[nx] = fin_mkt(e,t,i);
		      nx=nx+1;
		    }
		}

	      for(j=0; j<NC; j++)
		{	  
		  myf[nx] = mkt_clear_Y2(e,t,i,j);
		  nx=nx+1;
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

#ifdef CH_SIMPLE
uint eval_stoch_eqm_ch_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,j=0,t=NT,nx=0,ih=0;

  stoch_eqm_ch * se = &(sss_ch[tn]);
  stoch_params_ch * sp = &(sppp_ch[tn]);
  eqm *e, *e0, *e2;
  params *p, *p2;

  unstack_stoch_eqm_ch_vars(se,myx);

  uint t0 = TCH0;
  for(ih=0; ih<NHIST_CH; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(t=t0; t<(NT+1); t++)
	{
	  if(set_vars1(e,p,t,0))
	    {
	      fprintf(logfile,KRED "Error calling set_vars1!\n" RESET);
	      return 1;
	    }
	}
    }

  for(t=NT; t>=t0; t--)
    {
      for(ih=0; ih<NHIST_CH; ih++)
	{
	  e = &( (se->eee)[ih] );
	  p = &( (sp->ppp)[ih] );

	  if(t==TCH1)
	    {
	      e0 = &(se->eee[0]);
	      e2 = &(se->eee[1]);
	      if(set_stoch_vars2(e,e0,e2,p,t,pi_vote))
		{
		  return 1;
		}
	    }
	  else if(e != &(se->eee[0]) && t==TCH2)
	    {
	      e0 = &(se->eee[1]);
	      e2 = &(se->eee[2]);
	      if(set_stoch_vars2(e,e0,e2,p,t,pi_brexit))
		{
		  return 1;
		}
	    }
	  else
	    {
	      if(set_vars2(e,p,t,0))
		{
		  fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
		  return 1;
		}
	    }
	}
    }

  for(ih=0; ih<NHIST_CH; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(t=t0; t<(NT+1); t++)
	{
	  set_vars3(e,p,t,0);
	}
    }

  nx = 0;
  // first we do the equations for the periods after announcement but before reform
  // during these periods, we only have to evaluate equilibrium conditions for one history
  // (all histories are the same up until the reform occurs (or not)
  //
  // there is no uncertainty except in the period immediately before the reform, so we do all the usual equations
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);
  for(t=t0; t<(TCH1-1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      myf[nx] = noarb(p,e,t,i);
	      nx=nx+1;

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
	    }
	}
    }

  // in the period before the reform, we have to deal with the uncertainty about whether we will get
  // no reform (pi_vote) or reform (1-pi_vote)
  t = TCH1-1;
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);
  e2 = &(se->eee[1]);
  p2 = &(sp->ppp[1]);

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  for(i=0; i<NC; i++)
    {
      if(i!=NC-1)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      myf[nx] = mkt_clear_K(e,t,i);
      nx=nx+1;

      myf[nx] = pi_vote * euler(p,e,t,i) + (1.0-pi_vote) * euler(p2,e2,t,i);
      nx = nx+1;

      myf[nx] = pi_vote * noarb(p,e,t,i) + (1.0-pi_vote) * noarb(p2,e2,t,i);
      nx = nx+1;

      if(i>0)
	{
	  myf[nx] = fin_mkt(e,t,i);
	  nx=nx+1;
	}

      for(j=0; j<NC; j++)
	{	  
	  myf[nx] = mkt_clear_Y2(e,t,i,j);
	  nx=nx+1;
	}
    }

  // next we do the equations for the "no reform" path... no uncertainty here and it goes all the way to the steady state
  t0 = TCH1;
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);

  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      if(t<(NT-1))
		{
		  myf[nx] = noarb(p,e,t,i);
		  nx=nx+1;
		}

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
	    }
	}
    }

  // now we do the TCH2-TCH1 periods between vote and Brexit... here there is uncertainty again in the
  // period before reversion only
  t0 = TCH1;
  e = &(se->eee[1]);
  p = &(sp->ppp[1]);
  for(t=t0; t<(TCH2-1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      myf[nx] = noarb(p,e,t,i);
	      nx=nx+1;

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
	    }
	}
    }

  // in the period before reversion is possible, we have to deal with uncertainty
  t = TCH2-1;
  e = &(se->eee[1]);
  p = &(sp->ppp[1]);
  e2 = &(se->eee[2]);
  p2 = &(sp->ppp[2]);

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  for(i=0; i<NC; i++)
    {
      if(i!=NC-1)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      myf[nx] = mkt_clear_K(e,t,i);
      nx=nx+1;

      myf[nx] = pi_brexit * euler(p,e,t,i) + (1.0-pi_brexit) * euler(p2,e2,t,i);
      nx = nx+1;

      myf[nx] = pi_brexit * noarb(p,e,t,i) + (1.0-pi_brexit) * noarb(p2,e2,t,i);
      nx = nx+1;

      if(i>0)
	{
	  myf[nx] = fin_mkt(e,t,i);
	  nx=nx+1;
	}

      for(j=0; j<NC; j++)
	{	  
	  myf[nx] = mkt_clear_Y2(e,t,i,j);
	  nx=nx+1;
	}
    }

  // after reversion (or not), we do the histories separately
  t0 = TCH2;
  for(ih=1; ih<NHIST; ih++)
    {
      e = &(se->eee[ih]);
      p = &(sp->ppp[ih]);

      for(t=t0; t<(NT+1); t++)
	{
	  myf[nx] = price_norm(e,t);
	  nx=nx+1;

	  for(i=0; i<NC; i++)
	    {
	      if(i!=NC-1)
		{
		  myf[nx] = bop(p,e,t,i);
		  nx = nx+1;
		}

	      myf[nx] = muc_mul(p,e,t,i);
	      nx=nx+1;

	      myf[nx] = mkt_clear_K(e,t,i);
	      nx=nx+1;

	      if(t<NT)
		{
		  myf[nx] = euler(p,e,t,i);
		  nx = nx+1;

		  if(t<(NT-1))
		    {
		      myf[nx] = noarb(p,e,t,i);
		      nx=nx+1;
		    }

		  if(i>0)
		    {
		      myf[nx] = fin_mkt(e,t,i);
		      nx=nx+1;
		    }
		}

	      for(j=0; j<NC; j++)
		{	  
		  myf[nx] = mkt_clear_Y2(e,t,i,j);
		  nx=nx+1;
		}
	    }
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error evaluating stoch_eqm_ch eqns! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }
  
  return 0;
}
#else
uint eval_stoch_eqm_ch_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,j=0,t=NT,nx=0,ih=0;

  stoch_eqm_ch * se = &(sss_ch[tn]);
  stoch_params_ch * sp = &(sppp_ch[tn]);
  eqm *e, *e0, *e2;
  params *p, *p2;

  unstack_stoch_eqm_ch_vars(se,myx);

  uint t0 = TCH0;

  // setvars1
  for(ih=0; ih<NHIST_CH; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(t=t0; t<(NT+1); t++)
	{
	  if(set_vars1(e,p,t,0))
	    {
	      fprintf(logfile,KRED "Error calling set_vars1!\n" RESET);
	      return 1;
	    }
	}
    }

  // do setvars2 in chunks...
  // first the chunk where all branches are separate
  for(t=NT; t>=TCH2+1; t--)
    {
      for(ih=0; ih<NHIST_CH; ih++)
	{
	  e = &( (se->eee)[ih] );
	  p = &( (sp->ppp)[ih] );
	}
      if(set_vars2(e,p,t,0))
	{
	  fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
	  return 1;
	}
    }

  // next the part where the branches overlap
  for(t=TCH2; t>=t0; t--)
    {
      // before the period where reform can occur or not, all branches get the usual treatment
      if(t<TCH1)
	{
	  for(ih=0; ih<NHIST_CH; ih++)
	    {
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );   
	      if(set_vars2(e,p,t,0))
		{
		  fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
		  return 1;
		}
	    }
	}

      // in the period where the uncertainty about whether the reform happens or not gets resolved,
      // all branches are the same
      else if(t==TCH1)
	{
	  e0 = &(se->eee[0]);
	  e2 = &(se->eee[1]);
	  for(ih=0; ih<NHIST_CH; ih++)
	    {
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );    
	      if(set_stoch_vars2(e,e0,e2,p,t,pi_vote))
		{
		  return 1;
		}
	    }
	}

      // in subsequent periods, only some branches overlap
      else
	{
	  // count the branches that get dealth with to check...
	  uint hcnt=0;

	  // branch 0 is for sure separate
	  e = &( (se->eee)[0] );
	  p = &( (sp->ppp)[0] );
	  if(set_vars2(e,p,t,0))
	    {
	      fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
	      return 1;
	    }
	  hcnt=hcnt+1;

	  // other branches that have also separated already
	  for(ih=1; ih<=t-TCH1-1; ih++)
	    {
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );
	      if(set_vars2(e,p,t,0))
		{
		  fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
		  return 1;
		}
	      hcnt=hcnt+1;
	    }

	  // rest of the branches get the uncertainty treatment
	  ih = t-TCH1+1;
	  e0 = &(se->eee[ih]);
	  e2 = &(se->eee[ih-1]);
	  if(set_stoch_vars2(e0,e0,e2,p,t,pi_brexit))
	    {
	      fprintf(logfile,KRED "Error calling set_stoch_vars2!\n" RESET);
	      return 1;
	    }
	  hcnt=hcnt+1;
	  if(set_stoch_vars2(e2,e0,e2,p,t,pi_brexit))
	    {
	      fprintf(logfile,KRED "Error calling set_stoch_vars2!\n" RESET);
	      return 1;
	    }
	  hcnt=hcnt+1;

	  for(ih=t-TCH1+2; ih<NHIST_CH; ih++)
	    {	      
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );
	      if(set_stoch_vars2(e,e0,e2,p,t,pi_brexit))
		{
		  fprintf(logfile,KRED "Error calling set_stoch_vars2!\n" RESET);
		  return 1;
		}
	      hcnt=hcnt+1;
	    }

	  if(hcnt!=NHIST_CH)
	    {
	      fprintf(logfile,KRED "Error called set_stoch_vars2 wrong number of tines!\n");
	      return 1;
	    }
	}
    }
  
  // setvars3
  for(ih=0; ih<NHIST_CH; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      
      for(t=TCH0; t<(NT+1); t++)
	{
	  set_vars3(e,p,t,0);
	}
    }

  nx = 0;

  // first we do the equations for the periods after announcement but before reform
  // during these periods, we only have to evaluate equilibrium conditions for one history
  // (all histories are the same up until the reform occurs (or not)
  //
  // there is no uncertainty except in the period immediately before the reform, so we do all the usual equations
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);
  for(t=t0; t<(TCH1-1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      myf[nx] = noarb(p,e,t,i);
	      nx=nx+1;

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
	    }
	}
    }

  // in the period before the reform, we have to deal with the uncertainty about whether we will get
  // no reform (pi_vote) or reform (1-pi_vote)
  t = TCH1-1;
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);
  e2 = &(se->eee[1]);
  p2 = &(sp->ppp[1]);

  myf[nx] = price_norm(e,t);
  nx=nx+1;

  for(i=0; i<NC; i++)
    {
      if(i!=NC-1)
	{
	  myf[nx] = bop(p,e,t,i);
	  nx = nx+1;
	}

      myf[nx] = muc_mul(p,e,t,i);
      nx=nx+1;

      myf[nx] = mkt_clear_K(e,t,i);
      nx=nx+1;

      myf[nx] = pi_vote * euler(p,e,t,i) + (1.0-pi_vote) * euler(p2,e2,t,i);
      nx = nx+1;

      myf[nx] = pi_vote * noarb(p,e,t,i) + (1.0-pi_vote) * noarb(p2,e2,t,i);
      nx = nx+1;

      if(i>0)
	{
	  myf[nx] = fin_mkt(e,t,i);
	  nx=nx+1;
	}

      for(j=0; j<NC; j++)
	{	  
	  myf[nx] = mkt_clear_Y2(e,t,i,j);
	  nx=nx+1;
	}
    }

  // next we do the equations for the "no reform" path... no uncertainty here and it goes all the way to the steady state
  t0 = TCH1;
  e = &(se->eee[0]);
  p = &(sp->ppp[0]);

  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      nx = nx+1;

	      if(t<(NT-1))
		{
		  myf[nx] = noarb(p,e,t,i);
		  nx=nx+1;
		}

	      if(i>0)
		{
		  myf[nx] = fin_mkt(e,t,i);
		  nx=nx+1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;
	    }
	}
    }

  // next we do the equations where other branches have separated... no uncertainty here either
  for(ih=1; ih<NHIST_CH; ih++)
    {
      t0 = TCH1+ih;
      if(t0>TCH2)
	{
	  t0=TCH2;
	}
      e = &(se->eee[ih]);
      p = &(sp->ppp[ih]);
     
      for(t=t0; t<(NT+1); t++)
	{
	  myf[nx] = price_norm(e,t);
	  nx=nx+1;

	  for(i=0; i<NC; i++)
	    {
	      if(i!=NC-1)
		{
		  myf[nx] = bop(p,e,t,i);
		  nx = nx+1;
		}

	      myf[nx] = muc_mul(p,e,t,i);
	      nx=nx+1;

	      myf[nx] = mkt_clear_K(e,t,i);
	      nx=nx+1;

	      if(t<NT)
		{
		  myf[nx] = euler(p,e,t,i);
		  nx = nx+1;

		  if(t<(NT-1))
		    {
		      myf[nx] = noarb(p,e,t,i);
		      nx=nx+1;
		    }

		  if(i>0)
		    {
		      myf[nx] = fin_mkt(e,t,i);
		      nx=nx+1;
		    }
		}

	      for(j=0; j<NC; j++)
		{	  
		  myf[nx] = mkt_clear_Y2(e,t,i,j);
		  nx=nx+1;
		}
	    }
	}
    }

  // in the periods where reversion is possible, we have to deal with uncertainty
  for(t=TCH2-1; t>=TCH1; t--)
    {
      ih=t-TCH1+2;
      e = &(se->eee[ih]);
      p = &(sp->ppp[ih]);
      e2 = &(se->eee[ih-1]);
      p2 = &(sp->ppp[ih-1]);
      
      if(nx==12276){gdb_test_func();}
      myf[nx] = price_norm(e,t);
      nx=nx+1;

      for(i=0; i<NC; i++)
	{
	  if(i!=NC-1)
	    {
	      if(nx==12276){gdb_test_func();}
	      myf[nx] = bop(p,e,t,i);
	      nx = nx+1;
	    }

	  if(nx==12276){gdb_test_func();}
	  myf[nx] = muc_mul(p,e,t,i);
	  nx=nx+1;

	  if(nx==12276){gdb_test_func();}
	  myf[nx] = mkt_clear_K(e,t,i);
	  nx=nx+1;

	  if(nx==12276){gdb_test_func();}
	  myf[nx] = pi_brexit * euler(p,e,t,i) + 
	    (1.0-pi_brexit) * euler(p2,e2,t,i);
	  nx = nx+1;

	  if(nx==12276){gdb_test_func();}
	  myf[nx] = pi_brexit * noarb(p,e,t,i) + 
	    (1.0-pi_brexit) * noarb(p2,e2,t,i);
	  nx = nx+1;

	  if(i>0)
	    {
	      if(nx==12276){gdb_test_func();}
	      myf[nx] = fin_mkt(e,t,i);
	      nx=nx+1;
	    }

	  for(j=0; j<NC; j++)
	    {	  
	      if(nx==12276){gdb_test_func();}
	      myf[nx] = mkt_clear_Y2(e,t,i,j);
	      nx=nx+1;

	    }
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error evaluating stoch_eqm_ch eqns! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }
  
  return 0;
}
#endif

uint solve_bgp(double bb[NC-1])
{
  bbgp[0] = bb[0];
  bbgp[1] = bb[1];
  solver_n = nbgp;
  alloc_solver_mem();
  set_initial_bgp_guess();
  par=0;
  gsl_multiroot_function_fdf f = {&bgp_func_f,&bgp_func_df,&bgp_func_fdf,nbgp,NULL};
  uint status = find_root_deriv_mkl(&f);
  free_solver_mem();
  
  return status;
}

uint solve_eqm()
{
  char * sname;
  if(scenario < 6)
    {
      if(nokappa)
	{
	  sname = "output/seed_nokappa.bin";
	}
      else if(eqkappa)
	{
	  sname = "output/seed_eqkappa.bin";
	}
      else if(psi_sens)
	{
	  sname = "output/seed_spsi.bin";
	}
      else if(zeta_sens)
	{
	  sname = "output/seed_szeta.bin";
	}
      else
	{
	  sname = "output/seed_baseline.bin";
	}
    }
  else
    {
      if(nokappa)
	{
	  sname = "output/seed_ch_nokappa.bin";
	}
      else if(eqkappa)
	{
	  sname = "output/seed_ch_eqkappa.bin";
	}
      else if(psi_sens)
	{
	  sname = "output/seed_ch_spsi.bin";
	}
      else if(zeta_sens)
	{
	  sname = "output/seed_ch_szeta.bin";
	}
      else if(higher_ch_trd_costs)
	{
	  sname = "output/seed_ch_hi.bin";
	}
      else
	{
	  sname = "output/seed_ch_baseline.bin";
	}
    }

  // if we are solving for the no-Brexit (or no-reform in China-like exercise) counterfactual we must construct an initial guess for the solver
  if(scenario==0 || scenario==6)
    {
      if(read_seed==1)
	{
	  free_solver_mem();
	  solver_n = neqm;
	  alloc_solver_mem();

	  fprintf(logfile,KMAG "\n\tReading from seed file %s\n" RESET,sname);

	  if(scenario==6)
	    {
	      fprintf(logfile,KMAG "\n\tChina-like scenario: must construct guess first to reset initial conditions %s\n" RESET,sname);
	      if(set_initial_eqm_guess())
		{
		  fprintf(logfile,KRED "Error constructing equilibrium guess!\n" RESET);
		  free_solver_mem();
		  return 1;
		}
	    }
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
      stoch_eqm_ch * se_ch = &(sss_ch[0]);
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
#ifdef CH_SIMPLE
      if(scenario == 7)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1_ch[it]), &(eee0_ch[0])  );
	    }
	}
      if(scenario == 8)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee2_ch[it]), &(eee0_ch[0])  );
	    }
	}
#else
      if(scenario>=1000 && scenario<2000)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1_ch[scenario-1000][it]), &(eee0_ch[0])  );
	    }	  
	}
#endif
#ifdef CH_SIMPLE
      if(scenario == 10)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1_ch[it]), &((se_ch->eee)[0])  );
	    }
	}
      if(scenario == 11)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee2_ch[it]), &((se_ch->eee)[0])  );
	    }
	}
#else
      if(scenario>=2000)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1_ch[scenario-2000][it]), &((se_ch->eee)[0])  );
	    }	  
	}
#endif

      // now determine which eqm struct to use as the initial guess...
      // ...if we are in the optimistic scenario, use the first path...
      if(scenario == 2 || scenario == 4)
	{
	  e = &(eee0[0]);
	}
      // ...otherwise use the second path
      else if(scenario == 3 || scenario == 5)
	{
	  e = &(eee0[0]);
	}
      else if(scenario == 7 || scenario == 8)
	{
	  e = &(eee0_ch[0]);
	}
      else if(scenario == 10 || scenario == 11)
	{
	  e = &(eee0_ch[0]);
	}
      else
	{
	  e = &(eee0_ch[0]);
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
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      par=1;
      gsl_multiroot_function_fdf f = {&eqm_func_f,&eqm_func_df,&eqm_func_fdf,neqm,NULL};
      status = find_root_deriv_mkl(&f);
      if(status)
	{
	  fprintf(logfile,KRED "Error solving for equilibrium!\n" RESET);
	}

      if((scenario==0 || scenario==6) && write_seed==1 && !status)
	{
	  write_vec_bin(solver_x->data, neqm, sname);
	}
    }

  free_solver_mem();

  return status;
}

uint solve_stoch_eqm()
{  
  char * sname;
  if(nokappa)
    {
      sname = "output/stoch_seed_nokappa.bin";
    }
  else if(eqkappa)
    {
      sname = "output/stoch_seed_eqkappa.bin";
    }
  else if(low_pi)
    {
      sname = "output/stoch_seed_lowpi.bin";
    }
  else if(high_pi)
    {
      sname = "output/stoch_seed_highpi.bin";
    }
  else if(fin_aut)
    {
      sname = "output/stoch_seed_finaut.bin";
    }
  else if(comp_mkts)
    {
      sname = "output/stoch_seed_compmkts.bin";
    }
  else if(zeta_sens)
    {
      sname = "output/stoch_seed_szeta.bin";
    }
  else if(psi_sens)
    {
      sname = "output/stoch_seed_spsi.bin";
    }
  else
    {
      sname = "output/stoch_seed_baseline.bin";
    }


  // before we do anything, we need to copy all variables from the no-Brexit counterfactual to every 
  // thread-specific copy of every branch of the stochastic equilibrium structure...
  // i.e. copy everything from eee[0] to the stochastic equilibrium sss[it] for each
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

  // reallocate memory for solver
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

  // if we are going to read the initial guess from a seed file, do so
  if(read_stoch_seed)
    {
      if(read_vec_bin(solver_x->data, neqm, sname))
	{
	  fprintf(logfile,KRED "Error loading stochastic equilibrium guess from seed file!\n" RESET);
	  free_solver_mem();
	  return 1;
	}
    }
  // otherwise, we should use the solutions from the first round of deterministic exercises
  // (the ones from TREF onward) as the initial guess...
  // guess that the no-vote path follows the no-Brexit counterfactual, the yes-vote + soft-brexit
  // path follows the determinstic soft-Brexit equilibrium, and the yes-vote + hard-brexit path
  // follows the deterministic hard-Brexit equilibrium
  else
    {
      copy_vars( &( (&(sss[0]))->eee[0] ) , &(eee0[0]) );
      copy_vars( &( (&(sss[0]))->eee[1] ) , &(eee1[0]) );
      copy_vars( &( (&(sss[0]))->eee[2] ) , &(eee2[0]) );

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
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      par=1;
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

#ifdef CH_SIMPLE
uint solve_stoch_eqm_ch()
{  
  char * sname;
  if(nokappa)
    {
      sname = "output/stoch_seed_ch_nokappa.bin";
    }
  else if(eqkappa)
    {
      sname = "output/stoch_seed_ch_eqkappa.bin";
    }
  else if(low_pi)
    {
      sname = "output/stoch_seed_ch_lowpi.bin";
    }
  else if(high_pi)
    {
      sname = "output/stoch_seed_ch_highpi.bin";
    }
  else if(fin_aut)
    {
      sname = "output/stoch_seed_ch_finaut.bin";
    }
  else if(comp_mkts)
    {
      sname = "output/stoch_seed_ch_compmkts.bin";
    }
  else if(zeta_sens)
    {
      sname = "output/stoch_seed_ch_szeta.bin";
    }
  else if(psi_sens)
    {
      sname = "output/stoch_seed_ch_spsi.bin";
    }
  else if(higher_ch_trd_costs)
    {
      sname = "output/stoch_seed_ch_hi.bin";
    }
  else
    {
      sname = "output/stoch_seed_ch_baseline.bin";
    }

  // before we do anything, we need to copy all variables from counterfactual to every 
  // thread-specific copy of every branch of the stochastic equilibrium structure...
  uint it,ih;
  for(it=0; it<NTH; it++)
    {
      for(ih=0; ih<NHIST_CH; ih++)
	{
	  copy_vars( &( (&(sss_ch[it]))->eee[ih] ) , &(eee0_ch[0]) );
	}
    }

  // reallocate memory for solver
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

  // if we are going to read the initial guess from a seed file, do so
  if(read_stoch_seed)
    {
      if(read_vec_bin(solver_x->data, neqm, sname))
	{
	  fprintf(logfile,KRED "Error loading stochastic equilibrium guess from seed file!\n" RESET);
	  free_solver_mem();
	  return 1;
	}
    }
  // otherwise, we should use the solutions from the first round of deterministic exercises
  // (the ones from TCH0 onward) as the initial guess...
  else
    {
      copy_vars( &( (&(sss_ch[0]))->eee[0] ) , &(eee0_ch[0]) );
      copy_vars( &( (&(sss_ch[0]))->eee[1] ) , &(eee2_ch[0]) );
      copy_vars( &( (&(sss_ch[0]))->eee[2] ) , &(eee1_ch[0]) );

      if(stack_stoch_eqm_ch_vars(solver_x->data,&(sss_ch[0])))
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
      status = stoch_eqm_ch_func_f(solver_x,NULL,f0[0]);
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      par=1;
      gsl_multiroot_function_fdf f = {&stoch_eqm_ch_func_f,
				      &stoch_eqm_ch_func_df,
				      &stoch_eqm_ch_func_fdf,neqm,NULL};
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
#else
uint solve_stoch_eqm_ch()
{  
  char * sname;
  if(nokappa)
    {
      sname = "output/stoch_seed_ch_complex_nokappa.bin";
    }
  else if(eqkappa)
    {
      sname = "output/stoch_seed_ch_complex_eqkappa.bin";
    }
  else if(low_pi)
    {
      sname = "output/stoch_seed_ch_complex_lowpi.bin";
    }
  else if(high_pi)
    {
      sname = "output/stoch_seed_ch_complex_highpi.bin";
    }
  else if(fin_aut)
    {
      sname = "output/stoch_seed_ch_complex_finaut.bin";
    }
  else if(comp_mkts)
    {
      sname = "output/stoch_seed_ch_complex_compmkts.bin";
    }
  else if(zeta_sens)
    {
      sname = "output/stoch_seed_ch_complex_szeta.bin";
    }
  else if(psi_sens)
    {
      sname = "output/stoch_seed_ch_complex_spsi.bin";
    }
  else if(higher_ch_trd_costs)
    {
      sname = "output/stoch_seed_ch_complex_hi.bin";
    }
  else
    {
      sname = "output/stoch_seed_ch_complex_baseline.bin";
    }

  // before we do anything, we need to copy all variables from counterfactual to every 
  // thread-specific copy of every branch of the stochastic equilibrium structure...
  uint it,ih;
  for(it=0; it<NTH; it++)
    {
      for(ih=0; ih<NHIST_CH; ih++)
	{
	  copy_vars( &( (&(sss_ch[it]))->eee[ih] ) , &(eee0_ch[0]) );
	}
    }

  // reallocate memory for solver
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

  // if we are going to read the initial guess from a seed file, do so
  if(read_stoch_seed)
    {
      if(read_vec_bin(solver_x->data, neqm, sname))
	{
	  fprintf(logfile,KRED "Error loading stochastic equilibrium guess from seed file!\n" RESET);
	  free_solver_mem();
	  return 1;
	}
    }
  // otherwise, we should use the solutions from the first round of deterministic exercises
  // (the ones from TCH0 onward) as the initial guess...
  else
    {
      uint ix;
      copy_vars( &( (&(sss_ch[0]))->eee[0] ) , &(eee0_ch[0]) );

      for(ix=0; ix<NHIST_CH-1; ix++)
	{
	  copy_vars( &( (&(sss_ch[0]))->eee[ix+1] ) , &(eee1_ch[ix][0]) );
	}

      if(stack_stoch_eqm_ch_vars(solver_x->data,&(sss_ch[0])))
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
      status = stoch_eqm_ch_func_f(solver_x,NULL,f0[0]);
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      par=0;
      gsl_multiroot_function_fdf f = {&stoch_eqm_ch_func_f,
				      &stoch_eqm_ch_func_df,
				      &stoch_eqm_ch_func_fdf,neqm,NULL};
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
#endif

///////////////////////////////
void calc_welfare(eqm * e, const params * p)
{
  int t, i;
  for(i=0; i<NC; i++)
    {
      t=NT;
      e->welfare_t[t][i] = (1.0/(1.0-p->beta)) * 
	(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      
      for(t=(NT-1); t>=0; t--)
	{
	    e->welfare_t[t][i] = p->beta * e->welfare_t[t+1][i] + 
	      (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
	      
	    e->welfare_t[t+1][i] = pow(e->welfare_t[t+1][i],1.0/p->psi);
	}
      e->welfare_t[0][i] = pow(e->welfare_t[0][i],1.0/p->psi);
    }

  if(scenario == 0)
    {
      for(i=0; i<NC; i++)
	{
	  t=NT;
	  e->welfare_cost_t[t][i] = (1.0/(1.0-e->Q_t[t][i])) * e->P_t[t][i] * e->C_t[t][i];
	    
	  for(t=(NT-1); t>=0; t--)
	    {
	      e->welfare_cost_t[t][i] = e->P_t[t][i]*e->C_t[t][i] + e->Q_t[t][i] * e->welfare_cost_t[t+1][i];
	    }
	  for(t=0; t<NT; t++)
	    {
	      e->welfare_cost_t[t][i] = e->Q_t[t][i] * e->welfare_cost_t[t+1][i];
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
	  e->welfare2_t[t][i] = (1.0/(1.0-p->beta)) * 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      
	  for(t=(NT-1); t>=TBREXIT; t--)
	    {
	            e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		      (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
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
      e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
	(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
    }

  // last pre-Brexit period for "leave" branch
  t = TBREXIT-1;
  e = &( (se->eee)[1] );
  p = &( (sp->ppp)[1] );
  e2 = &( (se->eee)[2] );
  for(i=0; i<NC; i++)
    {
      e->welfare2_t[t][i] = (pi_brexit * p->beta * e->welfare2_t[t+1][i]) + 
	(1.0-pi_brexit) * (p->beta*e2->welfare2_t[t+1][i]) + 
	(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      e2->welfare2_t[t][i] = e->welfare2_t[t][i];

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
	            e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		      (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
		    e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	}
    }

  // vote period (all branches)
  // Brexit period for "leave" branch
  t = TVOTE-1;
  e0 = &( (se->eee)[0] );
  e2 = &( (se->eee)[1] );
  for(ih=0; ih<NHIST; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      
      for(i=0; i<NC; i++)
	{
	  e->welfare2_t[t][i] = (pi_vote * p->beta * e0->welfare2_t[t+1][i]) + 
	    (1.0-pi_vote) * (p->beta*e2->welfare2_t[t+1][i]) + 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * 
	     pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));	  
	}
    }
  for(ih=0; ih<NHIST; ih++)
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
	            e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		      (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
		    e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	  e->welfare2_t[0][i] = pow(e->welfare2_t[0][i],1.0/p->psi);
	}
    }
}

# ifdef CH_SIMPLE
void calc_stoch_ch_welfare()
{
  int t, i, ih;
  eqm *e, *e0, *e2;
  const params *p;

  stoch_eqm_ch * se = &(sss_ch[0]);
  const stoch_params_ch * sp = &(sppp_ch[0]);

  for(ih=0; ih<NHIST_CH; ih++)
    { 
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(i=0; i<NC; i++)
	{
	  t=NT;
	  e->welfare2_t[t][i] = (1.0/(1.0-p->beta)) * 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      
	  for(t=(NT-1); t>=TCH2; t--)
	    {
	            e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		      (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
		    e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	}
    }

  t=TCH2-1;
  e = &( (se->eee)[0] );
  p = &( (sp->ppp)[0] );
  for(i=0; i<NC; i++)
    {      
      e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
	(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
    }

  t = TCH2-1;
  e = &( (se->eee)[1] );
  p = &( (sp->ppp)[1] );
  e2 = &( (se->eee)[2] );
  for(i=0; i<NC; i++)
    {
      e->welfare2_t[t][i] = (pi_brexit * p->beta * e->welfare2_t[t+1][i]) + 
	(1.0-pi_brexit) * (p->beta*e2->welfare2_t[t+1][i]) + 
	(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      e2->welfare2_t[t][i] = e->welfare2_t[t][i];

      e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
      e2->welfare2_t[t+1][i] = pow(e2->welfare2_t[t+1][i],1.0/p->psi);
    }

  for(ih=0; ih<NHIST_CH; ih++)
    { 
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      
      for(i=0; i<NC; i++)
	{
	  for(t=TCH2-2; t>=TCH1; t--)
	    {
	            e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		      (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
		    e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	}
    }

  t = TCH1-1;
  e0 = &( (se->eee)[0] );
  e2 = &( (se->eee)[1] );
  for(ih=0; ih<NC; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      
      for(i=0; i<NC; i++)
	{
	  e->welfare2_t[t][i] = (pi_vote * p->beta * e0->welfare2_t[t+1][i]) + 
	    (1.0-pi_vote) * (p->beta*e2->welfare2_t[t+1][i]) + 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * 
	     pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
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

  for(ih=0; ih<NHIST_CH; ih++)
    { 
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(i=0; i<NC; i++)
	{      
	  for(t=(TCH1-2); t>=0; t--)
	    {
	            e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		      (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
		    e->welfare2_t[t+1][i] = pow(e->welfare2_t[t+1][i],1.0/p->psi);
	    }
	  e->welfare2_t[0][i] = pow(e->welfare2_t[0][i],1.0/p->psi);
	}
    }
}
#else
void calc_stoch_ch_welfare()
{
  int t, i, ih;
  eqm *e, *e0, *e2;
  const params *p;

  stoch_eqm_ch * se = &(sss_ch[0]);
  const stoch_params_ch * sp = &(sppp_ch[0]);

  // first do the parts for each branch after all uncertainty has been resolved
  for(ih=0; ih<NHIST_CH; ih++)
    { 
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );

      for(i=0; i<NC; i++)
	{
	  t=NT;
	  e->welfare2_t[t][i] = (1.0/(1.0-p->beta)) * 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
      
	  for(t=(NT-1); t>=TCH2; t--)
	    {
	      e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
	    }
	}
    }

  // now do the uncertainty during TCH1-TCH2
  for(t=TCH2-1; t>=TCH1; t--)
    {
      // in these periods, branch 0 has already split off
      e = &((se->eee)[0]);
      p = &( (sp->ppp)[0] );
      for(i=0; i<NC; i++)
	{
	  e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
	}

      // some branches have also already split off
      for(ih=1; ih<(t-TCH1+1); ih++)
	{
	  e = &((se->eee)[ih]);
	  p = &( (sp->ppp)[ih] );
	  for(i=0; i<NC; i++)
	    {
	      e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));	  
	    }
	}

      // now do the branches where we have to deal with the uncertainty
      ih=t-TCH1+2;
      e = &((se->eee)[ih]);
      p = &( (sp->ppp)[ih] );
      e2 = &((se->eee)[ih-1]);

      for(i=0; i<NC; i++)
	{
	  e->welfare2_t[t][i] = (pi_brexit * p->beta * e->welfare2_t[t+1][i]) + 
	    (1.0-pi_brexit) * (p->beta*e2->welfare2_t[t+1][i]) + 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
	  e2->welfare2_t[t][i] = e->welfare2_t[t][i];
	}
      
      // now do the
      for(ih=t-TCH1+3; ih<NHIST_CH; ih++)
	{
	  e2 = &((se->eee)[ih]);
	  for(i=0; i<NC; i++)
	    {
	      e2->welfare2_t[t][i] = e->welfare2_t[t][i];
	    }
	}
    }

  // now do the last bit where all branches are the same
  t = TCH1-1;
  e0 = &( (se->eee)[0] );
  e2 = &( (se->eee)[1] );
  for(ih=0; ih<NHIST_CH; ih++)
    {
      e = &( (se->eee)[ih] );
      p = &( (sp->ppp)[ih] );
      
      for(i=0; i<NC; i++)
	{
	  e->welfare2_t[t][i] = (pi_vote * p->beta * e0->welfare2_t[t+1][i]) + 
	    (1.0-pi_vote) * (p->beta*e2->welfare2_t[t+1][i]) + 
	    (pow(e->C_t[t][i],p->phi[i]*p->psi) * 
	     pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
	}
    }

  // now even earlier periods
  for(t=TCH1-2; t>=0; t--)
    {
      for(i=0; i<NC; i++)
	{
	  for(ih=0; ih<NHIST_CH; ih++)
	    { 
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );
	      e->welfare2_t[t][i] = p->beta * e->welfare2_t[t+1][i] + 
		(pow(e->C_t[t][i],p->phi[i]*p->psi) * pow(p->Lbar[i]-e->L_t[t][i],(1.0-p->phi[i])*p->psi));
	    }	  
	}
    }

  // now go through and do the 1/psi bit
  for(t=0; t<NT+1; t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(ih=0; ih<NHIST_CH; ih++)
	    {
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );
	      e->welfare2_t[t][i] = pow(e->welfare2_t[t][i],1.0/p->psi);
	    }
	}
    }
}
#endif
    
void calc_ce_welfare1()
{
  int t, i;

  stoch_eqm *se = &(sss[0]);
  const stoch_params *sp = &(sppp[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[1] );
  //const eqm *e_stay = &( (se->eee)[0] );
  const eqm *e_stay = &(eee0[0]);
  const eqm *e_soft = &(eee1[0]);
  const eqm * e_hard = &(eee2[0]);

  double C_stay, C_soft, C_hard, C;
  double L_stay, L_soft, L_hard, L;

  for(i=0; i<NC; i++)
    {
      t=NT;

      C_stay = e_stay->C_t[t][i];
      L_stay = e_stay->L_t[t][i];
      
      C_soft = e_soft->C_t[t][i];
      L_soft = e_soft->L_t[t][i];
      
      C_hard = e_hard->C_t[t][i];
      L_hard = e_hard->L_t[t][i];

      C = pi_vote*C_stay +
	(1.0-pi_vote)*pi_brexit*C_soft + 
	(1.0-pi_vote)*(1.0-pi_brexit)*C_hard;
      
      L = pi_vote*L_stay +
	(1.0-pi_vote)*pi_brexit*L_soft + 
	(1.0-pi_vote)*(1.0-pi_brexit)*L_hard;
      
      e->welfare3_t[t][i] = (1.0/(1.0-p->beta)) * 
	(pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
      
      for(t=(NT-1); t>=TREF; t--)
	{
	  C_stay = e_stay->C_t[t][i];
	  L_stay = e_stay->L_t[t][i];

	  C_soft = e_soft->C_t[t][i];
	  L_soft = e_soft->L_t[t][i];

	  C_hard = e_hard->C_t[t][i];
	  L_hard = e_hard->L_t[t][i];

	  C = pi_vote*C_stay +
	    (1.0-pi_vote)*pi_brexit*C_soft + 
	    (1.0-pi_vote)*(1.0-pi_brexit)*C_hard;

	  L = pi_vote*L_stay +
	    (1.0-pi_vote)*pi_brexit*L_soft + 
	    (1.0-pi_vote)*(1.0-pi_brexit)*L_hard;

	  e->welfare3_t[t][i] = p->beta * e->welfare3_t[t+1][i] + 
	    (pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
	  
	  e->welfare3_t[t+1][i] = pow(e->welfare3_t[t+1][i],1.0/p->psi);
	}
      e->welfare3_t[TREF][i] = pow(e->welfare3_t[TREF][i],1.0/p->psi);
    }
}

#ifdef CH_SIMPLE
void calc_ce_ch_welfare1()
{
  int t, i;

  stoch_eqm_ch *se = &(sss_ch[0]);
  const stoch_params_ch *sp = &(sppp_ch[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[2] );
  //const eqm *e_stay = &( (se->eee)[0] );
  const eqm *e_stay = &(eee0_ch[0]);
  const eqm *e_temp = &(eee2_ch[0]);
  const eqm * e_perm = &(eee1_ch[0]);

  double C_stay, C_temp, C_perm, C;
  double L_stay, L_temp, L_perm, L;

  for(i=0; i<NC; i++)
    {
      t=NT;

      C_stay = e_stay->C_t[t][i];
      L_stay = e_stay->L_t[t][i];
      
      C_temp = e_temp->C_t[t][i];
      L_temp = e_temp->L_t[t][i];
      
      C_perm = e_perm->C_t[t][i];
      L_perm = e_perm->L_t[t][i];

      C = pi_vote*C_stay +
	(1.0-pi_vote)*pi_brexit*C_temp + 
	(1.0-pi_vote)*(1.0-pi_brexit)*C_perm;
      
      L = pi_vote*L_stay +
	(1.0-pi_vote)*pi_brexit*L_temp + 
	(1.0-pi_vote)*(1.0-pi_brexit)*L_perm;
      
      e->welfare3_t[t][i] = (1.0/(1.0-p->beta)) * 
	(pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
      
      for(t=(NT-1); t>=TCH0; t--)
	{
	  C_stay = e_stay->C_t[t][i];
	  L_stay = e_stay->L_t[t][i];

	  C_temp = e_temp->C_t[t][i];
	  L_temp = e_temp->L_t[t][i];

	  C_perm = e_perm->C_t[t][i];
	  L_perm = e_perm->L_t[t][i];

	  C = pi_vote*C_stay +
	    (1.0-pi_vote)*pi_brexit*C_temp + 
	    (1.0-pi_vote)*(1.0-pi_brexit)*C_perm;

	  L = pi_vote*L_stay +
	    (1.0-pi_vote)*pi_brexit*L_temp + 
	    (1.0-pi_vote)*(1.0-pi_brexit)*L_perm;

	  e->welfare3_t[t][i] = p->beta * e->welfare3_t[t+1][i] + 
	    (pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
	  
	  e->welfare3_t[t+1][i] = pow(e->welfare3_t[t+1][i],1.0/p->psi);
	}
      e->welfare3_t[TCH0][i] = pow(e->welfare3_t[TCH0][i],1.0/p->psi);
    }
}
#else
void calc_ce_ch_welfare1()
{
  int t, i, ib;

  stoch_eqm_ch *se = &(sss_ch[0]);
  const stoch_params_ch *sp = &(sppp_ch[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[NHIST_CH-1] );
  const eqm *e_stay = &(eee0_ch[0]);

  double C_stay, Cbb, C;
  double L_stay, Lbb, L;

  for(i=0; i<NC; i++)
    {
      t=NT;

      C_stay = e_stay->C_t[t][i];
      L_stay = e_stay->L_t[t][i];
      
      Cbb = (&(eee1_ch[NHIST_CH-2][0]))->C_t[t][i];
      Lbb = (&(eee1_ch[NHIST_CH-2][0]))->L_t[t][i];
      for(ib=NHIST_CH-3; ib>=0; ib--)
	{
	  Cbb = pi_brexit*Cbb + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->C_t[t][i];
	  Lbb = pi_brexit*Lbb + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->L_t[t][i];
	}

      C = pi_vote*C_stay + (1.0-pi_vote)*Cbb;
      L = pi_vote*L_stay + (1.0-pi_vote)*Lbb;
      
      e->welfare3_t[t][i] = (1.0/(1.0-p->beta)) * 
	(pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
      
      for(t=(NT-1); t>=TCH0; t--)
	{
	  C_stay = e_stay->C_t[t][i];
	  L_stay = e_stay->L_t[t][i];

	  Cbb = (&(eee1_ch[NHIST_CH-2][0]))->C_t[t][i];
	  Lbb = (&(eee1_ch[NHIST_CH-2][0]))->L_t[t][i];
	  for(ib=NHIST_CH-3; ib>=0; ib--)
	    {
	      Cbb = pi_brexit*Cbb + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->C_t[t][i];
	      Lbb = pi_brexit*Lbb + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->L_t[t][i];
	    }

	  C = pi_vote*C_stay + (1.0-pi_vote)*Cbb;
	  L = pi_vote*L_stay + (1.0-pi_vote)*Lbb;

	  e->welfare3_t[t][i] = p->beta * e->welfare3_t[t+1][i] + 
	    (pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
	  
	  e->welfare3_t[t+1][i] = pow(e->welfare3_t[t+1][i],1.0/p->psi);
	}
      e->welfare3_t[TCH0][i] = pow(e->welfare3_t[TCH0][i],1.0/p->psi);
    }
}
#endif

void calc_ce_welfare2()
{
  int t, i;

  stoch_eqm *se = &(sss[0]);
  const stoch_params *sp = &(sppp[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[1] );
  const eqm *e_soft = &(eee1[0]);
  const eqm * e_hard = &(eee2[0]);

  double C_soft, C_hard, C;
  double L_soft, L_hard, L;

  for(i=0; i<NC; i++)
    {
      t=NT;

      C_soft = e_soft->C_t[t][i];
      L_soft = e_soft->L_t[t][i];

      C_hard = e_hard->C_t[t][i];
      L_hard = e_hard->L_t[t][i];

      C = pi_brexit*C_soft + (1.0-pi_brexit)*C_hard;
      L = pi_brexit*L_soft + (1.0-pi_brexit)*L_hard;

      e->welfare4_t[t][i] = (1.0/(1.0-p->beta)) * 
	(pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
      
      for(t=(NT-1); t>=TVOTE; t--)
	{
	  C_soft = e_soft->C_t[t][i];
	  L_soft = e_soft->L_t[t][i];

	  C_hard = e_hard->C_t[t][i];
	  L_hard = e_hard->L_t[t][i];

	  C = pi_brexit*C_soft + (1.0-pi_brexit)*C_hard;
	  L = pi_brexit*L_soft + (1.0-pi_brexit)*L_hard;

	    e->welfare4_t[t][i] = p->beta * e->welfare4_t[t+1][i] + 
	      (pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
	      
	    e->welfare4_t[t+1][i] = pow(e->welfare4_t[t+1][i],1.0/p->psi);
	}
      e->welfare4_t[TVOTE][i] = pow(e->welfare4_t[TVOTE][i],1.0/p->psi);
    }
}

#ifdef CH_SIMPLE
void calc_ce_ch_welfare2()
{
  int t, i;

  stoch_eqm_ch *se = &(sss_ch[0]);
  const stoch_params_ch *sp = &(sppp_ch[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[2] );
  const eqm *e_temp = &(eee2[0]);
  const eqm * e_perm = &(eee1[0]);

  double C_temp, C_perm, C;
  double L_temp, L_perm, L;

  for(i=0; i<NC; i++)
    {
      t=NT;

      C_temp = e_temp->C_t[t][i];
      L_temp = e_temp->L_t[t][i];

      C_perm = e_perm->C_t[t][i];
      L_perm = e_perm->L_t[t][i];

      C = pi_brexit*C_temp + (1.0-pi_brexit)*C_perm;
      L = pi_brexit*L_temp + (1.0-pi_brexit)*L_perm;

      e->welfare4_t[t][i] = (1.0/(1.0-p->beta)) * 
	(pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
      
      for(t=(NT-1); t>=TCH1; t--)
	{
	  C_temp = e_temp->C_t[t][i];
	  L_temp = e_temp->L_t[t][i];

	  C_perm = e_perm->C_t[t][i];
	  L_perm = e_perm->L_t[t][i];

	  C = pi_brexit*C_temp + (1.0-pi_brexit)*C_perm;
	  L = pi_brexit*L_temp + (1.0-pi_brexit)*L_perm;

	    e->welfare4_t[t][i] = p->beta * e->welfare4_t[t+1][i] + 
	      (pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
	      
	    e->welfare4_t[t+1][i] = pow(e->welfare4_t[t+1][i],1.0/p->psi);
	}
      e->welfare4_t[TCH1][i] = pow(e->welfare4_t[TCH1][i],1.0/p->psi);
    }
}
#else
void calc_ce_ch_welfare2()
{
  int t, i, ib;

  stoch_eqm_ch *se = &(sss_ch[0]);
  const stoch_params_ch *sp = &(sppp_ch[0]);
  const params *p = &( (sp->ppp)[0] );

  eqm *e = &( (se->eee)[NHIST_CH-1] );

  double C;
  double L;

  for(i=0; i<NC; i++)
    {
      t=NT;

      C = (&(eee1_ch[NHIST_CH-2][0]))->C_t[t][i];
      L = (&(eee1_ch[NHIST_CH-2][0]))->L_t[t][i];
      for(ib=NHIST_CH-3; ib>=0; ib--)
	{
	  C = pi_brexit*C + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->C_t[t][i];
	  L = pi_brexit*L + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->L_t[t][i];
	}

      e->welfare4_t[t][i] = (1.0/(1.0-p->beta)) * 
	(pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
      
      for(t=(NT-1); t>=TCH1; t--)
	{
	  C = (&(eee1_ch[NHIST_CH-2][0]))->C_t[t][i];
	  L = (&(eee1_ch[NHIST_CH-2][0]))->L_t[t][i];
	  for(ib=NHIST_CH-3; ib>=0; ib--)
	    {
	      C = pi_brexit*C + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->C_t[t][i];
	      L = pi_brexit*L + (1.0-pi_brexit)*(&(eee1_ch[ib][0]))->L_t[t][i];
	    }

	    e->welfare4_t[t][i] = p->beta * e->welfare4_t[t+1][i] + 
	      (pow(C,p->phi[i]*p->psi) * pow(p->Lbar[i]-L,(1.0-p->phi[i])*p->psi));
	      
	    e->welfare4_t[t+1][i] = pow(e->welfare4_t[t+1][i],1.0/p->psi);
	}
      e->welfare4_t[TCH1][i] = pow(e->welfare4_t[TCH1][i],1.0/p->psi);
    }
}
#endif

///////////////////////////////

int bgp_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
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

int stoch_eqm_ch_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_stoch_eqm_ch_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int stoch_eqm_ch_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&stoch_eqm_ch_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int stoch_eqm_ch_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(stoch_eqm_ch_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&stoch_eqm_ch_func_f, x, J, 0))
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
