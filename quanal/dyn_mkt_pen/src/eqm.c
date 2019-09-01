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
  if(scenario != 1 && scenario !=6)
    {
      if(scenario==0)
	{
	  nn = NT+1;
	}
      else if(scenario==2 || scenario==3 || scenario>=7)
	{
	  nn = NT+1 - TREF;
	}
      else if(scenario==4 || scenario==5)
	{
	  nn = NT+1 - TVOTE;
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
  else if(scenario==6)
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

      // along branches 1 and 3, wait till TREV to find out whether Brexit is permanent or not
      nn = TREV-TBREXIT;
      neqm += nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;
      neqm += nn*NC*3 + nn*NC*NC + nn*(NC-1) + nn*NC + nn*NC;

      // after TREV, we have 4 possible histories that can arise, represented by 4 distinct
      // deterministic sections
      nn = (NT+1) - TREV;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;  
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;
      neqm += nn*NC*2 + (nn-1)*NC + nn*NC*NC + (nn-1)*(NC-1) + (nn-1)*NC + (nn-1)*NC;
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
	      e->ev_t[t][i][ii].Z=0.0;
	      e->ev_t[t][i][ii].psi_mult=1.0;
	      e->ev_t[t][i][ii].n=0.0;
	      e->ev_t[t][i][ii].export_participation_rate=0.0;
	      init_policy(&(e->ev_t[t][i][ii]));
	      //alloc_splines(&(e->ev_t[t][i][ii]));
	    }

	  if(idio_uncertainty && i<2)
	    {
	      e->ev_good_t[t][i].Z=0.0;
	      e->ev_good_t[t][i].psi_mult=1.0;
	      e->ev_good_t[t][i].n=0.0;
	      e->ev_good_t[t][i].export_participation_rate=0.0;
	      init_policy(&(e->ev_good_t[t][i]));

	      e->ev_bad_t[t][i].Z=0.0;
	      e->ev_bad_t[t][i].psi_mult=1.0;
	      e->ev_bad_t[t][i].n=0.0;
	      e->ev_bad_t[t][i].export_participation_rate=0.0;
	      init_policy(&(e->ev_bad_t[t][i]));
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
  memcpy((double *)(e1->texrate_t),(const double *)(e0->texrate_t),sizeof(double)*(NT+1)*NC*(NC-1));
  memcpy((double *)(e1->mkt_pen_t),(const double *)(e0->mkt_pen_t),sizeof(double)*(NT+1)*NC*(NC-1));
  memcpy((double *)(e1->mkt_pen_w_t),(const double *)(e0->mkt_pen_w_t),sizeof(double)*(NT+1)*NC*(NC-1));


  uint t, i, ii;
  for(t=0; t<(NT+1); t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(ii=0; ii<2; ii++)
	    {
	      copy_exporter_vars(&(e1->ev_t[t][i][ii]),(const exporter_vars *)(&(e0->ev_t[t][i][ii])));
	    }

	  if(i<2 && idio_uncertainty)
	    {
	      copy_exporter_vars(&(e1->ev_good_t[t][i]),(const exporter_vars *)(&(e0->ev_good_t[t][i])));
	      copy_exporter_vars(&(e1->ev_bad_t[t][i]),(const exporter_vars *)(&(e0->ev_bad_t[t][i])));
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
  if(scenario == 2 || scenario == 3 || scenario >=7)
    {
      nn = NT+1-TREF;
      t0 = TREF;
    }
  else if(scenario == 4 || scenario == 5)
    {
      nn = NT+1-TVOTE;
      t0 = TVOTE;
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
  if(scenario == 2 || scenario == 3 || scenario >=7)
    {
      nn = NT+1-TREF;
      t0 = TREF;
    }
  else if(scenario == 4 || scenario == 5)
    {
      nn = NT+1-TVOTE;
      t0 = TVOTE;
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

uint stack_stoch_eqm_rb_vars(double * myx, const stoch_eqm_rb * se)
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

  // now stack the TREV-TBREXIT periods associated with the "soft" branch
  t0 = TBREXIT;
  nn = TREV-TBREXIT;
  e = &((se->eee)[1]);
  
  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TREV; t++)
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

  // now stack the TREV-TBREXIT periods associated with the "hard" branch
  t0 = TBREXIT;
  nn = TREV-TBREXIT;
  e = &((se->eee)[2]);
  
  COPY_SUBVECTOR_LOG(myx+nx,&(e->W_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->Y_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->X_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->P2_t[t0]),(nn)*NC*NC);
  nx = nx + (nn)*NC*NC;

  for(t=t0+1; t<=TREV; t++)
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

  // now stack the post-Brexit part for the 4 permament/temporary hard/soft branches
  t0 = TREV;
  nn = NT+1-TREV;

  uint ih;
  for(ih=1; ih<NHIST_RB; ih++)
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
      fprintf(logfile,KRED "Error stacking stoch_eqm_rb vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}

uint unstack_stoch_eqm_rb_vars(stoch_eqm_rb * se, const double * myx)
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
  for(ih=NHIST_RB-1; ih>=0; ih--)
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

  // now do the TBREXIT-TVOTE period... here we need to copy branches 1-4 since they
  // must be identical here
  t0 = TVOTE;
  nn = TBREXIT-TVOTE;
  uint nx2=nx;
  for(ih=(NHIST_RB-1); ih>=1; ih--)
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

  // now do the TREV-TBREXIT period for soft... here we need to copy branches 1 and 3 since they
  // must be identical here
  t0 = TBREXIT;
  nn = TREV-TBREXIT;
  nx2=nx;
  for(ih=1; ih<NHIST_RB; ih+=2)
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

      for(t=t0+1; t<=TREV; t++)
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

  // now do the TREV-TBREXIT period for hard... here we need to copy branches 1 and 3 since they
  // must be identical here
  t0 = TBREXIT;
  nn = TREV-TBREXIT;
  nx2=nx;
  for(ih=2; ih<NHIST_RB; ih+=2)
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

      for(t=t0+1; t<=TREV; t++)
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

  // now stack the post-reveral part for all four Brexit branches
  t0 = TREV;
  nn = NT+1-TREV;

  for(ih=1; ih<NHIST_RB; ih++)
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
      fprintf(logfile,KRED "Error unstacking stoch_eqm_rb vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  return 0;
}

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

  if(allexport==1)
    {
      fname3 = concat(fname2,"_no_export_costs.txt");
    }
  else if(habit_formation==1)
    {
      fname3 = concat(fname2,"_habits.txt");
    }
  else if(idio_uncertainty)
    {
      fname3 = concat(fname2,"_idio_uncertainty.txt");
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==0)
    {
      fname3 = concat(fname2,"_static_arkolakis.txt");
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==1)
    {
      fname3 = concat(fname2,"_static_melitz.txt");
    }
  else if(no_exporter_dyn==0 && full_mkt_pen==1)
    {
      fname3 = concat(fname2,"_dynamic_sunk_cost.txt");
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
  else if(zeta_sens)
    {
      fname3 = concat(fname2,"_szeta.txt");
    }
  else if(psi_sens)
    {
      fname3 = concat(fname2,"_spsi.txt");
    }
  else if(low_exit_rate_sens)
    {
      fname3 = concat(fname2,"_sxr.txt");
    }
  else if(psi_uncertainty_sens==1)
    {
      fname3 = concat(fname2,"_kappa_uncertainty.txt");
    }
  else if(psi_uncertainty_sens==2)
    {
      fname3 = concat(fname2,"_kappa_ntb_uncertainty.txt");
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
      fprintf(file,"UK mkt. pen. : %0.3f\t\t%0.3f\n",100*e->mkt_pen_t[t][i][0],100*e->mkt_pen_t[t][i][1]);
      i=1;
      fprintf(file,"EU mkt. pen. : %0.3f\t\t%0.3f\n",100*e->mkt_pen_t[t][i][0],100*e->mkt_pen_t[t][i][1]);
      i=2;
      fprintf(file,"RW mkt. pen. : %0.3f\t\t%0.3f\n",100*e->mkt_pen_t[t][i][0],100*e->mkt_pen_t[t][i][1]);

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

  if(allexport==1)
    {
      fname3 = concat(fname2,"_no_export_costs.csv");
    }
  else if(habit_formation==1)
    {
      fname3 = concat(fname2,"_habits.csv");
    }
  else if(low_dep==1)
    {
      fname3 = concat(fname2,"_sld.csv");
    }
  else if(high_interest_rate==1)
    {
      fname3 = concat(fname2,"_shighr.csv");
    }
  else if(TBREXIT>8)
    {
      fname3 = concat(fname2,"_long_prebrexit.csv");
    }
  else if(idio_uncertainty)
    {
      fname3 = concat(fname2,"_idio_uncertainty.csv");
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==0)
    {
      fname3 = concat(fname2,"_static_arkolakis.csv");
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==1)
    {
      fname3 = concat(fname2,"_static_melitz.csv");
    }
  else if(no_exporter_dyn==0 && full_mkt_pen==1)
    {
      fname3 = concat(fname2,"_dynamic_sunk_cost.csv");
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
  else if(zeta_sens)
    {
      fname3 = concat(fname2,"_szeta.csv");
    }
  else if(psi_sens)
    {
      fname3 = concat(fname2,"_spsi.csv");
    }
  else if(low_exit_rate_sens)
    {
      fname3 = concat(fname2,"_sxr.csv");
    }
  else if(psi_uncertainty_sens==1)
    {
      fname3 = concat(fname2,"_kappa_uncertainty.csv");
    }
  else if(psi_uncertainty_sens==2)
    {
      fname3 = concat(fname2,"_kappa_ntb_uncertainty.csv");
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
      fprintf(file,"period,rgdp,ngdp,k,lp,rva,y,x,lf,c,W,sW,cW,ceW1,ceW2");
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
	  fprintf(file,",rer%d,tot%d,nx%d,ex%d,im%d,exrate%d,texrate%d,mktpen%d,mktpenw%d,Z%d,tau%d,ntb%d",j,j,j,j,j,j,j,j,j,j,j,j);
	}
      fprintf(file,"\n");

      for(t=0;t<(NT+1);t++)
	{
	  fprintf(file,"%d,",t);
	  fprintf(file,"%0.16f,",e->RGDP_t[t][i]);
	  fprintf(file,"%0.16f,",e->NGDP_t[t][i]);
	  fprintf(file,"%0.16f,",e->K_t[t][i]);
	  fprintf(file,"%0.16f,",e->L_t[t][i]-e->Lf_t[t][i]);
	  fprintf(file,"%0.16f,",pow(e->K_t[t][i],p->alpha)*pow(e->L_t[t][i]-e->Lf_t[t][i],1.0-p->alpha));
	  fprintf(file,"%0.16f,",e->Y_t[t][i]);
	  fprintf(file,"%0.16f,",e->P_t[t][i] * e->X_t[t][i]);
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
	      fprintf(file,"%0.16f,",e->EX_t[t][i][j]);
	      fprintf(file,"%0.16f,",e->IM_t[t][i][j]);
	      fprintf(file,"%0.16f,",e->exrate_t[t][i][ii]);
	      fprintf(file,"%0.16f,",e->texrate_t[t][i][ii]);
	      fprintf(file,"%0.16f,",e->mkt_pen_t[t][i][ii]);
	      fprintf(file,"%0.16f,",e->mkt_pen_w_t[t][i][ii]);
	      fprintf(file,"%0.16f,",e->ev_t[t][i][ii].Z);
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
  uint i,ii,j;

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

  for(i=0; i<NC; i++)
    {
      double lam = 0.0;
      if(t==NT)
	{
	  lam = p->beta;
	}
      else
	{
	  lam = e->Lambda_t[t][i];
	}

      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  set_exporter_consts1(&(e->ev_t[t][i][ii]),
			       e->W_t[t][i],
			       lam,
			       e->MC_t[t][i],
			       p->psi_mult_ts[t][j][i]);

	  set_exporter_consts2(&(e->ev_t[t][i][ii]),
			       p->L0[j]/p->L0[0],
			       e->P2_t[t][j][i],
			       e->Y2_t[t][j][i],
			       p->Ybar2[j][i],
			       p->theta,
			       p->ep[i][ii].alpha,
			       1.0+p->ntb_ts[t][j][i]);

	  set_exporter_consts3(&(e->ev_t[t][i][ii]),
			       p->L0[i]/p->L0[0],
			       e->P2_t[t][i][i],
			       e->Y2_t[t][i][i],
			       p->Ybar2[i][i],
			       p->theta);
	}

      if(idio_uncertainty && i<2)
	{
	  set_exporter_consts1(&(e->ev_t[t][i][0]),
			       e->W_t[t][i],
			       lam,
			       e->MC_t[t][i],
			       p->psi_mult_ts[t][1-i][i]);
	  set_exporter_consts2(&(e->ev_t[t][i][0]),
			       p->L0[j]/p->L0[0],
			       e->P2_t[t][1-i][i],
			       e->Y2_t[t][1-i][i],
			       p->Ybar2[1-i][i],
			       p->theta,
			       p->ep[i][0].alpha,
			       1.0 + 1.0*p->ntb_ts[t][1-i][i] - 0.0*p->tau_ts[t][1-i][i]);
	  set_exporter_consts3(&(e->ev_t[t][i][ii]),
			       p->L0[i]/p->L0[0],
			       e->P2_t[t][i][i],
			       e->Y2_t[t][i][i],
			       p->Ybar2[i][i],
			       p->theta);

	  set_exporter_consts1(&(e->ev_good_t[t][i]),
			       e->W_t[t][i],
			       lam,
			       e->MC_t[t][i],
			       p->psi_mult_ts[t][1-i][i]);
	  set_exporter_consts2(&(e->ev_good_t[t][i]),
			       p->L0[j]/p->L0[0],
			       e->P2_t[t][1-i][i],
			       e->Y2_t[t][1-i][i],
			       p->Ybar2[1-i][i],
			       p->theta,
			       p->ep[i][0].alpha,
			       1.0 + 0.0*p->ntb_ts[t][1-i][i] - p->tau_ts[t][1-i][i]);
	  set_exporter_consts3(&(e->ev_good_t[t][i]),
			       p->L0[i]/p->L0[0],
			       e->P2_t[t][i][i],
			       e->Y2_t[t][i][i],
			       p->Ybar2[i][i],
			       p->theta);

	  set_exporter_consts1(&(e->ev_bad_t[t][i]),
			       e->W_t[t][i],
			       lam,
			       e->MC_t[t][i],
			       p->psi_mult_ts[t][1-i][i]);
	  set_exporter_consts2(&(e->ev_bad_t[t][i]),
			       p->L0[j]/p->L0[0],
			       e->P2_t[t][1-i][i],
			       e->Y2_t[t][1-i][i],
			       p->Ybar2[1-i][i],
			       p->theta,
			       p->ep[i][0].alpha,
			       1.0 + p->ntb_ts[t][1-i][i]*3.0 + p->tau_ts[t][1-i][i]*2.0);
	  set_exporter_consts3(&(e->ev_bad_t[t][i]),
			       p->L0[i]/p->L0[0],
			       e->P2_t[t][i][i],
			       e->Y2_t[t][i][i],
			       p->Ybar2[i][i],
			       p->theta);

	}
    }

  return 0;
}

uint set_vars2a(eqm * e, const params * p, uint t, uint bgp)
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
	      if(no_exporter_dyn==0 && allexport==0)
		{
		  if(solve_steady_state_policies(&(p->ep[i][ii]), &(e->ev_t[t][i][ii])))
		    {
		      fprintf(logfile,
			      KRED "Failed to solve for steady state policy function %d/%d!\n" RESET,
			      i,j);
		      return 1;

		    }
		}
	      else
		{
		  fprintf(logfile,
			  KRED "Set vars2 called in static model %d/%d!\n" RESET,
			  i,j);
		  return 1;
		}
	    }
	  else
	    {
	      if(no_exporter_dyn==0 && allexport==0)
		{
		  if(full_mkt_pen==0)
		    {
		      double tmp=0.0;
		      if(iterate_entrant_policy_fn(&(p->ep[i][ii]),
						   &(e->ev_t[t][i][ii]),
						   &(e->ev_t[t+1][i][ii]),
						   &tmp))
			{
			  fprintf(logfile,
				  KRED "Failed to update entrant policies %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		      if(iterate_incumbent_policy_fn(&(p->ep[i][ii]),
						     &(e->ev_t[t][i][ii]),
						     &(e->ev_t[t+1][i][ii]),
						     &tmp))
			{
			  fprintf(logfile,
				  KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		    }
		  else
		    {
		      double tmp=0.0;
		      if(iterate_vf_policy_full_mkt_pen(&(p->ep[i][ii]),
							&(e->ev_t[t][i][ii]),
							&(e->ev_t[t+1][i][ii]),
							&tmp))
			{
			  fprintf(logfile,
				  KRED "Failed to update vf and policy (full mkt pen) %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		    }
		}
	      else
		{
		  fprintf(logfile,
			  KRED "Set vars2 called in static model %d/%d!\n" RESET,
			  i,j);
		  return 1;
		}
	    }
	}

      if(idio_uncertainty && i<2)
	{
	  if(t==NT)
	    {
	      if(solve_steady_state_policies(&(p->ep[i][0]), &(e->ev_good_t[t][i])))
		{
		  fprintf(logfile,
			  KRED "Failed to solve for steady state policy function %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	      if(solve_steady_state_policies(&(p->ep[i][0]), &(e->ev_bad_t[t][i])))
		{
		  fprintf(logfile,
			  KRED "Failed to solve for steady state policy function %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	    }
	  else
	    {
	      double tmp=0.0;
	      if(iterate_entrant_policy_fn(&(p->ep[i][0]),
					   &(e->ev_good_t[t][i]),
					   &(e->ev_good_t[t+1][i]),
					   &tmp))
		{
		  fprintf(logfile,
			  KRED "Failed to update entrant policies %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	      if(iterate_incumbent_policy_fn(&(p->ep[i][0]),
					     &(e->ev_good_t[t][i]),
					     &(e->ev_good_t[t+1][i]),
					     &tmp))
		{
		  fprintf(logfile,
			  KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}

	      if(iterate_entrant_policy_fn(&(p->ep[i][0]),
					   &(e->ev_bad_t[t][i]),
					   &(e->ev_bad_t[t+1][i]),
					   &tmp))
		{
		  fprintf(logfile,
			  KRED "Failed to update entrant policies %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	      if(iterate_incumbent_policy_fn(&(p->ep[i][0]),
					     &(e->ev_bad_t[t][i]),
					     &(e->ev_bad_t[t+1][i]),
					     &tmp))
		{
		  fprintf(logfile,
			  KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	    }
	}
    }

  return 0;
}

uint set_stoch_vars2a(eqm * e, eqm * e_good, eqm * e_bad, const params * p, uint t, double pi)
{
  uint i, ii, j;
  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  if(no_exporter_dyn==0 && allexport==0)
	    {
	      if(full_mkt_pen==0)
		{
		  if(iterate_entrant_policy_fn_stoch(&(p->ep[i][ii]),
						     &(e->ev_t[t][i][ii]),
						     &(e_good->ev_t[t+1][i][ii]),
						     &(e_bad->ev_t[t+1][i][ii]),
						     pi))
		    {
		      fprintf(logfile,
			      KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
			      i,j);
		      return 1;

		    }
		  if(iterate_incumbent_policy_fn_stoch(&(p->ep[i][ii]),
						       &(e->ev_t[t][i][ii]),
						       &(e_good->ev_t[t+1][i][ii]),
						       &(e_bad->ev_t[t+1][i][ii]),
						       pi))
		    {
		      fprintf(logfile,
			      KRED "Failed to update stochastic incumbent policies %d/%d!\n" RESET,
			      i,j);
		      return 1;

		    }
		}
	      else
		{
		  if(iterate_vf_policy_full_mkt_pen_stoch(&(p->ep[i][ii]),
							  &(e->ev_t[t][i][ii]),
							  &(e_good->ev_t[t+1][i][ii]),
							  &(e_bad->ev_t[t+1][i][ii]),
							  pi))
		    {
		      fprintf(logfile,
			      KRED "Failed to update stochastic vf and policy (full mkt pen) %d/%d!\n" RESET,
			      i,j);
		      return 1;

		    }
		}
	    }
	  else
	    {
	      fprintf(logfile,
		      KRED "Set vars2 called in static model %d/%d!\n" RESET,
		      i,j);
	      return 1;
	    }
	}
    }

  return 0;

}

uint set_stoch_vars2a_idio(eqm * e, eqm * e_good, eqm * e_bad, const params * p, uint t, double pi)
{
  uint i, ii, j;
  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  if(idio_uncertainty && i<2 && ii==0)
	    {
	      exporter_vars * evp_list[6] = {&(e_good->ev_good_t[t+1][i]),
					     &(e_good->ev_t[t+1][i][ii]),
					     &(e_good->ev_bad_t[t+1][i]),
					     &(e_good->ev_good_t[t+1][i]),
					     &(e_good->ev_t[t+1][i][ii]),
					     &(e_good->ev_bad_t[t+1][i])};
	      
	      double probs[6] = {pi/3.0,
				 pi/3.0,
				 pi/3.0,
				 (1.0-pi)/3.0,
				 (1.0-pi)/3.0,
				 (1.0-pi)/3.0};
	      
	      if(iterate_entrant_policy_fn_stoch_idio(&(p->ep[i][ii]),
						      &(e->ev_t[t][i][ii]),
						      evp_list,
						      probs))
	      {
		fprintf(logfile,
			KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
			i,j);
		return 1;
			  
	      }
	      if(iterate_incumbent_policy_fn_stoch_idio(&(p->ep[i][ii]),
							&(e->ev_t[t][i][ii]),
							evp_list,
							probs))
	      {
		fprintf(logfile,
			KRED "Failed to update stochastic incumbent policies %d/%d!\n" RESET,
			i,j);
		return 1;

	      }

	      if(iterate_entrant_policy_fn_stoch_idio(&(p->ep[i][ii]),
						      &(e->ev_good_t[t][i]),
						      evp_list,
						      probs))
	      {
		fprintf(logfile,
			KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
			i,j);
		return 1;

	      }
	      if(iterate_incumbent_policy_fn_stoch_idio(&(p->ep[i][ii]),
							&(e->ev_good_t[t][i]),
							evp_list,
							probs))
	      {
		fprintf(logfile,
			KRED "Failed to update stochastic incumbent policies %d/%d!\n" RESET,
			i,j);
		return 1;

	      }

	      if(iterate_entrant_policy_fn_stoch_idio(&(p->ep[i][ii]),
						      &(e->ev_bad_t[t][i]),
						      evp_list,
						      probs))
	      {
		fprintf(logfile,
			KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
			i,j);
		return 1;

	      }
	      if(iterate_incumbent_policy_fn_stoch_idio(&(p->ep[i][ii]),
							&(e->ev_bad_t[t][i]),
							evp_list,
							probs))
	      {
		fprintf(logfile,
			KRED "Failed to update stochastic incumbent policies %d/%d!\n" RESET,
			i,j);
		return 1;

	      }
	    }
	  else
	    {
	      if(no_exporter_dyn==0 && allexport==0)
		{
		  if(full_mkt_pen==0)
		    {
		      if(iterate_entrant_policy_fn_stoch(&(p->ep[i][ii]),
							 &(e->ev_t[t][i][ii]),
							 &(e_good->ev_t[t+1][i][ii]),
							 &(e_bad->ev_t[t+1][i][ii]),
							 pi))
			{
			  fprintf(logfile,
				  KRED "Failed to update stochastic entrant policies %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		      if(iterate_incumbent_policy_fn_stoch(&(p->ep[i][ii]),
							   &(e->ev_t[t][i][ii]),
							   &(e_good->ev_t[t+1][i][ii]),
							   &(e_bad->ev_t[t+1][i][ii]),
							   pi))
			{
			  fprintf(logfile,
				  KRED "Failed to update stochastic incumbent policies %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		    }
		  else
		    {
		      if(iterate_vf_policy_full_mkt_pen_stoch(&(p->ep[i][ii]),
							      &(e->ev_t[t][i][ii]),
							      &(e_good->ev_t[t+1][i][ii]),
							      &(e_bad->ev_t[t+1][i][ii]),
							      pi))
			{
			  fprintf(logfile,
				  KRED "Failed to update stochastic vf and policy (full mkt pen) %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		    }
		}
	      else
		{
		  fprintf(logfile,
			  KRED "Set vars2 called in static model %d/%d!\n" RESET,
			  i,j);
		  return 1;
		}
	    }
	}
    }

  return 0;

}

uint set_vars2b(eqm * e, const params * p, uint t, uint bgp)
{
  uint i, ii, j;

  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];
	  if(t==0)
	    {
	      if(full_mkt_pen==0)
		{
		  memcpy((double *)(e->ev_t[t][i][ii].diste),
			 (const double *)(p->diste0[i][ii]),
			 NPHI*sizeof(double));
		  memcpy((double *)(e->ev_t[t][i][ii].dist),
			 (const double *)(p->dist0[i][ii]),
			 NPHI*NN*sizeof(double));
		}
	      else
		{
		  memcpy((double *)(e->ev_t[t][i][ii].diste),
			 (const double *)(p->diste0[i][ii]),
			 NPHI*sizeof(double));
		  memcpy((double *)(e->ev_t[t][i][ii].disti),
			 (const double *)(p->disti0[i][ii]),
			 NPHI*sizeof(double));
		}
	    }
	  else if(t==NT && bgp)
	    {
	      if(no_exporter_dyn==0 && allexport==0)
		{
		  if(solve_steady_state_dist(&(p->ep[i][ii]), &(e->ev_t[t][i][ii])))
		    {
		      fprintf(logfile,
			      KRED "Failed to solve for steady state policy function %d/%d!\n" RESET,
			      i,j);
		      return 1;

		    }
		}
	      else
		{
		  fprintf(logfile,
			  KRED "Set vars2 called in static model %d/%d!\n" RESET,
			  i,j);
		  return 1;
		}
	    }
	  else
	    {
	      if(no_exporter_dyn==0 && allexport==0)
		{
		  if(full_mkt_pen==0)
		    {
		      double tmp_diste[NPHI];
		      double tmp_dist[NPHI][NN];
		  
		      double tmp=0.0;
		      if(update_dist(&(p->ep[i][ii]),
				     &(e->ev_t[t-1][i][ii]),
				     &(e->ev_t[t][i][ii]),
				     tmp_diste,
				     tmp_dist,
				     &tmp))
			{
			  fprintf(logfile,
				  KRED "Failed to update distribution %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		      memcpy(&(e->ev_t[t][i][ii].diste),tmp_diste,NPHI*sizeof(double));
		      memcpy(&(e->ev_t[t][i][ii].dist),tmp_dist,NPHI*NN*sizeof(double));
		    }
		  else
		    {
		      double tmp_diste[NPHI];
		      double tmp_disti[NPHI];
		  
		      double tmp=0.0;
		      if(update_dist_full_mkt_pen(&(p->ep[i][ii]),
						  &(e->ev_t[t-1][i][ii]),
						  &(e->ev_t[t][i][ii]),
						  tmp_diste,
						  tmp_disti,
						  &tmp))
			{
			  fprintf(logfile,
				  KRED "Failed to update distribution %d/%d!\n" RESET,
				  i,j);
			  return 1;

			}
		      memcpy(&(e->ev_t[t][i][ii].diste),tmp_diste,NPHI*sizeof(double));
		      memcpy(&(e->ev_t[t][i][ii].disti),tmp_disti,NPHI*sizeof(double));

		    }
		}
	      else
		{
		  fprintf(logfile,
			  KRED "Set vars2 called in static model %d/%d!\n" RESET,
			  i,j);
		  return 1;
		}
	    }
	}
    }
  
  for(i=0; i<NC; i++)
    {   
      if(idio_uncertainty && i<2)
	{
	  if(t==0)
	    {
	      memcpy((double *)(e->ev_good_t[t][i].diste),
		     (const double *)(p->diste0[i][0]),
		     NPHI*sizeof(double));
	      memcpy((double *)(e->ev_good_t[t][i].dist),
		     (const double *)(p->dist0[i][0]),
		     NPHI*NN*sizeof(double));

	      memcpy((double *)(e->ev_bad_t[t][i].diste),
		     (const double *)(p->diste0[i][0]),
		     NPHI*sizeof(double));
	      memcpy((double *)(e->ev_bad_t[t][i].dist),
		     (const double *)(p->dist0[i][0]),
		     NPHI*NN*sizeof(double));
	    }
	  else if(t==NT && bgp)
	    {
	      if(solve_steady_state_dist(&(p->ep[i][0]), &(e->ev_good_t[t][i])))
		{
		  fprintf(logfile,
			  KRED "Failed to solve for steady state policy function %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}

	      if(solve_steady_state_dist(&(p->ep[i][0]), &(e->ev_bad_t[t][i])))
		{
		  fprintf(logfile,
			  KRED "Failed to solve for steady state policy function %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	    }
	  else
	    {
	      double tmp_diste[NPHI];
	      double tmp_dist[NPHI][NN];
		  
	      double tmp=0.0;
	      if(update_dist(&(p->ep[i][0]),
			     &(e->ev_good_t[t-1][i]),
			     &(e->ev_good_t[t][i]),
			     tmp_diste,
			     tmp_dist,
			     &tmp))
		{
		  fprintf(logfile,
			  KRED "Failed to update distribution %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	      memcpy(&(e->ev_good_t[t][i].diste),tmp_diste,NPHI*sizeof(double));
	      memcpy(&(e->ev_good_t[t][i].dist),tmp_dist,NPHI*NN*sizeof(double));

	      if(update_dist(&(p->ep[i][0]),
			     &(e->ev_bad_t[t-1][i]),
			     &(e->ev_bad_t[t][i]),
			     tmp_diste,
			     tmp_dist,
			     &tmp))
		{
		  fprintf(logfile,
			  KRED "Failed to update distribution %d/%d!\n" RESET,
			  i,j);
		  return 1;

		}
	      memcpy(&(e->ev_bad_t[t][i].diste),tmp_diste,NPHI*sizeof(double));
	      memcpy(&(e->ev_bad_t[t][i].dist),tmp_dist,NPHI*NN*sizeof(double));
	    }
	}
    }
  
  return 0;
}

uint set_vars3(eqm * e, const params * p, uint t, uint bgp)
{

  uint i,j,ii;
  for(i=0; i<NC; i++)
    {      
      if(leontief)
	{
	  e->Kd_t[t][i] = (p->gamma[i]/e->MC_t[t][i]) *
	    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha-1.0) *
	    p->Nbar[i] * e->Dh_t[t][i][i] * p->ep[i][0].Zi;
	      
	  e->L_t[t][i] = (p->gamma[i]/e->MC_t[t][i]) *
	    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha) *
	    p->Nbar[i] * e->Dh_t[t][i][i] * p->ep[i][0].Zi;
	      
	  e->M_t[t][i] = ((1.0-p->gamma[i])/e->MC_t[t][i]) * 
	    p->Nbar[i] * e->Dh_t[t][i][i] * p->ep[i][0].Zi;
	}
      else
	{
	  e->Kd_t[t][i] = (p->gamma[i]*p->alpha/e->R_t[t][i]) * p->Nbar[i] * e->Dh_t[t][i][i] * p->ep[i][0].Zi;
	  e->L_t[t][i] = (p->gamma[i]*(1.0-p->alpha)/e->W_t[t][i]) * p->Nbar[i] * e->Dh_t[t][i][i] * p->ep[i][0].Zi;
	  e->M_t[t][i] = ((1.0-p->gamma[i])/e->P_t[t][i]) * p->Nbar[i] * e->Dh_t[t][i][i] * p->ep[i][0].Zi;
	}

      e->Y2s_t[t][i][i] =  p->Ybar2[i][i] * pow(e->MC_t[t][i] * (p->theta/(p->theta-1.0)),-p->theta) * 
	pow(p->Nbar[i],p->theta/(p->theta-1.0)) * e->D_t[t][i][i] * pow(p->ep[i][0].Zi,p->theta/(p->theta-1.0));

      e->VA_t[t][i] = e->Db_t[t][i][i] * p->Nbar[i] * p->ep[i][0].Zi;
      e->Pi_t[t][i] = e->Dt_t[t][i][i] * p->Nbar[i] * p->ep[i][0].Zi;
     
      e->Lf_t[t][i] = 0.0;

      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  if(no_exporter_dyn==0 && allexport==0)
	    {
	      if(full_mkt_pen==0)
		{
		  if(compute_moments(&(p->ep[i][ii]), &(e->ev_t[t][i][ii])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		      return 1;
		    }
		}
	      else
		{
		  if(compute_moments_full_mkt_pen(&(p->ep[i][ii]), &(e->ev_t[t][i][ii])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (full mkt pen) (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		      return 1;
		    }
		}
	    }
	  else
	    {
	      if(compute_moments_static_version(&(p->ep[i][ii]), &(e->ev_t[t][i][ii])))
		{
		  fprintf(logfile,"Error! Failed to compute steady-state moments! (static version) (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		  return 1;
		}
	    }
	}

      if(idio_uncertainty && i<2)
	{
	  if(no_exporter_dyn==0 && allexport==0)
	    {
	      if(full_mkt_pen==0)
		{
		  if(compute_moments(&(p->ep[i][0]), &(e->ev_good_t[t][i])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		      return 1;
		    }

		  if(compute_moments(&(p->ep[i][0]), &(e->ev_bad_t[t][i])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		      return 1;
		    }
		}
	      else
		{
		  if(compute_moments_full_mkt_pen(&(p->ep[i][0]), &(e->ev_good_t[t][i])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (full mkt pen) (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		      return 1;
		    }

		  if(compute_moments_full_mkt_pen(&(p->ep[i][0]), &(e->ev_bad_t[t][i])))
		    {
		      fprintf(logfile,"Error! Failed to compute steady-state moments! (full mkt pen) (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		      return 1;
		    }
		}
	    }
	  else
	    {
	      if(compute_moments_static_version(&(p->ep[i][0]), &(e->ev_good_t[t][i])))
		{
		  fprintf(logfile,"Error! Failed to compute steady-state moments! (static version) (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		  return 1;
		}

	      if(compute_moments_static_version(&(p->ep[i][0]), &(e->ev_bad_t[t][i])))
		{
		  fprintf(logfile,"Error! Failed to compute steady-state moments! (static version) (t,i,ii)=(%d,%d,%d)\n\n",t,i,ii);
		  return 1;
		}
	    }
	}
    }

  for(i=0; i<NC; i++)
    {
      for(ii=0; ii<2; ii++)
	{
	  j=p->Ji[i][ii];

	  if(i<2 && idio_uncertainty && ii==0)
	    {
	      e->exrate_t[t][i][ii] = (0.25*e->ev_t[t][i][ii].export_participation_rate+
				       0.5*e->ev_good_t[t][i].export_participation_rate+
				       0.25*e->ev_bad_t[t][i].export_participation_rate);

	      e->exrate_t[t][i][ii] = (0.25*e->ev_t[t][i][ii].theoretical_export_participation_rate+
				       0.5*e->ev_good_t[t][i].theoretical_export_participation_rate+
				       0.25*e->ev_bad_t[t][i].theoretical_export_participation_rate);

	      e->mkt_pen_t[t][i][ii] = (0.25*e->ev_t[t][i][ii].n+
					0.5*e->ev_good_t[t][i].n+
					0.25*e->ev_bad_t[t][i].n);

	      e->mkt_pen_w_t[t][i][ii] = (0.25*e->ev_t[t][i][ii].n2+
					0.5*e->ev_good_t[t][i].n2+
					0.25*e->ev_bad_t[t][i].n2);
	      
	      double Z = (0.25*e->ev_t[t][i][ii].Z + 0.5*e->ev_good_t[t][i].Z + 0.25*e->ev_bad_t[t][i].Z);

	      if(leontief)
		{
		  e->Kd_t[t][i] += (p->gamma[i]/e->MC_t[t][i]) *
		    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha-1.0) *
		    p->Nbar[i] * e->Dh_t[t][j][i] * Z;
	      
		  e->L_t[t][i] += (p->gamma[i]/e->MC_t[t][i]) *
		    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha) *
		    p->Nbar[i] * e->Dh_t[t][j][i] * Z;
	      
		  e->M_t[t][i] += ((1.0-p->gamma[i])/e->MC_t[t][i]) * 
		    p->Nbar[i] * e->Dh_t[t][j][i] * Z;
		}
	      else
		{
		  e->Kd_t[t][i] += (p->gamma[i]*p->alpha/e->R_t[t][i]) * 
		    p->Nbar[i] * e->Dh_t[t][j][i] * Z;

		  e->L_t[t][i] += (p->gamma[i]*(1.0-p->alpha)/e->W_t[t][i]) * 
		    p->Nbar[i] * e->Dh_t[t][j][i] * Z;

		  e->M_t[t][i] += ((1.0-p->gamma[i])/e->P_t[t][i]) * p->Nbar[i] * 
		    e->Dh_t[t][j][i] * Z;
		}

	      if(!allexport)
		{
		  e->Lf_t[t][i] += p->Nbar[i] * (0.25*e->ev_t[t][i][ii].lf+
						 0.5*e->ev_good_t[t][i].lf+
						 0.25*e->ev_bad_t[t][i].lf);
		}

	      e->Y2s_t[t][j][i] =  p->Ybar2[j][i] * pow(e->MC_t[t][i] * (p->theta/(p->theta-1.0)),-p->theta) * 
		pow(p->Nbar[i],p->theta/(p->theta-1.0)) * pow(1.0+p->ntb_ts[t][j][i],-p->theta) * e->D_t[t][j][i] * 
		pow(Z,p->theta/(p->theta-1.0));

	      e->Pi_t[t][i] += e->Dt_t[t][j][i] * p->Nbar[i] * Z;
	      e->T_t[t][i] = p->tau_ts[t][i][j] * e->IM_t[t][i][j];
	      e->VA_t[t][i] += e->EX_t[t][i][j];
	    }
	  else
	    {
	      e->exrate_t[t][i][ii] = e->ev_t[t][i][ii].export_participation_rate;
	      e->texrate_t[t][i][ii] = e->ev_t[t][i][ii].theoretical_export_participation_rate;
	      e->mkt_pen_t[t][i][ii] = e->ev_t[t][i][ii].n;
	      e->mkt_pen_w_t[t][i][ii] = e->ev_t[t][i][ii].n2;

	      if(leontief)
		{
		  e->Kd_t[t][i] += (p->gamma[i]/e->MC_t[t][i]) *
		    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha-1.0) *
		    p->Nbar[i] * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;
	      
		  e->L_t[t][i] += (p->gamma[i]/e->MC_t[t][i]) *
		    pow(e->R_t[t][i]*(1.0-p->alpha)/e->W_t[t][i]/p->alpha,p->alpha) *
		    p->Nbar[i] * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;
	      
		  e->M_t[t][i] += ((1.0-p->gamma[i])/e->MC_t[t][i]) * 
		    p->Nbar[i] * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;
		}
	      else
		{
		  e->Kd_t[t][i] += (p->gamma[i]*p->alpha/e->R_t[t][i]) * 
		    p->Nbar[i] * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;

		  e->L_t[t][i] += (p->gamma[i]*(1.0-p->alpha)/e->W_t[t][i]) * 
		    p->Nbar[i] * e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;

		  e->M_t[t][i] += ((1.0-p->gamma[i])/e->P_t[t][i]) * p->Nbar[i] * 
		    e->Dh_t[t][j][i] * e->ev_t[t][i][ii].Z;
		}

	      if(!allexport)
		{
		  e->Lf_t[t][i] += p->Nbar[i] * e->ev_t[t][i][ii].lf;
		}

	      e->Y2s_t[t][j][i] =  p->Ybar2[j][i] * pow(e->MC_t[t][i] * (p->theta/(p->theta-1.0)),-p->theta) * 
		pow(p->Nbar[i],p->theta/(p->theta-1.0)) * pow(1.0+p->ntb_ts[t][j][i],-p->theta) * e->D_t[t][j][i] * 
		pow(e->ev_t[t][i][ii].Z,p->theta/(p->theta-1.0));

	      e->Pi_t[t][i] += e->Dt_t[t][j][i] * p->Nbar[i] * e->ev_t[t][i][0].Z;
	      e->T_t[t][i] = p->tau_ts[t][i][j] * e->IM_t[t][i][j];
	      e->VA_t[t][i] += e->EX_t[t][i][j];
	    }
	}

      e->L_t[t][i] += e->Lf_t[t][i];
      e->NGDP_t[t][i] = e->P_t[t][i]*(e->Y_t[t][i] - e->M_t[t][i]);
      e->RGDP_t[t][i] = e->Y_t[t][i] - e->M_t[t][i];
      e->XY_t[t][i] = 100.0*e->P_t[t][i]*e->X_t[t][i]/e->NGDP_t[t][i];
      e->C_t[t][i] = e->Y_t[t][i] - e->X_t[t][i] - e->M_t[t][i];

      if(habit_formation==1)
	{
	  if(t==0)
	    {
	      e->Chabit_t[t][i] = e->C_t[t][i];
	    }
	  else
	    {
	      e->Chabit_t[t][i] = 0.0*e->Chabit_t[t-1][i] + 1.0*e->C_t[t-1][i];
	    }
	  e->MUC_t[t][i] = muc(e->C_t[t][i]-0.9*e->Chabit_t[t][i],
			       e->L_t[t][i],
			       p->Lbar[i],
			       p->phi[i],
			       p->psi);
	}
      else
	{
	  e->Chabit_t[t][i] = 0.0;
	  e->MUC_t[t][i] = muc(e->C_t[t][i],
			       e->L_t[t][i],
			       p->Lbar[i],
			       p->phi[i],
			       p->psi);
	}

      e->VA_t[t][i] -= e->P_t[t][i]*e->M_t[t][i];

      e->Inc_t[t][i] = e->W_t[t][i]*e->L_t[t][i] + e->R_t[t][i]*e->K_t[t][i]
	+ e->Pi_t[t][i] + e->T_t[t][i];     

      e->Exp_t[t][i] = e->P_t[t][i]*(e->C_t[t][i] + e->X_t[t][i])
	+ sum(e->NX_t[t][i],NC);
    }
  
  return 0;
}

uint eval_bgp_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,j=0,t=NT,nx=0;
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);

  e->B_t[t][0] = bbgp[0];
  e->B_t[t][1] = bbgp[1];
  unstack_bgp_vars(e,myx);
  if(set_vars1(e,p,t,1))
    {
      fprintf(logfile,KRED "Error calling set_vars1!\n" RESET);
      return 1;
    }
  if(allexport==0 && no_exporter_dyn==0 && set_vars2a(e,p,t,1))
    {
      fprintf(logfile,KRED "Error calling set_vars2a!\n" RESET);
      return 1;
    }
  if(allexport==0 && no_exporter_dyn==0 && set_vars2b(e,p,t,1))
    {
      fprintf(logfile,KRED "Error calling set_vars2b!\n" RESET);
      return 1;
    }
  if(set_vars3(e,p,t,1))
    {
      fprintf(logfile,KRED "Error calling set_vars3!\n" RESET);
      return 1;
    }
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
  uint i=0,ii=0, j=0,nx=0;
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
  else if(scenario == 7)
    {
      t0 = TREF;
      e = &(eee3[tn]);
      p = &(ppp3[tn]);
    }
  else if(scenario == 8)
    {
      t0 = TREF;
      e = &(eee4[tn]);
      p = &(ppp4[tn]);
    }

  unstack_eqm_vars(e,myx);

  e->B_t[0][0] = p->B0[0];
  e->B_t[0][1] = p->B0[1];
  e->B_t[0][2] = p->B0[2];

  if(scenario==0)
    {
      for(i=0; i<NC; i++)
	{
	  e->K_t[0][i] = p->K0[i];

	  for(ii=0; ii<NC-1; ii++)
	    {
	      if(no_exporter_dyn==0)
		{
		  if(full_mkt_pen==0)
		    {
		      memcpy((double *)(e->ev_t[0][i][ii].diste),
			     (const double *)(p->diste0[i][ii]),
			     NPHI*sizeof(double));
		      memcpy((double *)(e->ev_t[0][i][ii].dist),
			     (const double *)(p->dist0[i][ii]),
			     NPHI*NN*sizeof(double));
		    }
		  else
		    {
		      memcpy((double *)(e->ev_t[0][i][ii].diste),
			     (const double *)(p->diste0[i][ii]),
			     NPHI*sizeof(double));
		      memcpy((double *)(e->ev_t[0][i][ii].disti),
			     (const double *)(p->disti0[i][ii]),
			     NPHI*sizeof(double));
		    }
		}
	    }
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

  if(allexport==0 && no_exporter_dyn==0)
    {
      for(t=NT; t>=t0; t--)
	{
	  if(set_vars2a(e,p,t,0))
	    {
	      fprintf(logfile,KRED "Error calling set_vars2a!\n" RESET);
	      return 1;
	    }
	}
      for(t=t0; t<NT+1; t++)
	{
	  if(set_vars2b(e,p,t,0))
	    {
	      fprintf(logfile,KRED "Error calling set_vars2b!\n" RESET);
	      return 1;
	    }
	}
    }

  for(t=t0; t<(NT+1); t++)
    {
      if(set_vars3(e,p,t,0))
	{
	  fprintf(logfile,KRED "Error calling set_vars3!\n" RESET);
	  return 1;
	}
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
  
  if(no_exporter_dyn==0)
    {
      for(t=NT; t>=t0; t--)
	{
	  for(ih=0; ih<NHIST; ih++)
	    {
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );

	      if(t==(TVOTE-1))
		{
		  e0 = &(se->eee[0]);
		  e2 = &(se->eee[1]);
		  if(idio_uncertainty && set_stoch_vars2a_idio(e,e0,e2,p,t,pi_vote))
		    {
		      return 1;
		    }
		  else if(!idio_uncertainty && set_stoch_vars2a(e,e0,e2,p,t,pi_vote))
		    {
		      return 1;
		    }
		}
	      else if(e != &(se->eee[0]) && t==(TBREXIT-1))
		{
		  e0 = &(se->eee[1]);
		  e2 = &(se->eee[2]);
		  if(idio_uncertainty && set_stoch_vars2a_idio(e,e0,e2,p,t,pi_brexit))
		    {
		      return 1;
		    }
		  else if(!idio_uncertainty && set_stoch_vars2a(e,e0,e2,p,t,pi_brexit))
		    {
		      return 1;
		    }
		}
	      else
		{
		  if(set_vars2a(e,p,t,0))
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
	      set_vars2b(e,p,t,0);
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

uint eval_stoch_eqm_rb_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,j=0,t=NT,nx=0,ih=0;

  stoch_eqm_rb * se = &(sss_rb[tn]);
  stoch_params_rb * sp = &(sppp_rb[tn]);
  eqm *e, *e0, *e2;
  params *p, *p2;

  unstack_stoch_eqm_rb_vars(se,myx);

  uint t0 = TREF;
  for(ih=0; ih<NHIST_RB; ih++)
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
  
  if(no_exporter_dyn==0)
    {
      for(t=NT; t>=t0; t--)
	{
	  for(ih=0; ih<NHIST_RB; ih++)
	    {
	      e = &( (se->eee)[ih] );
	      p = &( (sp->ppp)[ih] );

	      if(t==(TVOTE-1))
		{
		  e0 = &(se->eee[0]);
		  e2 = &(se->eee[1]);
		  if(set_stoch_vars2a(e,e0,e2,p,t,pi_vote))
		    {
		      return 1;
		    }
		}
	      else if(e != &(se->eee[0]) && t==(TBREXIT-1))
		{
		  e0 = &(se->eee[1]);
		  e2 = &(se->eee[2]);
		  if(set_stoch_vars2a(e,e0,e2,p,t,pi_brexit))
		    {
		      return 1;
		    }
		}
	      else if( (e==&(se->eee[1]) || e==&(se->eee[3])) && t==(TREV-1))
		{
		  e0 = &(se->eee[1]);
		  e2 = &(se->eee[3]);
		  if(set_stoch_vars2a(e,e0,e2,p,t,pi_rev))
		    {
		      return 1;
		    }
		}
	      else if( (e==&(se->eee[2]) || e==&(se->eee[4])) && t==(TREV-1))
		{
		  e0 = &(se->eee[2]);
		  e2 = &(se->eee[4]);
		  if(set_stoch_vars2a(e,e0,e2,p,t,pi_rev))
		    {
		      return 1;
		    }
		}
	      else
		{
		  if(set_vars2a(e,p,t,0))
		    {
		      fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
		      return 1;
		    }
		}
	    }
	}

      for(ih=0; ih<NHIST_RB; ih++)
	{
	  e = &( (se->eee)[ih] );
	  p = &( (sp->ppp)[ih] );

	  for(t=t0; t<(NT+1); t++)
	    {
	      set_vars2b(e,p,t,0);
	    }
	}
    }
      
  for(ih=0; ih<NHIST_RB; ih++)
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

  // now we do the TREV-TBREXIT periods between Brexit and reversal for soft Brexit...
  // here there is uncertainty again in the period before reversal only
  t0 = TBREXIT;
  e = &(se->eee[1]);
  p = &(sp->ppp[1]);
  for(t=t0; t<(TREV-1); t++)
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

  // in the period before the reversal, we have to deal with the uncertainty about whether
  // Brexit is temporary or permanent
  t = TREV-1;
  e = &(se->eee[1]);
  p = &(sp->ppp[1]);
  e2 = &(se->eee[3]);
  p2 = &(sp->ppp[3]);

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

      myf[nx] = pi_rev * euler(p,e,t,i) + (1.0-pi_rev) * euler(p2,e2,t,i);
      nx = nx+1;

      myf[nx] = pi_rev * noarb(p,e,t,i) + (1.0-pi_rev) * noarb(p2,e2,t,i);
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

  // now we do the TREV-TBREXIT periods between Brexit and reversal for hard Brexit...
  // here there is uncertainty again in the period before reversal only
  t0 = TBREXIT;
  e = &(se->eee[2]);
  p = &(sp->ppp[2]);
  for(t=t0; t<(TREV-1); t++)
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

  // in the period before the reversal, we have to deal with the uncertainty about whether
  // Brexit is temporary or permanent
  t = TREV-1;
  e = &(se->eee[2]);
  p = &(sp->ppp[2]);
  e2 = &(se->eee[4]);
  p2 = &(sp->ppp[4]);

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

      myf[nx] = pi_rev * euler(p,e,t,i) + (1.0-pi_rev) * euler(p2,e2,t,i);
      nx = nx+1;

      myf[nx] = pi_rev * noarb(p,e,t,i) + (1.0-pi_rev) * noarb(p2,e2,t,i);
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

  // after reversal, we have to worry about the separate possible histories
  // but there is no uncertainty to worry about... everything is deterministic
  // within a given history from reversal onwards
  t0 = TREV;
  for(ih=1; ih<NHIST_RB; ih++)
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

  if(allexport)
    {
      sname = "output/seed_no_export_costs.bin";
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==0)
    {
      sname = "output/seed_static_arkolakis.bin";
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==1)
    {
      sname = "output/seed_static_melitz.bin";
    }
  else if(no_exporter_dyn==0 && full_mkt_pen==1)
    {
      sname = "output/seed_dynamic_sunk_cost.bin";
    }
  else if(zeta_sens)
    {
      sname = "output/seed_szeta.bin";
    }
  else
    {
      sname = "output/seed_baseline.bin";
    }

  // if we are solving for the no-Brexit (or no-reform in China-like exercise) counterfactual we must construct an initial guess for the solver
  if(scenario==0)
    {
      if(read_seed==1)
	{
	  free_solver_mem();
	  solver_n = neqm;
	  alloc_solver_mem();

	  fprintf(logfile,KMAG "\n\tReading from seed file %s\n" RESET,sname);

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
      else if(scenario == 7)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee3[it]), &(eee0[0])  );
	    }
	}
      else if(scenario == 8)
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee4[it]), &(eee0[0])  );
	    }
	}

      // now determine which eqm struct to use as the initial guess...
      // ...if we are in the optimistic scenario, use the first path...
      if(scenario == 2 || scenario == 4 || scenario==7)
	{
	  e = &(eee0[0]);
	}
      // ...otherwise use the second path
      else if(scenario == 3 || scenario == 5 || scenario==8)
	{
	  e = &(eee0[0]);
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

      if((scenario==0) && write_seed==1 && !status)
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
  if(allexport==1)
    {
      sname = "output/stoch_seed_no_export_costs.bin";
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==0)
    {
      sname = "output/stoch_seed_static_arkolakis.bin";
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==1)
    {
      sname = "output/stoch_seed_static_melitz.bin";
    }
  else if(no_exporter_dyn==0 && full_mkt_pen==1)
    {
      sname = "output/stoch_seed_dynamic_sunk_cost.bin";
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
  else if(zeta_sens)
    {
      sname = "output/stoch_seed_szeta.bin";
    }
  else if(psi_sens)
    {
      sname = "output/stoch_seed_spsi.bin";
    }
  else if(psi_uncertainty_sens)
    {
      sname = "output/stoch_seed_kappa_uncertainty.bin";
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

uint solve_stoch_eqm_rb()
{  
  char * sname;
  if(allexport==1)
    {
      sname = "output/stoch_seed_rb_no_export_costs.bin";
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==0)
    {
      sname = "output/stoch_seed_rb_static_arkolakis.bin";
    }
  else if(no_exporter_dyn==1 && full_mkt_pen==1)
    {
      sname = "output/stoch_seed_rb_static_melitz.bin";
    }
  else if(no_exporter_dyn==0 && full_mkt_pen==1)
    {
      sname = "output/stoch_seed_rb_dynamic_sunk_cost.bin";
    }
  else if(low_pi)
    {
      sname = "output/stoch_rb_seed_lowpi.bin";
    }
  else if(high_pi)
    {
      sname = "output/stoch_rb_seed_highpi.bin";
    }
  else if(fin_aut)
    {
      sname = "output/stoch_rb_seed_finaut.bin";
    }
  else if(zeta_sens)
    {
      sname = "output/stoch_rb_seed_szeta.bin";
    }
  else if(psi_sens)
    {
      sname = "output/stoch_rb_seed_spsi.bin";
    }
  else if(psi_uncertainty_sens)
    {
      sname = "output/stoch_rb_seed_kappa_uncertainty.bin";
    }
  else
    {
      sname = "output/stoch_rb_seed_baseline.bin";
    }

  // before we do anything, we need to copy all variables from the no-Brexit counterfactual to every 
  // thread-specific copy of every branch of the stochastic equilibrium structure...
  // i.e. copy everything from eee[0] to the stochastic equilibrium sss[it] for each
  // thread it (we have to do it for all thread because otherwise those threads' variables for
  // t = 0 through t = TREF will never get set
  uint it,ih;
  for(it=0; it<NTH; it++)
    {
      for(ih=0; ih<NHIST_RB; ih++)
	{
	  copy_vars( &( (&(sss_rb[it]))->eee[ih] ) , &(eee0[0]) );
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
  else
    {
      copy_vars( &( (&(sss_rb[0]))->eee[0] ) , &(eee0[0]) );
      copy_vars( &( (&(sss_rb[0]))->eee[1] ) , &(eee1[0]) );
      copy_vars( &( (&(sss_rb[0]))->eee[2] ) , &(eee2[0]) );
      copy_vars( &( (&(sss_rb[0]))->eee[3] ) , &(eee3[0]) );
      copy_vars( &( (&(sss_rb[0]))->eee[4] ) , &(eee4[0]) );

      if(stack_stoch_eqm_rb_vars(solver_x->data,&(sss_rb[0])))
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
      status = stoch_eqm_rb_func_f(solver_x,NULL,f0[0]);
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      par=1;
      gsl_multiroot_function_fdf f = {&stoch_eqm_rb_func_f,
				      &stoch_eqm_rb_func_df,
				      &stoch_eqm_rb_func_fdf,
				      neqm,
				      NULL};
      status = find_root_deriv_mkl(&f);
      if(status)
	fprintf(logfile,KRED "Error solving for stochastic equilibrium with reversal!\n" RESET);

      if(!status && write_stoch_seed==1)
	{
	  write_vec_bin(solver_x->data, neqm, sname);
	}
    }

  free_solver_mem();

  return status;
}

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

/*void calc_stoch_welfare()
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
    }*/
    
/*void calc_ce_welfare1()
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
    }*/

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

int stoch_eqm_rb_func_f(const gsl_vector * x, void * data, gsl_vector * f)
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
  if(eval_stoch_eqm_rb_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int stoch_eqm_rb_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&stoch_eqm_rb_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int stoch_eqm_rb_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(stoch_eqm_rb_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&stoch_eqm_rb_func_f, x, J, 0))
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
