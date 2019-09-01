#ifndef __MAIN_C__
#define __MAIN_C__

#include "globals.h"
#include "calibrate.h"
#include "eqm.h"

uint parse_args(int argc, char **argv)
{
  int opt = 0;
  int cnt=0;
  int cnt2=0;

  recal=1;

  while((opt = getopt(argc, argv, "rsflhpzdmnkbixquv")) != -1){
    switch(opt){
    case 'v':
      low_dep=1;
      break;
    case 'u':
      high_interest_rate=1;
      break;
    case 'q':
      habit_formation=1;
      cnt++;
      break;
    case 'r':
      reversible_brexit=1;
      cnt++;
      break;
    case 'i':
      idio_uncertainty=1;
      cnt++;
      break;
    case 'x':
      low_exit_rate_sens=1;
      cnt++;
      break;
    case 'l':
      low_pi = 1;
      high_pi = 0;
      cnt++;
      break;
    case 'h':
      low_pi = 0;
      high_pi = 1;
      cnt++;
      break;
    case 'd':
      no_exporter_dyn=1;
      recal=1;
      cnt2++;
      break;
    case 'm':
      full_mkt_pen=1;
      recal=1;
      cnt2++;
      break;
    case 'n':
      allexport=1;
      no_exporter_dyn=1;
      recal=1;
      cnt++;
      break;
    case 'f':
      fin_aut = 1;
      cnt++;
      break;
    case 's':
      fixl = 2;
      cnt++;
      break;
    case 'p':
      psi_sens=1;
      recal=1;
      cnt++;
      break;
    case 'z':
      zeta_sens=1;
      recal=1;
      cnt++;
      break;
    case 'k':
      psi_uncertainty_sens=1;
      recal=1;
      cnt++;
      break;
    case 'b':
      psi_uncertainty_sens=2;
      recal=1;
      cnt++;
      break;
    case '?':
      fprintf(logfile,
	      "\nIncorrect command line option: %c. Possible options: -f -h -l -s -p -z -c -r -w -d -m -n -k -b\n",opt);
      fprintf(logfile,"\t-p: higher risk aversion\n");
      fprintf(logfile,"\t-d: no exporter dynamics\n");
      fprintf(logfile,"\t-m: full market penetration\n");
      fprintf(logfile,"\t-n: no exporting costs\n");
      fprintf(logfile,"\t-k: uncertainty in marketing costs instead of non-tariff barriers\n");
      fprintf(logfile,"\t-b: uncertainty in marketing costs and non-tariff barriers\n");
      fprintf(logfile,"\t-p: higher risk aversion\n");
      fprintf(logfile,"\t-d: no exporter dynamics\n");
      fprintf(logfile,"\t-m: full market penetration\n");
      fprintf(logfile,"\t-n: no exporting costs\n");
      fprintf(logfile,"\t-z: lower Armington elasticity\n");
      fprintf(logfile,"\t-r: reversible Brexit\n");
      fprintf(logfile,"\t-i: idiosyncratic uncertainty\n");
      fprintf(logfile,"\t-x: lower exit rate\n");

      return 1;
      break;
    }
  }

  if(cnt>1 || (cnt>0 && cnt2>0))
    {
      fprintf(logfile,KRED "Only one sensitivity analysis allowed (except -d + -m)!\n" RESET);
      return 1;
    }

  if(allexport)
    {
      fprintf(logfile,KBLU "\tNo exporting costs: kappa=0\n" RESET);
    }
  else if(full_mkt_pen==0 && no_exporter_dyn==0)
    {
      fprintf(logfile,KBLU "\tBaseline model with market penetation dynamics: kappa(n',n)\n" RESET);
    }
  else if(full_mkt_pen==0 && no_exporter_dyn==1)
    {
      fprintf(logfile,KBLU "\tArkolakis (2010) model with static market penetration decisions: kappa(n')\n" RESET);
    }
  else if(full_mkt_pen==1 && no_exporter_dyn==0)
    {
      fprintf(logfile,KBLU "\tDynamic export participation decisions with full market penetration: kappa(s)\n" RESET);
    }
  else if(full_mkt_pen==1 && no_exporter_dyn==1)
    {
      fprintf(logfile,KBLU "\tStatic export participation + full market penetration (Melitz, 2003): kappa\n" RESET);
    }

  if(psi_uncertainty_sens==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with increases in marketing costs instead of NTBs\n" RESET);
    }
  if(psi_uncertainty_sens==2)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with increases in marketing costs and NTBs\n" RESET);
    }
  if(psi_sens==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with higher risk aversion\n" RESET);
    }
  if(zeta_sens==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with lower Armington elasticity\n" RESET);
    }
  if(low_pi==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with lower probability of soft Brexit\n" RESET);
    }
  if(high_pi==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with lower probability of hard Brexit\n" RESET);
    }
  if(fin_aut==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with financial autarky\n" RESET);
    }
  if(fixl==2)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with sticky wages\n" RESET);
    }
  if(low_exit_rate_sens==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with lower exit rate\n" RESET);
    }
  if(reversible_brexit==1)
    {
      fprintf(logfile,KBLU "\tAdditional TPU: reversible Brexit\n" RESET);
    }
  if(idio_uncertainty==1)
    {
      fprintf(logfile,KBLU "\tAdditional TPU: firm-level tariff uncertainty\n" RESET);
    }
  if(habit_formation==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with habit formation\n" RESET);
    }
  if(high_interest_rate==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with higher interest rate\n" RESET);
    }
    if(TBREXIT>8)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with longer pre-Brexit period\n" RESET);
    }
    if(low_dep==1)
    {
      fprintf(logfile,KBLU "\tSensitivity analysis with lower customer depreciation\n" RESET);
    }
  
  return 0;
}

int quant_exercise()
{
  // ---------------------------------------------------------------------------------------------------
  // set up variable and parameter structures
  uint it;
  eqm * e;
  params * p;
  
  BREAK2;
  fprintf(logfile, KGRN "\nCalibrating model...\n");

  epsjac=1.0e-7;
  root_tol=1.0e-7;
  calc_all_moments=1;
  if(calibrate())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  epsjac=1.0e-8;
  root_tol=7.5e-8;
  calc_all_moments=0;

  // no-Brexit counterfactual
  BREAK2;
  scenario = 0;
  set_neqm();
  fprintf(logfile,KGRN "\nSolving for counterfactual no-Brexit equilibrium...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  calc_welfare(&(eee0[0]), &(ppp0[0]));
  write_eqm_vars(&(ppp0[0]),&(eee0[0]),"vars_nobrexit_uk",0);
  write_eqm_vars(&(ppp0[0]),&(eee0[0]),"vars_nobrexit_eu",1);
  write_eqm_vars(&(ppp0[0]),&(eee0[0]),"vars_nobrexit_rw",2);
  write_bgp_vars(&(eee0[0]),"bgp_nobrexit");

  // perfect foresight (soft Brexit, TREF onwards) 
  BREAK2;
  scenario = 2;
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp1[it]),scenario);
    }
  set_neqm();
  fprintf(logfile,KGRN "\nSolving for perfect-foresight soft Brexit equilibrium (referendum onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  calc_welfare(&(eee1[0]), &(ppp1[0]));
  write_eqm_vars(&(ppp1[0]),&(eee1[0]),"vars_det_soft_uk",0);
  write_eqm_vars(&(ppp1[0]),&(eee1[0]),"vars_det_soft_eu",1);
  write_eqm_vars(&(ppp1[0]),&(eee1[0]),"vars_det_soft_rw",2);
  write_bgp_vars(&(eee1[0]),"bgp_det_soft");

  // perfect foresight (hard Brexit, TREF onwards) 
  BREAK2;
  scenario = 3;
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp2[it]),scenario);
    }
  set_neqm();
  fprintf(logfile,KGRN "\nSolving for perfect-foresight hard Brexit equilibrium (referendum onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  calc_welfare(&(eee2[0]), &(ppp2[0]));
  write_eqm_vars(&(ppp2[0]),&(eee2[0]),"vars_det_hard_uk",0);
  write_eqm_vars(&(ppp2[0]),&(eee2[0]),"vars_det_hard_eu",1);
  write_eqm_vars(&(ppp2[0]),&(eee2[0]),"vars_det_hard_rw",2);
  write_bgp_vars(&(eee2[0]),"bgp_det_hard");

  if(reversible_brexit==0)
    {
      stoch_eqm * se;
      stoch_params * sp;

      // stochastic model
      BREAK2;
      scenario = 1;
      set_neqm();
      for(it=0; it<NTH; it++)
	{
	  set_tariffs(&(sppp[it]));
	}

      fprintf(logfile,KGRN "\nSolving for stochastic model with uncertainty about trade costs after Brexit...\n" RESET);
      if(solve_stoch_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      se = &(sss[0]);
      sp = &(sppp[0]);
 
      e = &((se->eee)[0]);
      p = &((sp->ppp)[0]);
      calc_welfare(e,p);

      e = &((se->eee)[1]);
      p = &((sp->ppp)[1]);
      calc_welfare(e,p);
  
      e = &((se->eee)[2]);
      p = &((sp->ppp)[2]);
      calc_welfare(e,p);

      e = &((se->eee)[0]);
      p = &((sp->ppp)[0]); 
      write_eqm_vars(p,e,"vars_stoch_stay_uk",0);
      write_eqm_vars(p,e,"vars_stoch_stay_eu",1);
      write_eqm_vars(p,e,"vars_stoch_stay_rw",2);
 
      e = &((se->eee)[1]);
      p = &((sp->ppp)[1]);
      write_eqm_vars(p,e,"vars_stoch_soft_uk",0);
      write_eqm_vars(p,e,"vars_stoch_soft_eu",1);
      write_eqm_vars(p,e,"vars_stoch_soft_rw",2);
  
      e = &((se->eee)[2]);
      p = &((sp->ppp)[2]);
      write_eqm_vars(p,e,"vars_stoch_hard_uk",0);
      write_eqm_vars(p,e,"vars_stoch_hard_eu",1);
      write_eqm_vars(p,e,"vars_stoch_hard_rw",2);
    }
  else
    {
      stoch_eqm_rb * se;
      stoch_params_rb * sp;

      // perfect foresight (temporary soft Brexit, TVOTE onwards) 
      BREAK2;
      scenario = 7;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp3[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KGRN "\nSolving for perfect-foresight version of TEMPORARY soft Brexit...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      calc_welfare(&(eee3[0]), &(ppp3[0]));
      write_eqm_vars(&(ppp3[0]),&(eee3[0]),"vars_det_temp_soft_uk",0);
      write_eqm_vars(&(ppp3[0]),&(eee3[0]),"vars_det_temp_soft_eu",1);
      write_eqm_vars(&(ppp3[0]),&(eee3[0]),"vars_det_temp_soft_rw",2);
      write_bgp_vars(&(eee3[0]),"bgp_det_temp_soft");
      
      // perfect foresight (temporary hard Brexit, TVOTE onwards) 
      BREAK2;
      scenario = 8;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp4[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KGRN "\nSolving for perfect-foresight version of TEMPORARY hard Brexit...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      calc_welfare(&(eee4[0]), &(ppp4[0]));
      write_eqm_vars(&(ppp4[0]),&(eee4[0]),"vars_det_temp_hard_uk",0);
      write_eqm_vars(&(ppp4[0]),&(eee4[0]),"vars_det_temp_hard_eu",1);
      write_eqm_vars(&(ppp4[0]),&(eee4[0]),"vars_det_temp_hard_rw",2);
      write_bgp_vars(&(eee4[0]),"bgp_det_temp_hard");

      // stochastic model
      BREAK2;
      scenario = 6;
      set_neqm();
      for(it=0; it<NTH; it++)
	{
	  set_tariffs_rb(&(sppp_rb[it]));
	}

      fprintf(logfile,KGRN "\nSolving for stochastic model with uncertainty about trade costs after Brexit AND about Brexit permanence...\n" RESET);
      if(solve_stoch_eqm_rb())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      se = &(sss_rb[0]);
      sp = &(sppp_rb[0]);
 
      e = &((se->eee)[0]);
      p = &((sp->ppp)[0]);
      calc_welfare(e,p);

      e = &((se->eee)[1]);
      p = &((sp->ppp)[1]);
      calc_welfare(e,p);
  
      e = &((se->eee)[2]);
      p = &((sp->ppp)[2]);
      calc_welfare(e,p);

      e = &((se->eee)[3]);
      p = &((sp->ppp)[3]);
      calc_welfare(e,p);

      e = &((se->eee)[4]);
      p = &((sp->ppp)[4]);
      calc_welfare(e,p);

      e = &((se->eee)[0]);
      p = &((sp->ppp)[0]); 
      write_eqm_vars(p,e,"vars_stoch_rb_stay_uk",0);
      write_eqm_vars(p,e,"vars_stoch_rb_stay_eu",1);
      write_eqm_vars(p,e,"vars_stoch_rb_stay_rw",2);
 
      e = &((se->eee)[1]);
      p = &((sp->ppp)[1]);
      write_eqm_vars(p,e,"vars_stoch_rb_perm_soft_uk",0);
      write_eqm_vars(p,e,"vars_stoch_rb_perm_soft_eu",1);
      write_eqm_vars(p,e,"vars_stoch_rb_perm_soft_rw",2);
  
      e = &((se->eee)[2]);
      p = &((sp->ppp)[2]);
      write_eqm_vars(p,e,"vars_stoch_rb_perm_hard_uk",0);
      write_eqm_vars(p,e,"vars_stoch_rb_perm_hard_eu",1);
      write_eqm_vars(p,e,"vars_stoch_rb_perm_hard_rw",2);

      e = &((se->eee)[3]);
      p = &((sp->ppp)[3]);
      write_eqm_vars(p,e,"vars_stoch_rb_temp_soft_uk",0);
      write_eqm_vars(p,e,"vars_stoch_rb_temp_soft_eu",1);
      write_eqm_vars(p,e,"vars_stoch_rb_temp_soft_rw",2);
  
      e = &((se->eee)[4]);
      p = &((sp->ppp)[4]);
      write_eqm_vars(p,e,"vars_stoch_rb_temp_hard_uk",0);
      write_eqm_vars(p,e,"vars_stoch_rb_temp_hard_eu",1);
      write_eqm_vars(p,e,"vars_stoch_rb_temp_hard_rw",2);
    } 

  return 0;

}

int main(int argc, char * argv[])
{
  int it;
  time_t start, stop;
  time(&start);

  psi_uncertainty_sens=0;
  par=0;
  psi_sens=0;
  zeta_sens=0;
  no_exporter_dyn=0;
  full_mkt_pen=0;
  leontief = 1;
  low_pi = 0;
  high_pi = 0;
  fin_aut = 0;
  recal = 0;
  eval_calfn_once=0;
  eval_eqm_once_flag=0;
  read_seed=0;
  write_seed=0;
  read_stoch_seed=0;
  write_stoch_seed=0;
  fixl=1;
  tfp_flag=0;
  logfile = stdout;
  small_gnewton_step_flag=0;
  reversible_brexit=0;
  idio_uncertainty=0;
  low_exit_rate_sens=0;
  high_interest_rate=0;
  
  BREAK2;
  fprintf(logfile, KGRN "\nBrexit and the Macroeconomic Impact of Trade Policy Uncertainty" RESET);
  fprintf(logfile, KGRN "\nJoseph Steinberg, University of Toronto" RESET);
  fprintf(logfile, KGRN "\nMarket penetration dynamics version" RESET);
  fprintf(logfile, KGRN "\nLast updated: August 2018" RESET);
  fprintf(logfile, KNRM "\n" RESET);

  BREAK2;
  fprintf(logfile, KGRN "\nSetting up environment...\n\n" RESET);

  if(parse_args(argc,argv))
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    };
  
  // set up parallel environment
#ifdef _OPENMP
  omp_set_num_threads(NTH);
  mkl_set_num_threads(NTH);
  uint nt = omp_get_max_threads();
  fprintf(logfile, KBLU "\n\tParallel processing using %d OMP threads\n" RESET,nt);
#pragma omp parallel num_threads(nt)
  {
    int it = omp_get_thread_num();
    fprintf(logfile,KBLU "\t\tHello from thread %d out of %d\n" RESET,
	    it, nt);
  }
  fprintf(logfile,"\n");
#endif

  fprintf(logfile, KBLU "\tInitializing equilibrium structure...\n");
  
  // initiate equilibrium structure
#pragma omp parallel for private(it) schedule(static)
  for(it=0; it<NTH; it++)
    {
      init_vars(&(eee0[it]));
      init_vars(&(eee1[it]));
      init_vars(&(eee2[it]));
      init_vars(&(eee3[it]));
      init_vars(&(eee4[it]));
      stoch_eqm * se = &(sss[it]);
      int ih;
      for(ih=0; ih<NHIST; ih++)
	{
	  init_vars( &((se->eee)[ih]) );
	}
      stoch_eqm_rb * se_rb = &(sss_rb[it]);
      for(ih=0; ih<NHIST_RB; ih++)
	{
	  init_vars( &((se_rb->eee)[ih]) );
	}
    }

  solver_verbose=1;

  pi_vote=0.75;
  pi_rev=0.5;
  if(low_pi==1)
    {
      pi_brexit=0.25;
    }
  else if(high_pi==1)
    {
      pi_brexit=0.75;
    }
  else
    {
      pi_brexit=0.5;
    }
  quant_exercise();

  BREAK2;
  fprintf(logfile, KGRN "\nProgram complete!\n" RESET);
  time(&stop);
  fprintf(logfile, KGRN "Total runtime: %0.2f\n\n" RESET,difftime(stop,start));

  return 0;
}

#endif
