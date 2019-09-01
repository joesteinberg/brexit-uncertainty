#ifndef __MAIN_C__
#define __MAIN_C__

#include "globals.h"
#include "calibrate.h"
#include "eqm.h"

#define BREAK1  fprintf(logfile, KNRM "---------------------------------------------------------------------------------------\n" RESET);

#define BREAK2 fprintf(logfile, KNRM "\n\n\n////////////////////////////////////////////////////////////////////////////////////////\n" RESET);

uint parse_args(int argc, char **argv)
{
  int opt = 0;

  while((opt = getopt(argc, argv, "crwsenflhptmzab")) != -1){
    switch(opt){
    case 'c':
      recal = 1;
      break;
    case 'r':
      read_seed = 1;
      read_stoch_seed = 1;
      break;
    case 'w':
      write_seed = 1;
      write_stoch_seed = 1;
      break;
    case 'l':
      low_pi = 1;
      high_pi = 0;
      break;
    case 'h':
      low_pi = 0;
      high_pi = 1;
      break;
    case 'f':
      fin_aut = 1;
      break;
    case 'm':
      comp_mkts = 1;
      break;
    case 's':
      fixl = 2;
      break;
    case 'e':
      eqkappa = 1;
      nokappa = 0;
      break;
    case 'n':
      eqkappa = 0;
      nokappa = 1;
      break;
    case 'p':
      psi_sens=1;
      break;
    case 'z':
      zeta_sens=1;
      break;
    case 'a':
      do_china=1;
      break;
    case 'b':
      higher_ch_trd_costs=1;
      break;
    case '?':
      fprintf(logfile,
	      "\nIncorrect command line option: %c. Possible options: -e -n -f -m -h -l -s -p -z -c -r -w -a -b \n",opt);

      fprintf(logfile,"\t-e: model with kappa0 = kappa1 > 0\n");
      fprintf(logfile,"\t-n: model with kappa0 = kappa1 = 0\n");
      fprintf(logfile,"\t-f: model with financial autarky\n");
      fprintf(logfile,"\t-h: model with high probability of soft Brexit\n");
      fprintf(logfile,"\t-l: model with low probability of soft Brexit\n");
      fprintf(logfile,"\t-s: sticky wages model\n");
      fprintf(logfile,"\t-p: higher risk aversion\n");
      fprintf(logfile,"\t-z: lower Armington elasticity\n");
      fprintf(logfile,"\t-c: redo calibration (otherwise load params from file)\n");
      fprintf(logfile,"\t-r: read initial guesses for solvers from seed files\n");
      fprintf(logfile,"\t-w: write equilibrium solutions to new seed files\n");
      fprintf(logfile,"\t-a: do China-like exercise (to investigate role of asymmetries in timing)\n");
      fprintf(logfile,"\t-b: higher trade costs in China exercise\n");

      return 1;
      break;
    }
  }
  
  return 0;
}

int quant_exercise()
{
  // ---------------------------------------------------------------------------------------------------
  // set up variable and parameter structures
  uint ih, it;
  eqm * e;
  stoch_eqm * se;
  stoch_eqm_ch * se_ch;
  stoch_params * sp;
  stoch_params_ch * sp_ch;
  params * p;

  if(calibrate())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  for(it=0; it<NTH; it++)
    {
      init_vars(&(eee0[it]));
      init_vars(&(eee0_ch[it]));

      se = &(sss[it]);      
      for(ih=0; ih<NHIST; ih++)
	{
	  init_vars( &((se->eee)[ih]) );
	}

      se_ch = &(sss_ch[it]);      
      for(ih=0; ih<NHIST_CH; ih++)
	{
	  init_vars( &((se_ch->eee)[ih]) );
	}
    }

  solver_verbose=1;

  if(!do_china)
    {
      // no-Brexit counterfactual
      BREAK1;
      scenario = 0;
      set_neqm();
      fprintf(logfile,KBLU "\nSolving for counterfactual no-Brexit equilibrium...\n" RESET);
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
      BREAK1;
      scenario = 2;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp1[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nSolving for perfect-foresight soft Brexit equilibrium (referendum onwards)...\n" RESET);
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
      BREAK1;
      scenario = 3;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp2[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nSolving for perfect-foresight hard Brexit equilibrium (referendum onwards)...\n" RESET);
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

      // stochastic model
      BREAK1;
      scenario = 1;
      set_neqm();
      for(it=0; it<NTH; it++)
	{
	  set_tariffs(&(sppp[it]));
	}

      fprintf(logfile,KBLU "\nSolving for stochastic model with uncertainty about trade costs after Brexit...\n" RESET);
      if(solve_stoch_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      se = &(sss[0]);
      sp = &(sppp[0]);
      calc_stoch_welfare();
      calc_ce_welfare1();
 
      e = &((se->eee)[0]);
      p = &((sp->ppp)[0]);
      calc_welfare(e,p);

      e = &((se->eee)[1]);
      p = &((sp->ppp)[1]);
      calc_welfare(e,p);
  
      e = &((se->eee)[2]);
      p = &((sp->ppp)[2]);
      calc_welfare(e,p);

      // perfect foresight (soft Brexit, TVOTE onwards) 
      BREAK1;
      scenario = 4;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp1[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nSolving for perfect-foresight soft Brexit equilibrium (vote onwards)...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      // perfect foresight (hard Brexit, TVOTE onwards) 
      BREAK1;
      scenario = 5;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp2[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nSolving for perfect-foresight hard Brexit equilibrium (vote onwards)...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      // write stochastic vars
      calc_ce_welfare2();

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

  // ---------------------------------------------------------------------------------------------------
  if(do_china)
    {

#ifdef CH_COMPLEX
      if((int)NHIST_CH != 2+TCH2-TCH1)
	{
	  fprintf(logfile,KRED "NHIST_CH != 2+TCH2-TCH1!\n" RESET);
	  return 1;
	}

      // rest pi_brexit so that hard Brexit equally likely as in simple China exercise
      if(TCH2==TBREXIT)
	{
	  pi_brexit = 0.7937;
	}
      else
	{
	  pi_brexit = 0.87055;
	}
#endif
      
      // China-like exercise, part 1: counterfactual where trade costs are high forever
      BREAK1;
      scenario=6;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp0_ch[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nChina-like exercise, part 1: no-reform counterfactual...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      calc_welfare(&(eee0_ch[0]), &(ppp0_ch[0]));

#ifdef CH_SIMPLE
      if((int)TCH2==(int)TBREXIT)
	{
	  if(higher_ch_trd_costs)
	    {
	      write_eqm_vars(&(ppp0_ch[0]),&(eee0_ch[0]),"vars_ch_counter_uk_hi",0);
	      write_eqm_vars(&(ppp0_ch[0]),&(eee0_ch[0]),"vars_ch_counter_eu_hi",1);
	      write_eqm_vars(&(ppp0_ch[0]),&(eee0_ch[0]),"vars_ch_counter_rw_hi",2);
	      write_bgp_vars(&(eee0_ch[0]),"bgp_ch_counter_hi");
	    }
	  else
	    {
	      write_eqm_vars(&(ppp0_ch[0]),&(eee0_ch[0]),"vars_ch_counter_uk",0);
	      write_eqm_vars(&(ppp0_ch[0]),&(eee0_ch[0]),"vars_ch_counter_eu",1);
	      write_eqm_vars(&(ppp0_ch[0]),&(eee0_ch[0]),"vars_ch_counter_rw",2);
	      write_bgp_vars(&(eee0_ch[0]),"bgp_ch_counter");
	    }
	}
#endif

#ifdef CH_SIMPLE
      // China-like exercise, part 2: perfect-foresight equilibrium where trade costs fall permanently
      BREAK1;
      scenario=7;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp1_ch[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nChina-like exercise, part 2: perfect-foresight equilibrium where trade costs fall permanently...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      calc_welfare(&(eee1_ch[0]), &(ppp1_ch[0]));

      if(higher_ch_trd_costs)
	{
	  write_eqm_vars(&(ppp1_ch[0]),&(eee1_ch[0]),"vars_ch_det_uk_hi",0);
	  write_eqm_vars(&(ppp1_ch[0]),&(eee1_ch[0]),"vars_ch_det_eu_hi",1);
	  write_eqm_vars(&(ppp1_ch[0]),&(eee1_ch[0]),"vars_ch_det_rw_hi",2);
	  write_bgp_vars(&(eee1_ch[0]),"bgp_ch_det_hi");
	}
      else
	{
	  write_eqm_vars(&(ppp1_ch[0]),&(eee1_ch[0]),"vars_ch_det_uk",0);
	  write_eqm_vars(&(ppp1_ch[0]),&(eee1_ch[0]),"vars_ch_det_eu",1);
	  write_eqm_vars(&(ppp1_ch[0]),&(eee1_ch[0]),"vars_ch_det_rw",2);
	  write_bgp_vars(&(eee1_ch[0]),"bgp_ch_det");
	}

      // China-like exercise, part 3: perfect-foresight equilibrium where trade costs fall then rise
      BREAK1;
      scenario=8;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp2_ch[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nChina-like exercise, part 3: perfect-foresight equilibrium where fall but rise again in TCH2...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      calc_welfare(&(eee2_ch[0]), &(ppp2_ch[0]));
#else
      // China-like exercise, part 2.x: perfect-foresight with reversion in TCH1+x
      uint is;
      for(is=0; is<NHIST_CH-1; is++)
	{
	  BREAK1;
	  scenario=1000+is;
	  for(it=0; it<NTH; it++)
	    {
	      set_tariffs2(&(ppp1_ch[is][it]),scenario);
	    }
	  set_neqm();
	  if(is<NHIST_CH-2)
	    {
	      fprintf(logfile,KBLU "\nChina-like exercise, part 2.%d: perfect-foresight equilibrium where fall and rise again in %d years...\n" RESET,is,is+1);
	    }
	  else
	    {
	      fprintf(logfile,KBLU "\nChina-like exercise, part 2.%d: perfect-foresight equilibrium where fall permanently...\n" RESET,is);	      
	    }
	  if(solve_eqm())
	    {
	      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	      return 1;
	    }
	  calc_welfare(&(eee1_ch[is][0]), &(ppp1_ch[is][0]));
	  
	  if(is==NHIST_CH-1)
	    {
	      if(higher_ch_trd_costs)
		{
		  write_eqm_vars(&(ppp1_ch[is][0]),&(eee1_ch[is][0]),"vars_ch_det_uk_hi",0);
		  write_eqm_vars(&(ppp1_ch[is][0]),&(eee1_ch[is][0]),"vars_ch_det_eu_hi",1);
		  write_eqm_vars(&(ppp1_ch[is][0]),&(eee1_ch[is][0]),"vars_ch_det_rw_hi",2);
		  write_bgp_vars(&(eee1_ch[is][0]),"bgp_ch_det_hi");
		}
	      else
		{
		  write_eqm_vars(&(ppp1_ch[is][0]),&(eee1_ch[is][0]),"vars_ch_det_uk",0);
		  write_eqm_vars(&(ppp1_ch[is][0]),&(eee1_ch[is][0]),"vars_ch_det_eu",1);
		  write_eqm_vars(&(ppp1_ch[is][0]),&(eee1_ch[is][0]),"vars_ch_det_rw",2);
		  write_bgp_vars(&(eee1_ch[is][0]),"bgp_ch_det");
		}
	    }
	}
#endif

      // stochastic model for China-like exercise
      BREAK1;
      scenario = 9;
      set_neqm();
      for(it=0; it<NTH; it++)
	{
	  set_tariffs_ch(&(sppp_ch[it]));
	}

#ifdef CH_SIMPLE
      fprintf(logfile,KBLU "\nChina-like exercise, part 4: stochastic model\n" RESET);
#else
      fprintf(logfile,KBLU "\nChina-like exercise, part 3: stochastic model\n" RESET);
#endif

      //eval_eqm_once_flag=1;
      if(solve_stoch_eqm_ch())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      se_ch = &(sss_ch[0]);
      sp_ch = &(sppp_ch[0]);
      calc_stoch_ch_welfare();

      uint ih;
      for(ih=0; ih<NHIST_CH; ih++)
	{
	  e = &((se_ch->eee)[ih]);
	  p = &((sp_ch->ppp)[ih]);
	  calc_welfare(e,p);
	}
      calc_ce_ch_welfare1();
      calc_ce_ch_welfare2();

#ifdef CH_SIMPLE
      // China-like exercise, part 4: perfect-foresight equilibrium with no reversion from TCH1
      BREAK1;
      scenario=10;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp1_ch[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nChina-like exercise, part 5: 2nd perfect-foresight equilibrium where trade costs fall permanently...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      // China-like exercise, part 5: perfect-foresight equilibrium with reversion from TCH1
      BREAK1;
      scenario=11;
      for(it=0; it<NTH; it++)
	{
	  set_tariffs2(&(ppp2_ch[it]),scenario);
	}
      set_neqm();
      fprintf(logfile,KBLU "\nChina-like exercise, part 6: 2rd perfect-foresight equilibrium where fall in TCH0 and rise again in TCH1...\n" RESET);
      if(solve_eqm())
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
#else
      // China-like exercise, part 4.x: perfect-foresight with reversion in TCH1+x
      for(is=0; is<NHIST_CH-1; is++)
	{
	  BREAK1;
	  scenario=2000+is;
	  for(it=0; it<NTH; it++)
	    {
	      set_tariffs2(&(ppp1_ch[is][it]),scenario);
	    }
	  set_neqm();
	  if(is<NHIST_CH-2)
	    {
	      fprintf(logfile,KBLU "\nChina-like exercise, part 4.%d: perfect-foresight equilibrium where fall in TCH0 and rise again in TCH1+%d...\n" RESET,is,is+1);
	    }
	  else
	    {
	      fprintf(logfile,KBLU "\nChina-like exercise, part 4.%d: perfect-foresight equilibrium where fall in TCH0 and never rise...\n" RESET,is);	      
	    }
	  if(solve_eqm())
	    {
	      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	      return 1;
	    }
	  calc_welfare(&(eee1_ch[is][0]), &(ppp1_ch[is][0]));
	}      
#endif

      // write stochastic vars
      calc_ce_ch_welfare2();
      e = &((se_ch->eee)[NHIST_CH-1]);
      p = &((sp_ch->ppp)[NHIST_CH-1]);

#ifdef CH_SIMPLE
      if(higher_ch_trd_costs)
	{
	  write_eqm_vars(p,e,"vars_stoch_ch_uk_hi",0);
	  write_eqm_vars(p,e,"vars_stoch_ch_eu_hi",1);
	  write_eqm_vars(p,e,"vars_stoch_ch_rw_hi",2);
	}
      else
	{
	  write_eqm_vars(p,e,"vars_stoch_ch_uk",0);
	  write_eqm_vars(p,e,"vars_stoch_ch_eu",1);
	  write_eqm_vars(p,e,"vars_stoch_ch_rw",2);
	}
#else
      if(higher_ch_trd_costs)
	{
	  write_eqm_vars(p,e,"vars_stoch_ch_complex_uk_hi",0);
	  write_eqm_vars(p,e,"vars_stoch_ch_complex_eu_hi",1);
	  write_eqm_vars(p,e,"vars_stoch_ch_complex_rw_hi",2);
	}
      else
	{
	  write_eqm_vars(p,e,"vars_stoch_ch_complex_uk",0);
	  write_eqm_vars(p,e,"vars_stoch_ch_complex_eu",1);
	  write_eqm_vars(p,e,"vars_stoch_ch_complex_rw",2);
	}
#endif
    }

  return 0;

}

int main(int argc, char * argv[])
{
  do_china=0;
  higher_ch_trd_costs=0;
  par=0;
  psi_sens=0;
  zeta_sens=0;
  comp_mkts=0;
  leontief = 1;
  low_pi = 0;
  high_pi = 0;
  fin_aut = 0;
  nokappa = 0;
  eqkappa = 0;
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

  BREAK2;
  fprintf(logfile, KGRN "\nBrexit and the Macroeconomic Impact of Trade Policy Uncertainty" RESET);
  fprintf(logfile, KGRN "\nJoseph Steinberg, University of Toronto" RESET);
  fprintf(logfile, KGRN "\nLast updated: May 2017" RESET);
  fprintf(logfile, KNRM "\n" RESET);

  BREAK1;
  fprintf(logfile, KBLU "\nSetting up environment...\n\n" RESET);

  if(parse_args(argc,argv))
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    };

  if(recal==1)
    {
      fprintf(logfile, KBLU "\tRecalibrating all parameters\n" RESET);      
    }
  else
    {
      fprintf(logfile, KBLU "\tLoading parameters from file\n" RESET);      
    }

  if(read_seed == 1 && read_stoch_seed == 1)
    {
      fprintf(logfile, KBLU "\tReading seed files\n" RESET);      
    }
  if(write_seed == 1 && write_stoch_seed == 1)
    {
      fprintf(logfile, KBLU "\tWriting new seed files\n" RESET);      
    }
  if(fixl == 1)
    {
      fprintf(logfile, KBLU "\tFixed labor supply\n" RESET);
    }
  if(fixl == 2)
    {
      fprintf(logfile, KBLU "\tSticky wages\n" RESET);
    }
  if(fixl == 0)
    {
      fprintf(logfile, KBLU "\tEndogenous labor supply\n" RESET);
    }
  if(eqkappa == 1)
    {
      fprintf(logfile, KBLU "\tModel with kappa0 = kappa1 > 0 (static export entry decision)\n" RESET);
    }
  else if(nokappa==1)
    {
      fprintf(logfile, KBLU "\tModel with kappa0 = kappa1 = 0 (no extensive margin in trade)\n" RESET);
    }
  else
    {
      fprintf(logfile, KBLU "\tModel with kappa0 > kappa1 > 0 (baseline with forward-looking export entry)\n" RESET);
    }

  if(fin_aut)
    {
      fprintf(logfile, KBLU "\tModel with financial autarky\n" RESET);
    }
  if(comp_mkts)
    {
      fprintf(logfile, KBLU "\tModel with complete markets\n" RESET);
    }

  if(psi_sens)
    {
      fprintf(logfile, KBLU "\tModel with higher risk aversion\n" RESET);
    }

  if(zeta_sens)
    {
      fprintf(logfile, KBLU "\tModel with lower trade elasticity\n" RESET);
    }

  if(low_pi == 1)
    {
      fprintf(logfile, KBLU "\tProbability of hard Brexit = 75pct\n" RESET);
    }
  else if(high_pi == 1)
    {
      fprintf(logfile, KBLU "\tProbability of hard Brexit = 25pct\n" RESET);
    }
  else
    {
      fprintf(logfile, KBLU "\tProbability of hard Brexit = 50pct\n" RESET);
    }
  
  if(do_china==1)
    {
#ifdef CH_SIMPLE
      fprintf(logfile, KBLU "\tDoing China-like exercise after Brexit exercises (simple version)\n" RESET);
#else
      fprintf(logfile, KBLU "\tDoing China-like exercise after Brexit exercises (complex version)\n" RESET);
#endif
      if(higher_ch_trd_costs==1)
	{
	  fprintf(logfile, KBLU "\tHigher trade costs in China-like exercise\n" RESET);
	}
    }

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

  BREAK2;
  fprintf(logfile,KGRN "Solving for equilibria...\n\n\n" RESET);
  pi_vote=0.75;
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

  BREAK1;
  fprintf(logfile, KGRN "\nProgram complete!\n\n" RESET);

  return 0;
}

#endif
