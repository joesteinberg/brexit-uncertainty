#ifndef __MAIN_C__
#define __MAIN_C__

#include "globals.h"
#include "calibrate.h"
#include "eqm.h"

uint parse_args(int argc, char **argv)
{
  int opt = 0;

  while((opt = getopt(argc, argv, "rwfmzlsxp")) != -1){
    switch(opt){
    case 'r':
      read_seed = 1;
      read_stoch_seed = 1;
      break;
    case 'w':
      write_seed = 1;
      write_stoch_seed = 1;
      break;
    case 'f':
      f_adj_cost = 1;
      m_adj_cost = 1;
      break;
    case 'z':
      was_flag = 1;
      break;
    case 'l':
      fixl = 1;
      break;
    case 'x':
      fixl = 2;
      break;
    case 'p':
      tfp_flag = 1;
      break;
    case 's':
      do_sensitivity = 1;
      break;
    case '?':
      fprintf(logfile,"\nIncorrect command line option: %c. Possible options: -r -w -f -z -l -x -s\n",opt);
      fprintf(logfile,"\t-r: read initial guesses for solvers from seed files\n");
      fprintf(logfile,"\t-w: write equilibrium solutions to new seed files\n");
      fprintf(logfile,"\t-f: model with import adjustment frictions\n");
      fprintf(logfile,"\t-z: model with import adjustment frictions, adj. cost. params set to zero\n");
      fprintf(logfile,"\t-l: fixed labor supply\n");
      fprintf(logfile,"\t-x: sticky wages model\n");
      fprintf(logfile,"\t-p: model with productivity losses\n");
      fprintf(logfile,"\t-s: conduct additional sensitivity analyses after solving baseline model\n");
      return 1;
      break;
    }
  }
  
  return 0;
}

int quant_exercise()
{
  // -----------------------------------------------------------------------------------------------------------
  // set up variable and parameter structures
  uint ih, it;
  stoch_eqm * se;
  stoch_params * sp;
  for(it=0; it<NTH; it++)
    {
      init_vars(&(eee0[it]));
      if(calibrate(&(ppp0[it])))
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      if(calibrate(&(ppp1[it])))
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}
      if(calibrate(&(ppp2[it])))
	{
	  fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	  return 1;
	}

      se = &(sss[it]);
      sp = &(sppp[it]);
      
      for(ih=0; ih<NHIST; ih++)
	{
	  init_vars( &((se->eee)[ih]) );
	  if(calibrate( &((sp->ppp)[ih]) ))
	    {
	      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
	      return 1;
	    }
	}
    }
  if(sensitivity == 0 && write_params(&(ppp0[0])))
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  // -------------------------------------------------------------------------------------------------------
  // no-Brexit counterfactual

  fprintf(logfile, KNRM "----------------------------------------------------------------------\n" RESET);

  scenario = 0;
  set_neqm();

  fprintf(logfile,KBLU "\nSolving for counterfactual no-Brexit equilibrium...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  calc_welfare(&(eee0[0]), &(ppp0[0]));

  write_eqm_vars(&(eee0[0]),"vars_nobrexit_uk",0);
  write_eqm_vars(&(eee0[0]),"vars_nobrexit_eu",1);
  write_eqm_vars(&(eee0[0]),"vars_nobrexit_rw",2);

  // ---------------------------------------------------------------------------------------------------------
  // stochastic model

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
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
  calc_stoch_welfare();
 
  eqm * e = &((se->eee)[0]);
  calc_welfare(e, &(ppp0[0]));
  //write_eqm_vars(e,"vars_stoch_stay_uk",0);
  //write_eqm_vars(e,"vars_stoch_stay_eu",1);
  //write_eqm_vars(e,"vars_stoch_stay_rw",2);
 
  e = &((se->eee)[1]);
  calc_welfare(e, &(ppp0[1]));
  //write_eqm_vars(e,"vars_stoch_opt_uk",0);
  //write_eqm_vars(e,"vars_stoch_opt_eu",1);
  //write_eqm_vars(e,"vars_stoch_opt_rw",2);
  
  e = &((se->eee)[2]);
  calc_welfare(e, &(ppp0[2]));
  //write_eqm_vars(e,"vars_stoch_pes_uk",0);
  //write_eqm_vars(e,"vars_stoch_pes_eu",1);
  //write_eqm_vars(e,"vars_stoch_pes_rw",2);

  // -------------------------------------------------------------------------------
  // perfect foresight model with soft Brexit (TREF onwards)

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  scenario = 2;
  set_neqm();
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp1[it]),scenario);
    }

  fprintf(logfile,KBLU "\nSolving for perfect foresight model with soft Brexit (TREF onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  calc_welfare(&(eee1[0]), &(ppp1[0]));
  write_eqm_vars(&(eee1[0]),"vars_det_opt_uk",0);
  write_eqm_vars(&(eee1[0]),"vars_det_opt_eu",1);
  write_eqm_vars(&(eee1[0]),"vars_det_opt_rw",2);

  // -------------------------------------------------------------------------------
  //  perfect foresight model with hard Brexit (TREF onwards)

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  scenario = 3;
  set_neqm();
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp2[it]),scenario);
    }

  fprintf(logfile,KBLU "\nSolving for perfect foresight model with hard Brexit (TREF onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  calc_welfare(&(eee2[0]), &(ppp2[0]));
  write_eqm_vars(&(eee2[0]),"vars_det_pes_uk",0);
  write_eqm_vars(&(eee2[0]),"vars_det_pes_eu",1);
  write_eqm_vars(&(eee2[0]),"vars_det_pes_rw",2);

  calc_ce_welfare1();

  // -------------------------------------------------------------------------------
  // perfect foresight model with soft Brexit (TVOTE onwards)

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  scenario = 4;
  set_neqm();
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp1[it]),scenario);
    }

  fprintf(logfile,KBLU "\nSolving for perfect foresight model with soft Brexit (TVOTE onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  calc_welfare(&(eee1[0]), &(ppp1[0]));
  write_eqm_vars(&(eee1[0]),"vars_det_opt2_uk",0);
  write_eqm_vars(&(eee1[0]),"vars_det_opt2_eu",1);
  write_eqm_vars(&(eee1[0]),"vars_det_opt2_rw",2);

  // -------------------------------------------------------------------------------
  // perfect foresight model with hard Brexit (TVOTE onwards)

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  scenario = 5;
  set_neqm();
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp2[it]),scenario);
    }

  fprintf(logfile,KBLU "\nSolving for perfect foresight model with hard Brexit (TVOTE onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  calc_welfare(&(eee2[0]), &(ppp2[0]));
  write_eqm_vars(&(eee2[0]),"vars_det_pes2_uk",0);
  write_eqm_vars(&(eee2[0]),"vars_det_pes2_eu",1);
  write_eqm_vars(&(eee2[0]),"vars_det_pes2_rw",2);

  calc_ce_welfare2();

  // -------------------------------------------------------------------------------
  // write stochastic vars

  e = &((se->eee)[0]); 
  write_eqm_vars(e,"vars_stoch_stay_uk",0);
  write_eqm_vars(e,"vars_stoch_stay_eu",1);
  write_eqm_vars(e,"vars_stoch_stay_rw",2);
 
  e = &((se->eee)[1]);
  write_eqm_vars(e,"vars_stoch_opt_uk",0);
  write_eqm_vars(e,"vars_stoch_opt_eu",1);
  write_eqm_vars(e,"vars_stoch_opt_rw",2);
  
  e = &((se->eee)[2]);
  write_eqm_vars(e,"vars_stoch_pes_uk",0);
  write_eqm_vars(e,"vars_stoch_pes_eu",1);
  write_eqm_vars(e,"vars_stoch_pes_rw",2);

  return 0;

}

int main(int argc, char * argv[])
{
  iceberg=1;
  f_adj_cost=0;
  m_adj_cost=0;
  eval_eqm_once_flag=0;
  read_seed=0;
  write_seed=0;
  read_stoch_seed=0;
  write_stoch_seed=0;
  fixl=0;
  do_sensitivity=0;
  tfp_flag=0;
  sensitivity=0;
  logfile = stdout;

  fprintf(logfile, KGRN "\nBrexit and the Macroeconomic Impact of Trade Policy Uncertainty" RESET);
  fprintf(logfile, KGRN "\nJoseph Steinberg, University of Toronto" RESET);
  fprintf(logfile, KNRM "\n" RESET);

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KBLU "\nSetting up environment...\n\n" RESET);

  if(parse_args(argc,argv))
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    };

  if(read_seed == 1 && read_stoch_seed == 1)
    {
      fprintf(logfile, KBLU "\tReading seed files\n" RESET);      
    }
  if(write_seed == 1 && write_stoch_seed == 1)
    {
      fprintf(logfile, KBLU "\tWriting new seed files\n" RESET);      
    }
  if(f_adj_cost == 0 && m_adj_cost == 0)
    {
      fprintf(logfile, KBLU "\tNo import adjustment costs\n" RESET);      
    }
  if(m_adj_cost == 1)
    {
      fprintf(logfile, KBLU "\tIntermediate import adjustment costs\n" RESET);      
    }
  if(f_adj_cost == 1)
    {
      fprintf(logfile, KBLU "\tFinal import adjustment costs\n" RESET);
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
  if(tfp_flag == 1)
    {
      fprintf(logfile, KBLU "\tTFP losses in addition to trade costs increases\n" RESET);
    }
  if(do_sensitivity == 1)
    {
      fprintf(logfile, KBLU "\tConducting additional sensitivity analyses\n" RESET);
    }
  if(was_flag == 1)
    {
      fprintf(logfile, KBLU "\tAdj. cost equations but zero cost params. for seed file construction\n" RESET);
    }

  // ---------------------------------------------------------------------------------
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

  fprintf(logfile, KNRM "\n\n\n////////////////////////////////////////////////////////////////////////\n" RESET);
  fprintf(logfile,KGRN "Baseline quantitative exercise\n\n\n" RESET);
  pi_vote=0.75;
  pi_brexit=0.5;
  quant_exercise();

  if(do_sensitivity)
    {
      // 25% chance of good state
      fprintf(logfile, KNRM "\n\n\n////////////////////////////////////////////////////////////////////////\n" RESET);
      fprintf(logfile,KGRN "Sensitivity analysis 1: pi = 0.25\n\n\n" RESET);
      sensitivity=1;
      pi_brexit=0.25;
      quant_exercise();

      // 75% chance of good state
      fprintf(logfile, KNRM "\n\n\n////////////////////////////////////////////////////////////////////////\n" RESET);
      fprintf(logfile,KGRN "Sensitivity analysis 2: pi = 0.75\n\n\n" RESET);
      sensitivity=2;
      pi_brexit = 0.75;
      quant_exercise();
      
      pi_brexit=0.5;
      
      // risk aversion  = 5
      fprintf(logfile, KNRM "\n\n\n////////////////////////////////////////////////////////////////////////\n" RESET);
      fprintf(logfile,KGRN "Sensitivity analysis 3: risk aversion  = 5\n\n\n" RESET);
      sensitivity=3;
      quant_exercise();

      // risk aversion  = 10
      fprintf(logfile, KNRM "\n\n\n////////////////////////////////////////////////////////////////////////\n" RESET);
      fprintf(logfile,KGRN "Sensitivity analysis 4: Cobb-Douglas Armington aggregators\n\n\n" RESET);
      sensitivity=4;
      quant_exercise();
    }


  // ---------------------------------------------------------------------------------
  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KGRN "\nProgram complete!\n\n" RESET);

  return 0;
}

#endif
