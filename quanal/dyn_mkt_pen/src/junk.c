1,  1 sec, max(f) = +2.2e-01, min(f) = -6.9e+00, max(dx) = +2.3e+00, min(dx) = -8.2e-03
2,  2 sec, max(f) = +1.8e-03, min(f) = -6.1e-02, max(dx) = +1.1e-01, min(dx) = -2.3e-04
3,  1 sec, max(f) = +9.6e-08, min(f) = -3.1e-06, max(dx) = +1.4e-03, min(dx) = -3.2e-06
4,  2 sec, max(f) = +1.4e-11, min(f) = -3.8e-11, max(dx) = +4.9e-08, min(dx) = -1.8e-10


  // perfect foresight (soft Brexit, TVOTE onwards) 
  /*BREAK2;
  scenario = 4;
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp1[it]),scenario);
    }
  set_neqm();
  fprintf(logfile,KGRN "\nSolving for perfect-foresight soft Brexit equilibrium (vote onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }

  // perfect foresight (hard Brexit, TVOTE onwards) 
  BREAK2;
  scenario = 5;
  for(it=0; it<NTH; it++)
    {
      set_tariffs2(&(ppp2[it]),scenario);
    }
  set_neqm();
  fprintf(logfile,KGRN "\nSolving for perfect-foresight hard Brexit equilibrium (vote onwards)...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
      }*/

  // write stochastic vars
  //calc_ce_welfare2();




// death probability
double death_prob(double phi, double xi, double omega)
{
  double tmp = 1.0-exp(-xi*phi);
  tmp -= omega;
  tmp = fmax(tmp,0.0);
  tmp = fmin(1.0,tmp);
  //return 1.0-tmp;
  return omega;
}

// survival probability grid associated with a given productivity grid
/*void set_surv_probs(exporter_params * ep)
{
  uint ip=0;
  for(ip=0; ip<NPHI; ip++)
    {
      ep->xi_grid[ip] = 1.0-death_prob(ep->phi_grid[ip],ep->xi,ep->omega);
    }
  return;
  }*/


// analytical entry cutoff
void set_entry_cutoff(exporter_params * ep, exporter_vars * ev)
{
  //double tmp=0.0;
  /*tmp = ev->W*ev->La*( (1.0/ep->psi) - 
				 ev->Q*ep->xi*(ep->omega-ep->psi)/ep->psi + 
				 ev->Q*ep->xi*(ep->omega-ep->psi-(1.0-ep->delta))*ep->psi );
				 tmp = tmp / ev->B;*/

  ev->cutoff = pow(tmp,1.0/(ep->theta-1.0));
}

// marketing cost
double mp_f(const exporter_params * ep, double np, double n, double L)
{
  //double top = ep->psi + n*(ep->omega - ep->psi) - np;
  //double bottom = ep->psi + n*(ep->omega - ep->psi - (1.0 - ep->delta));

  //double retval=(pow(L,ep->alpha)/(1.0-ep->beta)) * (1.0-pow(frac,1.0-ep->beta));

  return retval;
}

// derivative of marketing cost with respect to n_{t+1}
double mp_f1(const exporter_params * ep, double np, double n, double La)
{
  //double tmp1 = (ep->psi + n*(ep->omega - ep->psi - (1.0 - ep->delta)));
  //double tmp2 = (ep->psi + n*(ep->omega - ep->psi)-np);

  //double retval=(La * pow(tmp1,ep->beta-1.0) * pow(tmp2,-ep->beta);
  
  return retval;
}

// derivative of marketing cost with respect to n_t
double mp_f2(const exporter_params * ep, double np, double n, double La)
{
  //double tmp1 = (ep->psi + n*(ep->omega - ep->psi - (1.0-ep->delta)));
  //double tmp2 = (ep->psi + n*(ep->omega - ep->psi)-np);
  //double retval= -La * (ep->omega - ep->psi) * pow(tmp2,-ep->beta) * pow(tmp1,ep->beta-1.0) 
  //  + La * (ep->omega - ep->psi-(1.0-ep->delta)) * pow(tmp1,ep->beta) * pow(tmp2,1.0-ep->beta);

  return retval;
}


// Euler equation
/*double foc(double np, void * params)
{
  foc_params * params2 = (foc_params *)params;

  exporter_params * ep = params2->ep;
  exporter_vars * ev = params2->ev;
  double W = ev->W;
  double B = ev->B;
  double Q = ev->Q;
  double La = ev->La;
  int ip = params2->ip;
  int in = params2->in;

  double phi_pow = ep->phi_pow_grid[ip];
  double n = ep->n_grid[in];

  double npp = gsl_spline_eval(ev->spline[ip],np,ev->acc[ip]);
  double lhs = W*mp_f1(ep,np,n,La);
  double rhs = B*phi_pow;

  rhs = rhs - Q*ep->xi*W*mp_f2(ep,npp,np,La);

  return lhs-rhs;
  }

// root finder
int find_root_1d(gsl_function * f, double xlo, double xhi, double * x)
{
  int status = 0;
  int iter = 0;
  const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;
  gsl_root_fsolver * s = gsl_root_fsolver_alloc(T);
  
  status = gsl_root_fsolver_set(s,f,xlo,xhi);
  if(status)
    {
      printf("Error initializing root-finder!\n");
    }
  else
    {
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate(s);
	  if(status)
	    {
	      printf("Error iterating root-finder!\n");
	      break;
	    }
	  *x = gsl_root_fsolver_root(s);
	  xlo = gsl_root_fsolver_x_lower(s);
	  xhi = gsl_root_fsolver_x_upper(s);
	  status = gsl_root_test_interval(xlo,xhi,root_tol_abs,root_tol_rel);
	}while(status==GSL_CONTINUE && iter<max_root_iter);
    }

  gsl_root_fsolver_free(s);

  return status;
  
}*/

// policy function iteration driver
/*
int iterate_policy_fn(exporter_params * ep, exporter_vars * ev, double * supnorm)
{
  int ip, inp;
  *supnorm = 0.0;

  for(ip=0; ip<NPHI; ip++)
    {
      //gsl_spline_init(ev->spline[ip],ep->n_grid,ev->hn[ip],NN);
      //gsl_interp_accel_reset(ev->acc[ip]);

      for(inp=0; inp<NN; inp++)
	{
	  foc_params params = {phi_pow_grid[ip],n_grid[in],spline[ip],vf_spline[ip],acc[ip]};
	  foc_params params = {ep,ev,ip,in};
	  double lower_bound=(1.0-ep->delta)*ep->n_grid[in]+1.0e-8;
	  double upper_bound = fmin(ep->psi*(1.0-ep->n_grid[in])+ep->n_grid[in]-1.0e-8,
	  ep->n_grid[NN-1]-1.0e-8);
	  
	  // first, check if we are at a corner; namely, if
	  // marginal cost of marketing when n'=(1-delta)*n is greater than
	  // marginal benefit
	  if(foc(lower_bound,&params)>0.0)
	   {
	     ev->hn[ip][in]=(1.0-ep->delta)*ep->n_grid[in];
	    }
	  
	  //else if(foc(upper_bound,&params)<0.0)
	    {
	      ev->hn[ip][in] = upper_bound;
	    }	  
	  
	    // if not, find interior solution
	  else
	   {
	     if(in>0)
		{
		  lower_bound = fmax( (1.0-ep->delta)*ep->n_grid[in]+1.0e-8,
				      ev->hn[ip][in-1]-1.0e-4 );
		}
	      if(ip>0)
		{
		  lower_bound = fmax(lower_bound,ev->hn[ip-1][in]-1.0e-4);
		}

	      gsl_function f;
	      f.function = &foc;	      
	      f.params=&params;
	      if(find_root_1d(&f,lower_bound,upper_bound,&(ev->hn[ip][in])))
		{
		  fprintf(logfile, KRED"\nError solving first-order condition! (ip,in) = (%d,%d)\n" RESET,ip,in);
		  return 1;
		}
		}

	  // check distance from last iteration's policy function
	  double hn_old = ev->spline[ip]->y[in];
	  if(fabs(hn_old - ev->hn[ip][in])> *supnorm)
	    {
	      *supnorm = fabs(hn_old - ev->hn[ip][in]);
	      }

	}
      
    }

  return 0;
  }*/


  // pareto
  /*double phi_lo=1.0;
  double phi_hi=2.0*ep->kappa;
  linspace(phi_lo,phi_hi,NPHI,ep->phi_grid);
  double sum = 0.0;
  for(i=1; i<NPHI; i++)
    {
      ep->phi_probs[i] = pareto_cdf(ep->phi_grid[i],ep->kappa)-pareto_cdf(ep->phi_grid[i-1],ep->kappa);
      sum += ep->phi_probs[i];
    }
    ep->phi_probs[0] = 1.0 - sum;*/

  // adda and cooper method for loglinear
  /*int i;
  double m[NPHI-1];
  for(i=0; i<(NPHI-1); i++)
    {
      m[i] = gsl_cdf_ugaussian_Pinv( (i+1.0)/NPHI ) * ep->kappa;
    }

  ep->phi_grid[0] = -ep->kappa * NPHI * gsl_ran_ugaussian_pdf(m[0]/ep->kappa);
  for(i=1; i<(NPHI-1.0); i++)
    {
      ep->phi_grid[i] = -ep->kappa * NPHI * (gsl_ran_ugaussian_pdf(m[i]/ep->kappa) - 
					     gsl_ran_ugaussian_pdf(m[i-1]/ep->kappa));
    }
    ep->phi_grid[NPHI-1] = ep->kappa * NPHI * gsl_ran_ugaussian_pdf(m[NPHI-2]/ep->kappa);*/
