#ifndef __SOLVER_H__
#define __SOLVER_H__

#ifdef _STATIC
#include "globals_static.h"
#else
#include "globals.h"
#endif

//#include "calibrate.h"
//#include "eqm.h"

#include "gnewton.h"

uint par;
uint fcnt;
uint jcnt;
uint solver_mem_alloced;
uint solver_n;
uint solver_verbose;
gsl_vector * solver_x;
gsl_vector * xh[NTH];
gsl_vector * f0[NTH];
gsl_vector * fh[NTH];
extern const double epsjac;
extern const double root_tol;
extern const uint max_root_iter;

int jacobian(
	     int (* F) (const gsl_vector * x, void * data, gsl_vector * f),
	     const gsl_vector * x,
	     gsl_matrix * J,
	     uint set_f0
	     );

uint find_root_deriv_mkl(gsl_multiroot_function_fdf * f);

static inline void free_solver_mem()
{
  if(solver_mem_alloced)
    {
      uint i;
      gsl_vector_free(solver_x);
      for(i=0; i<NTH; i++)
	{
	  gsl_vector_free(xh[i]);
	  gsl_vector_free(f0[i]);
	  gsl_vector_free(fh[i]);
	}
    }
  solver_mem_alloced = 0;
}

static inline void alloc_solver_mem()
{
  uint i;
  if(solver_mem_alloced)
    {
      free_solver_mem();
    }
  solver_x = gsl_vector_calloc(solver_n);
  for(i=0; i<NTH; i++)
    {
      xh[i] = gsl_vector_calloc(solver_n);
      f0[i] = gsl_vector_calloc(solver_n);
      fh[i] = gsl_vector_calloc(solver_n);
    }
  solver_mem_alloced = 1;
}

#endif
