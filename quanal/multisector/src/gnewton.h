#ifndef __GNEWTON_H__
#define __GNEWTON_H__

#ifdef _STATIC
#include "globals_static.h"
#else
#include "globals.h"
#endif

typedef struct
{
  gsl_multiroot_function_fdf * FDF;
  int n;
  double phi;
  gsl_vector * x;
  gsl_vector * x_trial;
  gsl_vector * f;
  gsl_vector * dx;
  gsl_matrix * J;
  gsl_matrix * LU;
}
gnewton_state_t;

gnewton_state_t * gnewton_alloc (gsl_multiroot_function_fdf * FDF);

int gnewton_set (gnewton_state_t * state, const gsl_vector * x);

int gnewton_iterate (gnewton_state_t * state);

void gnewton_free (gnewton_state_t * state);

#endif
