#ifndef __GNEWTON_C__
#define __GNEWTON_C__

#include "gnewton.h"

gnewton_state_t * gnewton_alloc (gsl_multiroot_function_fdf * FDF)
{
  gnewton_state_t * state = (gnewton_state_t * )malloc(sizeof(gnewton_state_t));

  int n = FDF -> n;
  state -> FDF = FDF;
  state -> n = n;
  
  gsl_vector * x;
  gsl_vector * x_trial;
  gsl_vector * f;
  gsl_matrix * J;
  gsl_matrix * LU;
  gsl_vector * dx;

  J = gsl_matrix_calloc (n,n);
  if (J == 0) 
    {
      GSL_ERROR ("failed to allocate space for J", GSL_ENOMEM);
    }
  state->J = J;

  LU = gsl_matrix_calloc (n,n);
  if (LU == 0) 
    {
      gsl_matrix_free(J);
      GSL_ERROR ("failed to allocate space for J", GSL_ENOMEM);
    }
  state->LU = LU;

  x = gsl_vector_calloc (n);
  if (x == 0)
    {
      gsl_matrix_free(LU);
      gsl_matrix_free(J);

      GSL_ERROR ("failed to allocate space for x", GSL_ENOMEM);
    }
  state->x = x;

  x_trial = gsl_vector_calloc (n);
  if (x_trial == 0)
    {
      gsl_matrix_free(LU);
      gsl_vector_free(x);
      gsl_matrix_free(J);

      GSL_ERROR ("failed to allocate space for x", GSL_ENOMEM);
    }
  state->x_trial = x_trial;

  f = gsl_vector_calloc (n);
  if (f == 0)
    {
      gsl_matrix_free(LU);
      gsl_vector_free(x);
      gsl_matrix_free(J);
      gsl_vector_free(x_trial);

      GSL_ERROR ("failed to allocate space for f", GSL_ENOMEM);
    }
  state->f = f;

  dx = gsl_vector_calloc (n);
  if (dx == 0)
    {
      gsl_vector_free(x);
      gsl_vector_free(x_trial);
      gsl_matrix_free(J);
      gsl_matrix_free(LU);
      gsl_vector_free(f);

      GSL_ERROR ("failed to allocate space for dx", GSL_ENOMEM);
    }
  state->dx = dx;

  return state;
}

int gnewton_set (gnewton_state_t * state, const gsl_vector * x)
{
  size_t i, n = state->n ;

  gsl_vector_memcpy(state->x,x);

  GSL_MULTIROOT_FN_EVAL_F_DF (state->FDF, state->x, state->f, state->J);

  for (i = 0; i < n; i++)
    {
      gsl_vector_set (state->dx, i, 0.0);
    }

  state->phi = enorm(state->f);

  return GSL_SUCCESS;
}

int gnewton_iterate (gnewton_state_t * state)
{
  double t = 0.0, phi0 = 0.0, phi1 = 0.0;

  size_t i = 0;

  MKL_INT n = state->n;
  MKL_INT nrhs = 1;
  MKL_INT lda = n;
  MKL_INT ipiv[n];
  MKL_INT ldb = 1;
  MKL_INT info = 0;

  gsl_matrix_memcpy (state->LU, state->J);
  gsl_vector_memcpy (state->dx, state->f);

  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, (state->LU)->data, lda, ipiv,
                        (state->dx)->data, ldb);

  if( info > 0 )
    {
      fprintf(logfile,KRED "The diagonal element of the triangular factor of A,\n" RESET);
      fprintf(logfile,KRED "U(%i,%i) is zero, so that A is singular;\n" RESET, info, info );
      fprintf(logfile,KRED "the solution could not be computed.\n" RESET);
      return GSL_EDOM;
    }

  t = 1;
  if(small_gnewton_step_flag)
    {
      t=0.05;
    }

  phi0 = state->phi;

 new_step:

  for (i = 0; i < n; i++)
    {
      double di = gsl_vector_get (state->dx, i);
      double xi = gsl_vector_get (state->x, i);
      gsl_vector_set (state->x_trial, i, xi - t*di);
      if(gsl_isnan(xi-t*di))
	{
	  return GSL_EBADFUNC;
	}
    }
  
  { 
    int status = GSL_MULTIROOT_FN_EVAL_F (state->FDF, state->x_trial, state->f);
    
    if (status != GSL_SUCCESS)
      {
        return GSL_EBADFUNC;
      }
  }
  
  phi1 = enorm (state->f);

  if(gsl_isnan(phi1) || gsl_isinf(phi1))
    {
      t *= 0.1;
      goto new_step;
    }
  if (phi1 > phi0 && t > GSL_DBL_EPSILON)  
    {
      /* full step goes uphill, take a reduced step instead */

      double theta = phi1 / phi0;
      double u = (sqrt(1.0 + 6.0 * theta) - 1.0) / (3.0 * theta);

      t *= u ;
     
      goto new_step;
    }

  /* copy x_trial into x */

  gsl_vector_memcpy (state->x, state->x_trial);

  for (i = 0; i < n; i++)
    {
      double di = gsl_vector_get (state->dx, i);
      gsl_vector_set (state->dx, i, -t*di);
    }

  { 
    int status = GSL_MULTIROOT_FN_EVAL_DF (state->FDF, state->x, state->J);
    
    if (status != GSL_SUCCESS)
      {
        return GSL_EBADFUNC;
      }
  }

  state->phi = phi1;

  mkl_free_buffers();

  return GSL_SUCCESS;
}

void gnewton_free (gnewton_state_t * state)
{
  gsl_vector_free(state->dx);
  gsl_vector_free(state->x);
  gsl_vector_free(state->f);
  gsl_vector_free(state->x_trial);
  gsl_matrix_free(state->J);
  gsl_matrix_free(state->LU);
  free(state);
}


#endif
