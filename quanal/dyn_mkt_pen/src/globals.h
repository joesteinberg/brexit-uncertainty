#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <mkl.h>
#include <mkl_lapacke.h>

#define NC 3 // number of countries
#define NT 50 // number of time periods
#define NPHI 300 //150 // number of productivity grid nodes
#define NN 50 // number of market penetration grid points

#define NHIST 3 // number of branches in the uncertainty tree
#define NHIST_RB 5 // number of branches in reversible Brexit case
#define TREF 4 // 2015: year in which referendum bill introduced to parliament
#define TVOTE 5 // 2016: vote takes place midway through 2016
//#define TBREXIT 13 // 2019: brexit takes place in 2019
#define TBREXIT 8 // 2019: brexit takes place in 2019
#define TREV 13 // brexit reversal: 5 years after Brexit

#define TINY 1.0e-12
#define TINYSQ 1.0e-8

#ifdef _OPENMP
#define NTH  36
#include <omp.h>
#else
#define NTH 1
#endif

#define _CONSOLE
#ifdef _CONSOLE
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define RESET "\033[0m"
#else
#define KNRM  ""
#define KRED  ""
#define KGRN  ""
#define KYEL  ""
#define KBLU  ""
#define KMAG  ""
#define KCYN  ""
#define KWHT  ""
#define RESET ""
#endif

#define SUM_3D_DIM1(v,i2,i3) ((v)[0][(i2)][(i3)] + (v)[1][(i2)][(i3)] + v[2][(i2)][(i3)])
#define SUM_4D_DIM12(v,i2,i3) ((v)[0][0][(i2)][(i3)] + (v)[0][1][(i2)][(i3)] + (v)[1][0][(i2)][(i3)] + (v)[1][1][(i2)][(i3)] + (v)[2][0][(i2)][(i3)] + (v)[2][1][(i2)][(i3)])

#define BREAK1  fprintf(logfile, KNRM "\n\t--------------------------------------------------------------------------------\n" RESET);

#define BREAK2 fprintf(logfile, KNRM "\n\n////////////////////////////////////////////////////////////////////////////////////////\n" RESET);

typedef unsigned int uint;
FILE * logfile;
uint low_exit_rate_sens;
uint higher_ch_trd_costs;
uint psi_sens;
uint zeta_sens;
uint allexport;
uint no_exporter_dyn;
uint full_mkt_pen;
uint leontief;
uint recal;
uint eval_eqm_once_flag;
uint eval_calfn_once;
uint read_seed;
uint read_stoch_seed;
uint write_seed;
uint write_stoch_seed;
uint fixl;
uint tfp_flag;
uint low_pi;
uint high_pi;
uint fin_aut;
uint comp_mkts;
uint small_gnewton_step_flag;
uint psi_uncertainty_sens;
uint reversible_brexit;
uint idio_uncertainty;
uint habit_formation;
uint high_interest_rate;
uint low_dep;

static inline void set_all_v(double * v, uint n, double x)
{
  int i;
  for(i=0; i<n; i++)
    {
      v[i] = x;
    }
}
#define SET_ALL_V(v,n,x) set_all_v( (double *)(v), (n), (x) )

static inline double sum(const double * v, uint n)
{
  uint i;
  double mysum = 0.0;
  for(i=0; i<n; i++)
    {
      mysum = mysum + v[i];
    }
  return mysum;
}
#define SUM(v,n) (sum ( (double *)(v), (n) ))

static inline double dot_prod(const double * v1, const double * v2, uint n)
{
  uint i;
  double sum = 0.0;
  for(i=0; i<n; i++)
    {
      sum = sum + v1[i]*v2[i];
    }
  return sum;
}
#define DOT_PROD(v1,v2,n) (dot_prod( (double *)(v1), (double *)(v2), (n) ))

static inline double dot_prod_ex(const double * v1, const double * v2, uint n, double ex)
{
  uint i;
  double sum = 0.0;
  for(i=0; i<n; i++)
    {
      sum = sum + pow(v1[i],ex)*v2[i];
    }
  return sum;  
}
#define DOT_PROD_EX(v1,v2,n,ex) (dot_prod_ex( (double *)(v1), (double *)(v2), (n), (ex) ))

static inline double dot_prod_ex2(const double * v1, const double * v2, uint n, double ex1, double ex2)
{
  uint i;
  double sum = 0.0;
  for(i=0; i<n; i++)
    {
      sum = sum + pow(v1[i],ex1)*pow(v2[i],ex2);
    }
  return sum;  
}
#define DOT_PROD_EX2(v1,v2,n,ex1,ex2) (dot_prod_ex2( (double *)(v1), (double *)(v2), (n), (ex1), (ex2) ))

static inline void copy_subvector(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = src[i];
    }
}
#define COPY_SUBVECTOR(d,s,n) ( copy_subvector( (double *)(d), (const double *)(s), (n) ) )

static inline void copy_subvector_log(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = log(src[i]);
    }
}
#define COPY_SUBVECTOR_LOG(d,s,n) ( copy_subvector_log( (double *)(d), (const double *)(s), (n) ) )

static inline void copy_subvector_exp(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = exp(src[i]);
    }
}
#define COPY_SUBVECTOR_EXP(d,s,n) ( copy_subvector_exp( (double *)(d), (const double *)(s), (n) ) )

static inline void write_array(double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"wb");
  uint i;    
  for(i=0; i<n; i++)
    {
      fprintf(file,"%.15f\n",v[i]);
    }    
  fclose(file);
}

static inline void write_vec_bin(const double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"wb");
  fwrite(v,n*sizeof(double),1,file);
  fclose(file);
}

static inline void write_vec_txt(const double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"wb");
  uint i;
  for(i=0; i<n; i++)
    {
      fprintf(file,"%0.15f\n",v[i]);
    }
  fclose(file);
}

static inline uint read_vec_bin(double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"rb");

  if(file)
    {
      uint bytes_read = fread(v,n*sizeof(double),1,file);
      fclose(file);
      if(bytes_read == 0)
        {   
	  return 1;
        }

      uint i;
      for(i=0; i<n; i++)
	{
	  if(!isfinite(v[i]) || isnan(v[i]))
	    {
	      fprintf(logfile,"Inf or NaN dected in binary data!\n");
	      return 1;
	    }
	}
      return 0;
    }
  else
    {
      fprintf(logfile,"File %s doesn't exist!\n",fname);
      return 1;
    }
}

static inline void linspace(double lo, double hi, uint n, double * v)
{
  double d=(hi-lo)/(n-1.0);  
  v[0]=lo;
  int i;
  for(i=1;i<n;i++)
    {
      v[i] = v[i-1]+d;
    }
  return;
}
#define LINSPACE(l,h,n,v) ( linspace((l),(h),(n),(double *)(v)) )

static inline void expspace(double lo, double hi, int n, double ex, double * v)
{
  linspace(0.0,pow(hi-lo,1.0/ex),n,v);
  int i;
  for(i=0;i<n;i++)
    {
      v[i] = pow(v[i],ex)+lo;
    }
  return;
}
#define EXPSPACE(l,h,n,v) ( expspace((l),(h),(n),(double *)(v)) )

static inline void linspace_2d(const double * lo, const double * hi, uint n1, uint n2, double * v)
{
  uint j;
  for(j=0; j<n2; j++)
    {
      double d=(hi[j]-lo[j])/(n1-1.0);
      v[0*n2+j]=lo[j];
      int i;
      for(i=1;i<n1;i++)
	{
	  v[i*n2+j] = v[(i-1)*n2+j]+d;
	}
    }
  return;
}
#define LINSPACE_2D(l,h,n,m,v) ( linspace_2d( ( (double *)l),((double *)h),(n),(m),(double *)(v) ) )

static inline void logspace(double lo, double hi, uint n, double * v)
{
  uint i;
  linspace(log(lo),log(hi),n,v);
  for(i=0; i<n; i++)
    {
      v[i] = exp(v[i]);
    }
  return;
}
#define LOGSPACE(l,h,n,v) ( logspace((l),(h),(n),(double *)(v)) )

static inline void logspace_2d(const double * lo, const double * hi, uint n1, uint n2, double * v)
{
  uint i,j;

  double llo[n2], lhi[n2];
  for(j=0; j<n2; j++)
    {
      llo[j] = log(lo[j]);
      lhi[j] = log(hi[j]);
    }
  linspace_2d(llo,lhi,n1,n2,v);
  
  for(j=0; j<n1; j++)
    {
      for(i=0; i<n2; i++)
	{
	  v[i*n2+j] = exp(v[i*n2+j]);
	}
    }
  return;
}
#define LOGSPACE_2D(l,h,n,m,v) ( linspace_2d( ( (double *)l),((double *)h),(n),(m),(double *)(v) ) )

static inline char* concat(const char *s1, const char *s2)
{
  char *result = malloc(strlen(s1)+strlen(s2)+1);
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}

static inline double enorm(const gsl_vector * f) 
{
  double e2 = 0 ;
  size_t i, n = f->size ;
  for (i = 0; i < n ; i++) {
    double fi= gsl_vector_get(f, i);
    e2 += fi * fi ;
  }
  return sqrt(e2);
}

static inline double vsum(const gsl_vector * f) 
{
  double sum = 0 ;
  size_t i, n = f->size ;
  for (i = 0; i < n ; i++) {
    double fi= gsl_vector_get(f, i);
    sum += fabs(fi) ;
  }
  return sum;
}

static inline void gdb_test_func()
{
  return;
}

// pareto distribution
static inline double pareto_cdf(double x, double kappa)
{
  return 1.0 - pow(x,-kappa);
}

static inline double pareto_pdf(double x, double kappa)
{
  return kappa*pow(x,-kappa-1.0);
}

static inline double pareto_cdf_inv(double P, double kappa)
{
  return pow(1.0-P,-1.0/kappa);
}

static inline double pareto_unconditional_moment(double kappa, double n)
{
  return kappa/(kappa-n);
}

static inline double pareto_conditional_moment_ub(double ub, double kappa, double n)
{
  return (kappa/(n-kappa))*pow(ub,n-kappa) - (kappa/(n-kappa));
}

static inline double pareto_conditional_moment_lb(double lb, double kappa, double n)
{
  return (kappa/(kappa-n))*pow(lb,n-kappa);
}


#endif
