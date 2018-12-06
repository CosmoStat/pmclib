/*
 *  optimize.h
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 09/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#undef LBFGS_FLOAT
#define LBFGS_FLOAT 64

#ifdef __PLANCK__
#include "HL2_likely/pmclib/montecarlo2.h"
#include "HL2_likely/tools/maths.h"
#undef epsilon
#if 0
#include "src/lbfgs/lbfgs.h"
#endif
#else
#include "distribution.h"
#include "pmc.h"
#include "maths.h"
#undef epsilon
#if 0
#include "lbfgs.h"
#endif
#endif



#ifndef _PMC_OPTIMIZE_
#define _PMC_OPTIMIZE_

#if 0
typedef struct {
  distribution* distrib;
  error *err;
  lbfgs_parameter_t params;
  double *data;
  double *x,*df,*scales,*pars,*df0,*x0,*hess,*bs;
} lbfgs_payload;

lbfgsfloatval_t lbfgs_evaluate(void *instance,
                         lbfgsfloatval_t *x,
                         lbfgsfloatval_t *g,
                         const int n,
                         const lbfgsfloatval_t step
                         );

int lbfgs_progress(
                   void *instance,
                   const lbfgsfloatval_t *x,
                   const lbfgsfloatval_t *g,
                   const lbfgsfloatval_t fx,
                   const lbfgsfloatval_t xnorm,
                   const lbfgsfloatval_t gnorm,
                   const lbfgsfloatval_t step,
                   int n,
                   int k,
                   int ls
                   );

lbfgs_payload *lbfgs_payload_init(distribution *dist, double *pars, double* scales,error **err);
void lbfgs_payload_free(void** dist);
void lbfgs_optimize(lbfgs_payload *pay,error **err);
char* lbfgs_error_string(int lbfgs_err);
double* lbfgs_get_result(lbfgs_payload *pay,error **err);
double* lbfgs_get_variance(lbfgs_payload *pay,error **err);
void lbfgs_set_linesearch(lbfgs_payload *pay,int ls, error **err);
void lbfgs_set_minmaxstep(lbfgs_payload *pay,double mins, double maxs, error **err);
void lbfgs_set_convergence(lbfgs_payload *pay, double epsilon, int past, double delta, error **err);
#endif

// high level

typedef struct _opts optimize_struct;
typedef void optimize_func(optimize_struct*, error **err);
typedef double*  optimize_get_func(optimize_struct*, error **err);

struct _opts {
  int ndim;
  distribution *target;
  void *data;
  optimize_func *optimize;
  optimize_get_func *get_result;
  optimize_get_func *get_variance;
  posterior_log_free *free;
  int cvg;
  void* dlhandle;
};
  
optimize_struct * optimize_init(distribution *target,void *data ,optimize_func *opt, optimize_get_func *get_p, optimize_get_func *get_v,posterior_log_free *fropt, error **err);
void optimize(optimize_struct *opt,error **err);
double* optimize_get_best_pars(optimize_struct *opt, error **err);
double* optimize_get_variance(optimize_struct *opt, error **err);
void optimize_free(void **opt);

#ifdef _WITH_RC_
#include "pmc_rc.h"
#if 0
void *rcinit_lbfgs(confFile *rc, char * root, error **err);
#endif
void *init_optimize_from_rc(confFile *rc, char * root, error **err);
void *rcinit_bydir(confFile *rc, char * root, error **err);
void *rcinit_bfgs(confFile *rc, char * root, error **err);
void *rcinit_secant_linesearch(confFile *rc, char * root, error **err);
void *init_linesearch_from_rc(confFile *rc, char * root, error **err);
#endif

typedef struct _lstr linesearch_struct;

typedef double ls_func(linesearch_struct *lst, distribution *dist, double *PIM, double *pk, double *gk, error **err);

struct _lstr {
  int ndim;
  void *data;
  ls_func * ls_fnc;
  double *PIM;
  posterior_log_free *ls_free;
  void* dlhandle;
}; 
  
typedef struct {
  double step;
  int maxstep;
  double tol;
  int ndim;
  int golden;
} linesearch_secant_struct;

typedef struct {
  double tol;
  int maxstep;
  double *PIM,*pk,*scale;
  linesearch_struct *ls;
  distribution *dist;
} bydir_optimize_struct;

linesearch_struct* linesearch_init(int ndim, void* data, ls_func* fnc, posterior_log_free * frr, error **err);
double linesearch(linesearch_struct *ls, distribution *dist, double *PIM, double *pk, double *gk,error **err);
void free_linesearch(void **pls);

linesearch_struct* secant_linesearch_init(int ndim, double step, int maxstep, double tol, error **err);
void secant_free(void** psc);
void addalphap(int ndim, double *PIM, double alpha, double* pk);
double deriv_along(distribution *dist, double *PIM, double *pk, double *gk,int update,error **err);
double secant_linesearch(linesearch_struct *lst, distribution *dist, double *PIM, double *pk, double *gk,error **err);

optimize_struct* bydir_init(distribution *dist, double* guess, double* scale,linesearch_struct *ls, double tol, int maxstep, error **err);
void bydir_func(optimize_struct *opt, error **err);
void bydir_free(void **pbd);
double *bydir_get_p(optimize_struct *opt, error **err);

#ifdef _WITH_LAPACK_
typedef struct {
  double *bk,*bk_buf,*gk,*sk,*pk,*PIM,*yk,*scale;
  distribution *dist;
  int maxstep;
  double tol;
  linesearch_struct *ls;
} bfgs_optimize_struct;
optimize_struct* bfgs_init(distribution *dist, double* guess, double* scale,linesearch_struct *ls, double tol, int maxstep, error **err);
void bfgs_func(optimize_struct *opt, error **err);
void bfgs_free(void **pbd);
double *bfgs_get_p(optimize_struct *opt, error **err);
#endif

#endif