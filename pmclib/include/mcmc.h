/*
 *  mcmc.h
 *  likely
 *
 *  Created by Karim Benabed on 12/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __MCMC_H
#define __MCMC_H

#include "parabox.h"

#include "pmctools/errorlist.h"
#include "pmctools/mvdens.h"
#include "allmc.h"
#include "distribution.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_int.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>

/* Errors */
#define mcmc_base       -4000
#define mcmc_allocate   -1 + mcmc_base
#define mcmc_serialize  -2 + mcmc_base
#define mcmc_outOfBound -3 + mcmc_base
#define mcmc_badComm    -4 + mcmc_base
#define mcmc_negWeight  -5 + mcmc_base
#define mcmc_cholesky   -6 + mcmc_base
#define mcmc_negative   -7 + mcmc_base
#define mcmc_undef      -8 + mcmc_base
#define mcmc_file       -9 + mcmc_base
#define mcmc_io        -10 + mcmc_base
#define mcmc_tooManySteps -11 + mcmc_base
#define mcmc_dimension -12 + mcmc_base
#define mcmc_type      -13 + mcmc_base
#define mcmc_negHatCl  -14 + mcmc_base
#define mcmc_infnan    -15 + mcmc_base
#define mcmc_unknown   -16 + mcmc_base
#define mcmc_prior     -17 + mcmc_base
#define mcmc_ded       -18 + mcmc_base

#define MC_AF_ALL      -1
#define MC_AF_ACCEPT    1
#define MC_AF_REJECT    0
#define MC_AF_OUTOFBOX -2

#define MC_INF 1e99

struct mcmc_run_object;
struct mcmc_kernel_object;
 
typedef void transition_func(double*, void*, gsl_rng*, error**);
typedef double transition_ratio_func(double*,double*, void*, error**);
typedef void adapt_func(struct mcmc_kernel_object *,struct mcmc_run_object*, error **);

typedef  struct mcmc_kernel_object {
	size_t ndim;
	transition_func *sample;
	transition_ratio_func *transition_ratio; 
	void* extra;
  posterior_log_free *free;
  adapt_func *adapt;
  int nadapt;
  int adaptifaccepted;
  void* dlhandle;
  
} mc_law;
typedef mc_law mc_kernel;

mc_law * mc_law_init(size_t ndim, transition_func *sample, transition_ratio_func *log_pdf, void* extra, error ** err);
void mc_law_free(mc_law** self, void (*extra_free)(void**));

int mcmc_step(double *pstate,double* nstate, double *log_pdf_value, posterior_log_pdf_func *pdf_func, 
	      void *lkl_extra, mc_law* law, parabox *pb,gsl_rng* rng,error **err);

mc_law *mc_mvdens_init(mvdens* mself,error **err);
void mc_mvdens_free(mc_law** self);
void mc_mvdens_sample(double* state, void* self, gsl_rng* rng, error** err);
double trivial_transition_func(double* pstate,double* nstate, void* self, error** err);

void estimate_param_covar(size_t ndim, size_t nsamples, size_t psamples, const double* params, 
                          double *pmean, double *pvar, error **err);
void estimate_param_covar_weight(size_t ndim, size_t nsamples, size_t psamples, const double* params,
				 const double *weight, double *pmean, double *pvar, error **err);


typedef struct mcmc_run_object {
  mc_kernel *kernel;
  distribution *target;
  parabox* pb;
  
  int ndim;
  int nded;
  int ntot;
  
  int nsample;
  int current;
  int naccepted;
  
  int nbatch;
  int lastbatch;
  
  int nadapt;
  int adaptifaccepted;
  int cntadapt;
  int lastadapt;
  
  double* pars;
  double* deds;
  double* loglkl;
} mcmc_run;

mcmc_run *init_mcmc_run(int nsample, int nbatch, mc_law* kernel, distribution *target, double* pars0, parabox * pb, error **err);
int mcmc_run_batch(mcmc_run *mcr, gsl_rng *r, error **err);
void free_mcmc_run(void** pmcr);
void free_mc_kernel(void** pmcl);
mc_kernel *mc_kernel_init(size_t ndim, transition_func *sample,transition_ratio_func *log_pdf, void* data, posterior_log_free* free, int nadapt, int adaptifaccept, adapt_func *adapt, error ** err);

typedef struct  {
  mvdens *gauss;
  double k_damp,c_enlarge;
  double *buf,*mu_t,*sig_t;
}  adaptMH;
  
mc_kernel * adaptMH_init(int ndim, double *sig0, int nadapt, int adaptifaccepted, double c_enlarge, double k_damp,error **err);
void adaptMH_sample(double* state, void* self, gsl_rng* rng, error** err);
void adaptMH_free(void** pmh);
void adaptMH_adapt(mc_kernel *kernel, mcmc_run* mc, error **err);

#endif

