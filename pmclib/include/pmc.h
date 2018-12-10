/*
 *  pmc.h
 *  likely
 *
 *  Created by Karim Benabed on 10/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PMC_H
#define __PMC_H

#include "errorlist.h"
#include "mvdens.h"
#include "parabox.h"
#include "tools.h"
#include "distribution.h"
#include "allmc.h"

#define _PRSTP_ 20. 

/* errors */
#define pmc_base       -6000
#define pmc_allocate   -1 + pmc_base
#define pmc_serialize  -2 + pmc_base
#define pmc_outOfBound -3 + pmc_base
#define pmc_badComm    -4 + pmc_base
#define pmc_negWeight  -5 + pmc_base
#define pmc_cholesky   -6 + pmc_base
#define pmc_negative   -7 + pmc_base
#define pmc_undef      -8 + pmc_base
#define pmc_file       -9 + pmc_base
#define pmc_io        -10 + pmc_base
#define pmc_tooManySteps -11 + pmc_base
#define pmc_dimension -12 + pmc_base
#define pmc_type      -13 + pmc_base
#define pmc_negHatCl  -14 + pmc_base
#define pmc_infnan    -15 + pmc_base
#define pmc_incompat  -16 + pmc_base
#define pmc_nosamplep -17 + pmc_base
#define pmc_sort      -18 + pmc_base
#define pmc_infinite  -19 + pmc_base
#define pmc_isLog     -20 + pmc_base

#define PI         3.141592653589793
#define LOGSQRT2PI 0.918938533204673
#define DYNMAX 500.0


/* Output of parameters in case an error occured */
#define ParameterError(err,par,ndim) ParameterErrorVerb(err,par,0,ndim)
#define ParameterErrorVerb(err,par,quiet,ndim) {	\
     int j;						\
     if (isError(err)) {					\
	if (! quiet) { \
	   fprintf(stderr, "*** Error, parameter ***\n");		\
	   for (j=0; j<ndim; j++) {					\
	      fprintf(stderr, " %10g", par[j]);				\
	   }								\
	   fprintf(stderr, "\n*** -- ***\n");				\
	   printError(stderr,err);					\
	} \
	purgeError(&(err));						\
	continue;							\
     }									\
}


/* Minimum number of sample points sampled from one component required for that *
 * component to be kept alive.						        */
#define MINCOUNT 20

#define MC_NORM 1
#define MC_UNORM 0

typedef struct _pmc_simu_struct_ {
  // data section
  long nsamples;
  int ndim, n_ded;
  void *data;
  double *X, *X_ded, *weights;
  double *log_rho;
  short *flg;
  size_t *indices;
  int isLog;
  double logSum;
  double maxW,maxR;
  
  // parabox 
  parabox *pb;
  
  // target specific section
  distribution *target;
  
  // importance specific section
  distribution *proposal;
  int prop_print_step;
  
  // pmc specific section
  void* filter_data;
  filter_func *pmc_filter;
  update_func* pmc_update;
  int retry;
  
  //mpi section
  int mpi_rank;
  int mpi_size;
} pmc_simu;

typedef struct {
  long nsamples;
  int ndim,n_ded;
  int isLog;
  double logSum;
} pmc_simu_head;


// high level interface
pmc_simu* pmc_simu_init(long nsample, int ndim, error **err);
pmc_simu* pmc_simu_init_plus_ded(long nsample, int ndim, int n_ded, error **err);
void pmc_simu_free(pmc_simu **self);
void pmc_simu_init_target(pmc_simu* psim, 
                          distribution* target,
                          parabox *pb,error **err);
void pmc_simu_init_proposal(pmc_simu* psim, 
                            distribution *proposal,
                            int prop_print_step,error **err);

void pmc_simu_init_pmc(pmc_simu* psim, 
                       void* filter_data,
                       filter_func* pmc_filter,
                       update_func* pmc_update);
size_t pmc_simu_importance(pmc_simu *psim, gsl_rng *r,error **err);
double pmc_simu_pmc_step(pmc_simu *psim, gsl_rng *r,error **err);
void pmc_simu_init_classic_importance(pmc_simu* psim,distribution* target, parabox *pb,
                                      void *prop_data,int prop_print_step,error **err);
void pmc_simu_init_classic_pmc(pmc_simu* psim,distribution* target, parabox *pb,
                               void *prop_data,int prop_print_step,void* filter_data,
                               filter_func* pmc_filter,error **err);



void out_pmc_simu(const char *name, const pmc_simu *psim, error **err);
void pmc_simu_realloc(pmc_simu *psim,long newsamples,error **err);
void pmc_simu_print(FILE *where, pmc_simu *psim,double nrm);
void pmc_simu_dump(FILE *where,pmc_simu *psim,error **err);
size_t *sample_from_pmc_simu(const pmc_simu *psim, const gsl_rng * rng, error **err);

long simulate_mix_mvdens_void(struct _pmc_simu_struct_ *psim, void *proposal, gsl_rng * r,
			      parabox *pb, error **err);
long simulate_mix_mvdens(pmc_simu *psim, mix_mvdens *m, gsl_rng * r,parabox *pb,error **err);

size_t get_importance_weight(pmc_simu *psim, mix_mvdens *m,
                                 posterior_log_pdf_func *posterior_log_pdf, 
                                 void *extra, error **err);
size_t generic_get_importance_weight(pmc_simu *psim, void *m,
                             posterior_log_pdf_func *proposal_log_pdf,
                             posterior_log_pdf_func *posterior_log_pdf,
                             void *extra, error **err);
size_t get_importance_weight_plus_ded(pmc_simu *psim, mix_mvdens *m,
				      posterior_log_pdf_func *posterior_log_pdf, 
				      retrieve_ded_func *retrieve_ded, void *extra, error **err);
size_t generic_get_importance_weight_and_deduced(pmc_simu *psim, void *m,
                             posterior_log_pdf_func *proposal_log_pdf,
                             posterior_log_pdf_func *posterior_log_pdf,
                             retrieve_ded_func *retrieve_ded, void *extra, error **err);
size_t generic_get_importance_weight_and_deduced_verb(pmc_simu *psim, void *proposal_data,
						      posterior_log_pdf_func *proposal_log_pdf,
						      posterior_log_pdf_func *posterior_log_pdf,
						      retrieve_ded_func *retrieve_ded, void *target_data, double t_iter_tempered,
						      int quiet, error **err);


void filter_importance_weight(pmc_simu *psim,double nsig,error **err);

double normalize_importance_weight(pmc_simu *psim, error **err);

void update_prop(mix_mvdens *m,pmc_simu *psim, error **err);
void update_prop_rb_void(void *m, pmc_simu *psim, error **err);
void update_prop_rb(mix_mvdens *m, pmc_simu *psim, error **err);
void cleanup_after_update(mix_mvdens *m,size_t *count,double minwgt,pmc_simu *psim,double *buf,error **err);
void prepare_update(mix_mvdens *m, size_t *count,pmc_simu *psim,double **pbuf,error **err);

double perplexity(pmc_simu* psim,int normalize, error **err);
double perplexity_and_ess(pmc_simu* psim,int normalize,double * ess, error **err);
double evidence(pmc_simu *psim, double *ln_evi,error **err);

double mean_from_psim(double *X, double *w, short *flg, int nsamples, int ndim, int a);
int double_cmp(const void *av, const void *bv);
void clip_weights(pmc_simu *psim, int n, FILE *OUT, error **err);

pmc_simu *pmc_simu_from_file(FILE *PMCSIM, int this_nsamples, int npar, int n_ded,
			     mix_mvdens *proposal, int nclipw, error **err);
double fill_read_pmc_simu(pmc_simu *psim, mix_mvdens *proposal, error **err);

unsigned int *resample_residual(const pmc_simu *psim, const gsl_rng *rng, error **err);


size_t* sort_indices(double *list, size_t size, size_t *in_ind,size_t *buf,int order,error **err);

int* filter_histogram_init(int nbins,int clipleft,int clipright, error **err);
void filter_histogram(pmc_simu *psim ,void* filter_data, error **err);

long simulate_mix_mvdens(struct _pmc_simu_struct_ *psim,mix_mvdens *proposal, gsl_rng * r,parabox *pb,error **err);
distribution* mix_mvdens_distribution(int ndim, void *mxmv, error **err);

#endif
