/*
 *  distribution.h
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 04/11/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


#ifndef __MC_DIST
#define __MC_DIST
#include <gsl/gsl_randist.h>

#include "pmctools/errorlist.h"
#include "pmctools/mvdens.h"
#include "pmctools/maths.h"

#include "parabox.h"
#include "allmc.h"
#include <dlfcn.h>

 typedef char _char_name[1024];

typedef struct _distribution_struct_ {
  int ndim, n_ded,ndef;
  double *pars;
  int * def;
  void* data;
  posterior_log_pdf_func* log_pdf;
  retrieve_ded_func* retrieve;
  posterior_log_free* free;
  simulate_func *simulate;
  mpi_exchange_func *broadcast_mpi;
  first_derivative_func *f_der;
  second_derivative_func *d_der;
  
  void* dlhandle;
  _char_name *name;
   
} distribution;

distribution* init_distribution_full(int ndim,
                                     void* data, 
                                     posterior_log_pdf_func* log_pdf,
                                     posterior_log_free* freef,
                                     simulate_func *simulate,
                                     int nded,
                                     retrieve_ded_func* retrieve,
                                     error **err);

distribution* init_simple_distribution(int ndim,
                               void* data, 
                               posterior_log_pdf_func* log_pdf,
                               posterior_log_free* freef,
                               error **err);

distribution* init_distribution(int ndim,
                                void* data, 
                                posterior_log_pdf_func* log_pdf,
                                posterior_log_free* freef,
                                simulate_func *simulate,
                                error **err);
void free_distribution(distribution **pdist) ;

int distribution_get_name(distribution *dist,char* name,error **err);
int* distribution_get_names(distribution *dist,int nname,char** name,int includeded,error **err);
void distribution_set_names(distribution *dist,char** name, error **err);

void distribution_set_broadcast(distribution* dist, mpi_exchange_func* broadcast, error **err);


double distribution_lkl(void* pdist, const double* pars, error **err);
#define distribution_log_pdf distribution_lkl

void distribution_retrieve(const void* pdist, double* pars, error **err);

void distribution_set_default(distribution *dist, int ndef, int* idef, double* vdef,error **err);
void distribution_set_default_name(distribution *dist, int ndef, char** idef, double* vdef,error **err);

void distribution_set_default(distribution *dist, int ndef, int* idef, double* vdef,error **err);

void distribution_set_derivative_func(distribution *dist, first_derivative_func* f_der, second_derivative_func* d_der, error **err);
double distribution_first_derivative(distribution *dist, int idir, double* par, error **err);
double distribution_second_derivative(distribution *dist, int idir, int jdir, double* par, error **err);
double* distribution_second_derivative_matrix(distribution *dist, double *par, error **err);
double distribution_deriv_along(distribution *dist, double *pars, double *along, error **err);

distribution * combine_distribution_init(int ndim, int nded, error **err);
distribution * combine_distribution_simple_init(int ndim, error **err);

double combine_lkl(void *pcbd, const double* pars, error **err);
void combine_retrieve(const void *pcbd, double* pded, error **err);
void add_to_combine_distribution(distribution *comb, distribution *addon, int *dim_idx, int *ded_idx, error **err);

void add_to_combine_distribution_name(distribution *comb, distribution *addon, error **err);


void combine_free(void **pcbd);

typedef struct {
  double *pars,*pded,*dummy;
  int *ded_from, **dim,**ded;
  distribution **dist;
  int ndim,nded,ndist,ndummy;
} comb_dist_data;

distribution *add_gaussian_prior(distribution *orig, int ndim, int *idim, double* loc, double *var, error **err);
distribution *add_gaussian_prior_2(distribution *orig, int ndim, int *idim, double* loc, double *var, error **err);
distribution *add_gaussian_prior_name(distribution *orig, int ndim, char**idim, double* loc, double *var, error **err);
distribution *add_gaussian_prior_2_name(distribution *orig, int ndim, char**idim, double* loc, double *var, error **err);



#define dist_base       -6700
#define dist_undef       -1 + dist_base 
#define dist_type        -2 + dist_base 

#endif
