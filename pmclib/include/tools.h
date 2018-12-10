/*
 *  tools.h
 *  likely
 *
 *  Created by Karim Benabed on 13/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef __TOOLS_H
#define __TOOLS_H

#include "errorlist.h"
#include "io.h"
#include "maths_base.h"
#include "mcmc.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>

#define tls_base       -8000
#define tls_allocate   -1 + tls_base
#define tls_serialize  -2 + tls_base
#define tls_outOfBound -3 + tls_base
#define tls_badComm    -4 + tls_base
#define tls_negWeight  -5 + tls_base
#define tls_cholesky   -6 + tls_base
#define tls_negative   -7 + tls_base
#define tls_undef      -8 + tls_base
#define tls_file       -9 + tls_base
#define tls_io        -10 + tls_base
#define tls_tooManySteps -11 + tls_base
#define tls_dimension -12 + tls_base
#define tls_type      -13 + tls_base
#define tls_negHatCl  -14 + tls_base
#define tls_Nparam    -15 + tls_base
#define tls_overflow  -16 + tls_base
/* The following error causes cosmo_pmc to exceptionally stop */
#define tls_cosmo_par -17 + tls_base



/* ============================================================ *
 * Functions to print error messages.				*
 * ============================================================ */
int print_if_err(error **err, FILE *F1, FILE *F2);
void out_err_log(FILE *FLOG, const char* str, ...);

/* ============================================================ *
 * MK: Reading and writing MCM chains/PMC simulations.		*
 * ============================================================ */

ssize_t getdelim(char **lineptr, size_t *n, int delimiter, FILE *fp);
ssize_t getline(char **lineptr, size_t *n, FILE *stream);
void read_header_lc(const char *line, int npar, int n_ded, error **err);
int read_step_lc(char **line, int npar, int n_ded, double *param, double *param_ded, double *logL,
		 double *accept, error **err);
int read_next_mkmc_step(FILE *CHPRE, int npar, int n_ded, double *pstate, double *pstate_ded,
			double *logL, double *accept, error **err);
double *read_mkmc_chain(int npar, int n_ded, int nchainmax, FILE *CHACC, long *naccepted,
			fpos_t *headpos, double *param_ded, double *weight, size_t *indices, error **err);
void print_parameter(FILE *where, size_t npar, const double *params);

/* ============================================================ *
 * KB: Reading and writing MCM chains/PMC simulations.		*
 * ============================================================ */

void print_step(FILE* where,int accept, double loglkl, size_t ndim, double *params);
size_t read_chain(size_t *accepted,int accept_flag, double *loglk, double* params, int samples, int ndim, FILE* where, error **err);

/* ============================================================ *
 * To calculate the program run time.				*
 * ============================================================ */
time_t start_time(FILE *FOUT);
void end_time(time_t t_start, FILE *FOUT);

void estimate_param_covar(size_t ndim, size_t nsamples, size_t psamples, const double* params, 
                          double *pmean, double *pvar, error **err);
void estimate_param_covar_weight(size_t ndim, size_t nsamples, size_t psamples, const double* params,
				 const double *weight, double *pmean, double *pvar, error **err);


/* ============================================================ *
 * reparametrize parameters arbitrarily                 				*
 * ============================================================ */
 
typedef void repar_func_type(void*, double*,double*,error **);

typedef struct {
  int nin,nout;
  double* rpars;
  void *payload, *repar_func_data;
  posterior_log_pdf_func *payload_pdf;
  repar_func_type* repar_func;
  posterior_log_free *repar_func_free, *payload_free;
} repar;


repar* repar_init(int nin,int nout,void* payload, posterior_log_pdf_func *payload_pdf, posterior_log_free* payload_free,void* repar_func_data, repar_func_type* rfnt, posterior_log_free* repar_func_free, error **err);
double repar_lkl(void* pelf,double* pars, error **err);
void repar_free(repar** pelf);

#define REPAR_SELECT_MODE_EXPAND 1
#define REPAR_SELECT_MODE_REDUCE 2

typedef struct {
  int nin;
  int nout;
  int* select;
  double* base;
  int mode;
} repar_select;
repar_select* repar_select_init(int nin, int nout, int* select,double* base,error **err);
void repar_select_func(void* pelf, double* parin, double* parout, error **err);
void repar_select_free(void **pelf);

#endif
