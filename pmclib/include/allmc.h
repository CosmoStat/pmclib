/*
 *  allmc.h
 *  ecosstat
 *
 *  Created by Karim Benabed on 01/09/08.
 *  Copyright 2008 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifndef __ALLMC_H
#define __ALLMC_H

#include <gsl/gsl_randist.h>

#ifdef __PLANCK__
#include "HL2_likely/tools/errorlist.h"
#include "HL2_likely/pmclib/parabox.h"
#else
#include "errorlist.h"
#include "parabox.h"
#endif

struct _pmc_simu_struct_;
struct _distribution_struct_;

typedef double posterior_log_pdf_func(void *, const double *, error **);
typedef posterior_log_pdf_func log_pdf_func;

typedef void retrieve_ded_func(const void *, double *, error **);
typedef void (posterior_log_free)(void**);
typedef posterior_log_free free_func;

typedef long simulate_func(struct _pmc_simu_struct_ *, void *, gsl_rng *, parabox *, error **);
typedef void filter_func(struct _pmc_simu_struct_* , void *, error **);
typedef void update_func(void *, struct _pmc_simu_struct_ *, error **);
typedef void* mpi_exchange_func(void *, error **);
typedef double first_derivative_func(void*, int , const double*, error **err);
typedef double second_derivative_func(void*, int ,int, const double*, error **err);


#endif
