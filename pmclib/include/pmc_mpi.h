/*
 *  pmc_mpi.h
 *  likely
 *
 *  Created by Karim Benabed on 13/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PMC_MPI_H
#define __PMC_MPI_H

#include <mpi.h>

#ifdef __PLANCK__
#include "HL2_likely/pmclib/pmc.h"
#include "HL2_likely/tools/mvdens.h"
#include "HL2_likely/tools/errorlist.h"
#else
#include "pmc.h"
#include "mvdens.h"
#include "errorlist.h"
#endif

/* Set DEBUG_MPI to 1 for debug outputs concerning communication between master and slaves */
#define DEBUG_MPI 0
#define fprintfDEBUG(file,str,...) if (DEBUG_MPI) fprintf(file, str, __VA_ARGS__)


#undef isErrorExit

#define isErrorExit(err,F) {                      \
  if(isError(err)) {                              \
    printError(F,err);                            \
    MPI_Abort(MPI_COMM_WORLD,getErrorValue(err)); \
  }                                               \
}

#define exitAllOnError(myid,err) exitOnError(err,stderr);

/* communication tags */
#define mc_tag_data_size          1
#define mc_tag_nok                2
#define mc_tag_data               3
#define mc_tag_flag               4
#define mc_tag_mix_mvdens_size    5
#define mc_tag_mix_mvdens_buffer  6
#define mc_tag_simulation         7
#define mc_tag_simulation_size    8
#define mc_tag_flg                9
#define mc_tag_mix_mvt_size       10
#define mc_tag_mix_mvt_buffer     11
#define mc_tag_logrho             12
#define mc_tag_maxw               13
#define mc_tag_maxr               14
#define mc_tag_data_ded           15

pmc_simu* pmc_simu_init_mpi(long nsamples, int ndim, int nded,error **err);

size_t pmc_simu_importance_mpi(pmc_simu *psim, gsl_rng *r,error **err);
double pmc_simu_pmc_step_mpi(pmc_simu *psim, gsl_rng *r,error **err);
void* mvdens_broadcast_mpi(void*,error **err);
void pmc_simu_realloc_mpi(pmc_simu *psim,long newsamples,error **err);

void pmc_simu_init_classic_importance_mpi(pmc_simu* psim,
                                          distribution *target, parabox *pb,
                                          void *prop_data,int prop_print_step,error **err);
void pmc_simu_init_classic_pmc_mpi(pmc_simu* psim,
                                   distribution *target, parabox *pb,
                                   void *prop_data,int prop_print_step,void* filter_data,
                                   filter_func* pmc_filter,error **err);
distribution* init_distribution_full_mpi(int ndim,
                                         void* data, 
                                         posterior_log_pdf_func* log_pdf,
                                         posterior_log_free* freef,
                                         simulate_func *simulate,
                                         int nded,
                                         retrieve_ded_func* retrieve,
                                         mpi_exchange_func *broadcast_mpi, error **err);

distribution* init_distribution_mpi(int ndim,
                                    void* data, 
                                    posterior_log_pdf_func* log_pdf,
                                    posterior_log_free* freef,
                                    simulate_func *simulate,
                                    mpi_exchange_func *broadcast_mpi, error **err);
void distribution_broadcast(distribution *dist,error **err);

distribution *mix_mvdens_distribution_mpi(int ndim,void *mxmv, error **err);

size_t send_simulation(pmc_simu *psim, int nproc, error **err);
void receive_simulation(pmc_simu *psim, int basetag,int myid, error **err);

void send_mix_mvdens(mix_mvdens *m, int nproc,error **err);
mix_mvdens *receive_mix_mvdens(int myid, int nproc,error **err);

void send_importance_weight(int myid, int basetag, pmc_simu *psim,size_t nok,error **err);
size_t receive_importance_weight(pmc_simu *psim,int nproc,int master_cnt,int 
                                 master_sample,error **err);

#ifdef _WITH_RC_               

#include "readConf.h"                  
#include "pmc_rc.h"                  
mix_mvdens* rc_mix_mvdens_mpi(confFile* rc, error **err);
distribution* rcinit_mix_mvdens_mpi(confFile* rc, error **err);                               
pmc_simu* init_importance_mpi_from_rc(confFile* rc,char* root,error **err);
pmc_simu* init_pmc_mpi_from_rc(confFile *rc,char *root,error **err);
#endif

#endif
