/*
 *  pmc_rc.h
 *  ecosstat
 *
 *  Created by Karim Benabed on 01/09/08.
 *  Copyright 2008 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifndef __PMCRC_H
#define __PMCRC_H

#ifdef __PLANCK__
#include "HL2_likely/tools/readConf.h"
#include "HL2_likely/pmclib/pmc.h"
#else
#include "readConf.h"
#include "pmc.h"
#endif

#include <dlfcn.h>

#ifdef HAS_RTLD_DEFAULT
#define DLHANDLE_INIT(dlhandle) dlhandle = RTLD_DEFAULT;
#define DLHANDLE_SET(dlhandle,target) if (dlhandle != RTLD_DEFAULT) {target->dlhandle = dlhandle;}
#else
#define DLHANDLE_INIT(dlhandle) dlhandle = NULL;
#define DLHANDLE_SET(dlhandle,target) {target->dlhandle = dlhandle;}
#endif

#define INIT_FROM_RC(rca, target) {                    \
  char *funcname,*distributionname,buffer[2000];       \
  char *sofile;                                        \
  rcinit_smthng* funcinit;                             \
  void* dlhandle;                                      \
                                                       \
  funcname = rc_safeget_string(rca,"initfunc","",err); \
  forwardError(*err,__LINE__,NULL);                    \
                                                       \
  if (funcname[0]=='\0') {                             \
    distributionname = rc_get_string(rca,"name",err);  \
    forwardError(*err,__LINE__,NULL);                  \
    funcname=buffer;                                   \
    sprintf(funcname,"rcinit_%s",distributionname);    \
  }                                                    \
                                                       \
  DLHANDLE_INIT(dlhandle);                             \
                                                       \
  sofile = rc_safeget_string(rca,"libpath","",err);    \
  forwardError(*err,__LINE__,NULL);                    \
   if (sofile[0]!='\0') {                              \
    dlhandle = dlopen(sofile,RTLD_NOW | RTLD_GLOBAL);  \
    testErrorRetVA(dlhandle==NULL,dl_err,"dl error : %s",*err,__LINE__,NULL,dlerror()); \
  }                                                    \
                                                       \
  funcinit = dlsym(dlhandle,funcname);                 \
  testErrorRetVA(funcinit==NULL,dl_err,"dl error : %s",*err,__LINE__,NULL,dlerror());  \
                                                       \
  target = funcinit(rca,"",err);                       \
  forwardError(*err,__LINE__,NULL);                    \
                                                       \
  DLHANDLE_SET(dlhandle,target);                       \
}


typedef void* (rcinit_smthng)(confFile*, char*, error **err);
distribution * init_distribution_from_rc(confFile* rc,char* root, error **err);
parabox * parabox_from_rc(confFile *rc, char* root,error **err);
mix_mvdens* rc_mix_mvdens(confFile* rc, error **err);
mvdens* rc_mvdens(confFile* rc, mvdens* mv,error **err);
distribution* rcinit_mix_mvdens(confFile* rc, char* root,error **err);
pmc_simu* init_importance_from_rc(confFile *rc,char *root,error **err);
gsl_rng * rc_gsl_rng(confFile *rc, char *root,error **err);
pmc_simu* init_pmc_from_rc(confFile *rc,char *root,error **err);
char* get_itername(confFile* rc,char *key, int iter,error **err);

mcmc_run *init_mcmc_from_rc(confFile *rc,char *root,error **err);
mc_kernel * init_mc_kernel_from_rc(confFile* rc,char* root, error **err);
mc_kernel * rcinit_adaptMH(confFile* rc,char* root, error **err);

#ifdef _WITH_HDF5_
#include "hdf5.h"

void pmc_simu_hdfdump(pmc_simu*,char*,error **);
void mcmc_simu_hdfdump(mcmc_run* mcr,char* fname,error **err);

#endif

#define dl_err -101010
#endif
