/*
 *  pmc_rc.c
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 04/11/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifdef __PLANCK__
#include "HL2_likely/pmclib/pmc_rc.h"
#else
#include "pmc_rc.h"
#endif

int _rc_dist_set_names(confFile *rca, distribution* target,error **err) {
  int hn;
  int ndum;
  
  hn = rc_has_key(rca,"parameter_name",err);
  forwardError(*err,__LINE__,0);
  if (hn==1) {
    char **strname;
    
    ndum = rc_get_string_array(rca,"parameter_name",&strname,err);
    forwardError(*err,__LINE__,0);

    testErrorRetVA(ndum!=target->ndim+target->n_ded,dist_undef,"not enough parameter names (got %d expected %d)",*err,__LINE__,NULL,ndum,target->ndim+target->n_ded);
    
    distribution_set_names(target,strname,err);
    forwardError(*err,__LINE__,0);
  }
  return hn;
}
distribution * init_distribution_from_rc(confFile* rc,char* root, error **err) {
  distribution *target;
  confFile *rca;
  int hk,ndef,ndum,*idef,i,hn;
  char **cdef;
  double *vdef;

  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);

  INIT_FROM_RC(rca,target);
  
  hn = _rc_dist_set_names(rca, target,err);
  forwardError(*err,__LINE__,NULL);
  
  hk = rc_has_key(rca,"conditional",err);  
  forwardError(*err,__LINE__,NULL);
  
  if (hk==1) {
    ndum = rc_get_real_array(rca,"conditional.value",&vdef,err);
    forwardError(*err,__LINE__,NULL);
    
    ndef = rc_get_integer_array_asint(rca,"conditional.id",&idef,err);
    if (isError(*err)) {
      _DEBUGHERE_("","");
      printError(stderr,*err);
      purgeError(err);
      ndef = rc_get_string_array(rca,"conditional.id",&cdef,err);
      forwardError(*err,__LINE__,NULL);      
      testErrorRetVA(ndum!=ndef,dist_base,
                     "conditional.id and conditional.value have different sizes (got %d and %d)",
                     *err,__LINE__,NULL,ndef,ndum);
      distribution_set_default_name(target, ndef, cdef, vdef,err);
      forwardError(*err,__LINE__,NULL);
    } else {
      testErrorRetVA(ndum!=ndef,dist_base,
                     "conditional.id and conditional.value have different sizes (got %d and %d)",
                     *err,__LINE__,NULL,ndef,ndum);
      distribution_set_default(target, ndef, idef, vdef,err);
      forwardError(*err,__LINE__,NULL);
  
    }
  }
  
  rc_close(&rca);
  return target;
}

void * rcinit_combine_distributions(confFile* rc,char* root, error **err) {
  confFile *rca;
  int ndim,nded,ndist,i,j,isar,hn;
  distribution *comb,*addon;
  int *dim_idx,*ded_idx;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  ndim = rc_get_integer(rca,"ndim",err);
  forwardError(*err,__LINE__,NULL);
  
  nded = rc_safeget_integer(rca,"nded",0,err);
  forwardError(*err,__LINE__,NULL);
  
  comb = combine_distribution_init(ndim,nded,err);
  forwardError(*err,__LINE__,NULL);
  
  hn = _rc_dist_set_names(rca, comb,err);
  forwardError(*err,__LINE__,NULL);
  
  ndist = rc_array_size(rca,"distributions",err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<ndist;i++) {
    char dname[500];
    long *ldim,*lded;
    int ndum;
    int doname;
    
    sprintf(dname,"distributions[%d].distribution",i+1);
    addon = init_distribution_from_rc(rca,dname,err);
    forwardError(*err,__LINE__,NULL);
    
    dim_idx=NULL;
    sprintf(dname,"distributions[%d].dim_idx",i+1);
    isar = rc_has_key(rca,dname,err);
    forwardError(*err,__LINE__,NULL);
    doname=1;
    if (isar==1) {
      doname=0;
      dim_idx = malloc_err((addon->ndim)*sizeof(int),err);
      forwardError(*err,__LINE__,NULL);
      ndum = rc_get_integer_array(rca,dname,&ldim,err);
      forwardError(*err,__LINE__,NULL);
      testErrorRetVA(ndum!=addon->ndim,dist_base,"bad size for %s (got %d expected %d)",
                     *err,__LINE__,NULL,dname,ndum,addon->ndim);
      for(j=0;j<ndum;j++) {
        dim_idx[j] = ldim[j];
      }  
    }
    
    ded_idx=NULL;
    sprintf(dname,"distributions[%d].ded_idx",i+1);
    isar = rc_has_key(rca,dname,err);
    forwardError(*err,__LINE__,NULL);
    if (isar==1) {
      doname=0;
      ded_idx = malloc_err((addon->n_ded)*sizeof(int),err);
      forwardError(*err,__LINE__,NULL);
      ndum = rc_get_integer_array(rca,dname,&lded,err);
      forwardError(*err,__LINE__,NULL);
      testErrorRetVA(ndum!=addon->n_ded,dist_base,"bad size for %s (got %d expected %d)",
                     *err,__LINE__,NULL,dname,ndum,addon->n_ded);
      for(j=0;j<ndum;j++) {
        ded_idx[j] = lded[j];
      }  
    }
    
    if (doname==1 && hn==1) {
      add_to_combine_distribution_name(comb,addon,err);
      forwardError(*err,__LINE__,NULL);      
    }
    add_to_combine_distribution(comb, addon, dim_idx, ded_idx, err);
    forwardError(*err,__LINE__,NULL);

    if (dim_idx!=NULL) {
      free(dim_idx);  
    }
    if (ded_idx!=NULL) {
      free(ded_idx);  
    }

  }
  rc_close(&rca);
  return comb;
}


void * rcinit_reparametrise_distribution(confFile* rc,char* root, error **err) {
  confFile *rca;
  int ndim,nded,ndist,i,j,isar,ndum;
  distribution *comb,*addon;
  int *dim_idx,*ded_idx;
  long *ldim,*lded;

  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  addon = init_distribution_from_rc(rca,"distribution",err);
  forwardError(*err,__LINE__,NULL);

  comb = combine_distribution_init(addon->ndim,addon->n_ded,err);
  forwardError(*err,__LINE__,NULL);

  dim_idx=NULL;
  isar = rc_has_key(rca,"dim_idx",err);
  forwardError(*err,__LINE__,NULL);
  if (isar==1) {
    dim_idx = malloc_err((addon->ndim)*sizeof(int),err);
    forwardError(*err,__LINE__,NULL);
    ndum = rc_get_integer_array(rca,"dim_idx",&ldim,err);
    forwardError(*err,__LINE__,NULL);
    testErrorRetVA(ndum!=addon->ndim,dist_base,"bad size for %s (got %d expected %d)",
                   *err,__LINE__,NULL,"dim_idx",ndum,addon->ndim);
    for(j=0;j<ndum;j++) {
      dim_idx[j] = ldim[j];
    }  
  }
  
  ded_idx=NULL;
  isar = rc_has_key(rca,"ded_idx",err);
  forwardError(*err,__LINE__,NULL);
  if (isar==1) {
    ded_idx = malloc_err((addon->n_ded)*sizeof(int),err);
    forwardError(*err,__LINE__,NULL);
    ndum = rc_get_integer_array(rca,"ded_idx",&lded,err);
    forwardError(*err,__LINE__,NULL);
    testErrorRetVA(ndum!=addon->n_ded,dist_base,"bad size for %s (got %d expected %d)",
                   *err,__LINE__,NULL,"ded_idx",ndum,addon->n_ded);
    for(j=0;j<ndum;j++) {
      ded_idx[j] = lded[j];
    }  
  }
  
  add_to_combine_distribution(comb, addon, dim_idx, ded_idx, err);
  forwardError(*err,__LINE__,NULL);

  if (dim_idx!=NULL) {
    free(dim_idx);  
  }
  if (ded_idx!=NULL) {
    free(ded_idx);  
  }

  rc_close(&rca);
  return comb;
}

void * rcinit_add_gaussian_prior(confFile* rc,char* root, error **err) {
  confFile *rca;
  int ndim,i,isar,ndum;
  distribution *comb,*addon;
  int *dim_idx;
  long *ldim;
  double *mn,*var;
  int hn;
  char **strname;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  
  addon = init_distribution_from_rc(rca,"distribution",err);
  forwardError(*err,__LINE__,NULL);
  
  // get the mean
  ndim = rc_get_real_array(rca,"mean",&mn,err);
  forwardError(*err,__LINE__,NULL);
  
  // get the var
  rc_get_real_assquarematrix(rca,"var",&var,ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  dim_idx = NULL;
  
  hn = rc_has_key(rca,"parameter_name",err);
  forwardError(*err,__LINE__,0);
  
  if (hn==1) {
    int ndum;
    ndum = rc_get_string_array(rca,"parameter_name",&strname,err);
    forwardError(*err,__LINE__,0);

    testErrorRetVA(ndum!=ndum,dist_undef,"not enough parameter names (got %d expected %d)",*err,__LINE__,NULL,ndum,ndim);

    comb = add_gaussian_prior_2_name(addon, ndim, strname,mn, var, err);    
    forwardError(*err,__LINE__,NULL);
    rc_close(&rca);

    return comb;
  }
   
  isar = rc_has_key(rca,"prior_idx",err);
  forwardError(*err,__LINE__,NULL);
  if (isar==1) {
    
    ndum = rc_get_integer_array(rca,"prior_idx",&ldim,err);
    forwardError(*err,__LINE__,NULL);
    testErrorRetVA(ndum!=ndim,dist_base,"bad size for %s (got %d expected %d)",
                   *err,__LINE__,NULL,"prior_idx",ndum,ndim);
  
    dim_idx = malloc_err(sizeof(int)*ndum,err);
    forwardError(*err,__LINE__,NULL);
           
    for(i=0;i<ndum;i++) {
      dim_idx[i] = ldim[i];
    }
  }
  
  comb = add_gaussian_prior_2(addon, ndim, dim_idx,mn, var, err);
  forwardError(*err,__LINE__,NULL);

  if(isar==1) {
    free(dim_idx);
  }
  
  rc_close(&rca);
  
  return comb;
  
}

mc_kernel * init_mc_kernel_from_rc(confFile* rc,char* root, error **err) {
  mc_kernel *kernel;
  confFile *rca;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);

  INIT_FROM_RC(rca,kernel);
  
  rc_close(&rca);
  return kernel;
}

mc_kernel * rcinit_adaptMH(confFile* rc,char* root, error **err) {
  mc_kernel *kernel;
  confFile *rca;
  int ndim,nadapt,adaptifaccepted;
  double c_enlarge,k_damp;
  double *sig2;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);

  ndim = rc_get_integer(rca,"ndim",err);
  forwardError(*err,__LINE__,NULL);
  
  rc_get_real_assquarematrix(rca,"sig0",&sig2,ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  nadapt = rc_safeget_integer(rca,"nadapt",0,err);
  forwardError(*err,__LINE__,NULL);
  
  adaptifaccepted = rc_safeget_integer(rca,"adaptifaccepted",1,err);
  forwardError(*err,__LINE__,NULL);
  
  k_damp = rc_safeget_real(rca,"k_damp",0,err);
  forwardError(*err,__LINE__,NULL);
  
  c_enlarge = rc_safeget_real(rca,"c_enlarge",2.38,err);
  forwardError(*err,__LINE__,NULL);
  
  kernel = adaptMH_init(ndim, sig2,nadapt,adaptifaccepted,c_enlarge,k_damp,err);
  forwardError(*err,__LINE__,NULL);
  
  rc_close(&rca);
  
  return kernel;
}

parabox * parabox_from_rc(confFile *rc, char* root,error **err) {
  confFile *rca;
  parabox *pb;
  int ndim,ndum,i;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  ndim = rc_array_size(rca,"min",err);
  forwardError(*err,__LINE__,NULL);
  ndum = rc_array_size(rca,"max",err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(ndim!=ndum,pb_base,
                 "min and max arrays have different sizes (%d and %d)",
                 *err, __LINE__,NULL,ndim,ndum);
  
  pb = init_parabox(ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<ndim;i++) {
    double rmin,rmax;
    char par[100];
    sprintf(par,"min[%d]",i+1);
    rmin = rc_get_real(rca,par,err);
    forwardError(*err,__LINE__,NULL);
    sprintf(par,"max[%d]",i+1);
    rmax = rc_get_real(rca,par,err);
    forwardError(*err,__LINE__,NULL);
    add_slab(pb,i,rmin,rmax,err);
    forwardError(*err,__LINE__,NULL);
  }
  rc_close(&rca);
  return pb;
}

mvdens* rc_mvdens(confFile *rc, mvdens* mv, error **err) {
  int ndim,ndum;
  double *mn;
  mvdens *mmv;
  
  // get the mean
  ndum = rc_get_real_array(rc,"mean",&mn,err);
  forwardError(*err,__LINE__,NULL);
  
  if (mv==NULL) {
    ndim = ndum;
    mmv = mvdens_alloc(ndim,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    mmv = mv;
    ndim = mmv->ndim;
  }
  testErrorRetVA(ndum!=ndim,pmc_dimension,
                 "'mean' in %s does not have the right number elements (got %d expected %d)",
                 *err,__LINE__,NULL,rc->alias_prefix,ndum,ndim);
  memcpy(mmv->mean,mn,sizeof(double)*ndim);

  // get the var
  rc_get_real_assquarematrix(rc,"var",&mn,ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  memcpy(mmv->std,mn,sizeof(double)*ndim*ndim);
  
  // get the other parameters
  mmv->chol = rc_safeget_integer(rc,"chol",0,err);
  forwardError(*err,__LINE__,NULL);
  if (mmv->chol) {
    mmv->detL = determinant(mmv->std, ndim);
  }
  
  mmv->band_limit = rc_safeget_integer(rc,"bdl",ndim,err);
  forwardError(*err,__LINE__,NULL);

  mmv->df = rc_safeget_integer(rc,"df",-1,err);
  forwardError(*err,__LINE__,NULL);      
  
  return mmv;
}

double* rc_get_weight(confFile *rc,int sz,error **err) {
  double *w;
  double sw;
  int i;
  int hw;
  
  hw = rc_has_key(rc,"weight",err);
  forwardError(*err,__LINE__,NULL);
  
  if (hw==0) {
    w = rc_add_malloc(rc,sizeof(double)*sz,err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<sz;i++) {
      w[i] = 1./sz;
    }
    return w;
  }
  rc_get_real_asarray(rc,"weight",&w,sz,err);
  forwardError(*err,__LINE__,NULL);
  
  sw = 0;
  for(i=0;i<sz;i++) {
    testErrorRetVA(w[i]<0,pmc_negWeight,"weight %d is negative (%g)",*err,__LINE__,NULL,i,w[i]);
    sw += w[i];
  }
  testErrorRet(sw==0,pmc_negWeight,"all weights are 0",*err,__LINE__,NULL);
  for(i=0;i<sz;i++) {
    w[i] = w[i]/sw;
  }
  return w;
}

mix_mvdens* rc_mix_mvdens(confFile* rc, error **err) {
  int flg,ic,ndim,ncomp,ndum;
  mix_mvdens *mmv;
  confFile *rca;
  double *w;
  
  flg = rc_has_key(rc,"fromfile",err);
  forwardError(*err,__LINE__,NULL);
  
  if (flg==1) {
#ifdef _WITH_HDF5_
    char *hdfname;
    
    hdfname = rc_get_string(rc,"fromfile",err);
    forwardError(*err,__LINE__,NULL);
    
    mmv = mix_mvdens_hdfdwnp(hdfname,err);
    forwardError(*err,__LINE__,NULL);
    
    flg = rc_has_key(rc,"extract_component",err);
    forwardError(*err,__LINE__,NULL);
    if (flg==1) {
      long* lcomp;
      int ic;
      mix_mvdens *rmv;
      double *w;
      

      ncomp = rc_get_integer_array(rc,"extract_component",&lcomp,err);
      forwardError(*err,__LINE__,NULL);
      
      rmv = mix_mvdens_alloc(ncomp,mmv->ndim,err);
      forwardError(*err,__LINE__,NULL);
      
      for(ic=0;ic<ndum;ic++) {
        int oc;
        oc = lcomp[ic];
        testErrorRetVA(oc>=mmv->ncomp,pmc_dimension,"asked for too many component (got %d max %d)",*err,__LINE__,NULL,oc,mmv->ncomp);
        rmv->wght[ic] = mmv->wght[oc];
        memcpy(rmv->comp[ic]->mean,mmv->comp[oc]->mean,sizeof(double)*mmv->ndim);
        memcpy(rmv->comp[ic]->std,mmv->comp[oc]->std,sizeof(double)*mmv->ndim*mmv->ndim);
        rmv->comp[ic]->band_limit = rmv->comp[oc]->band_limit;
        rmv->comp[ic]->df = rmv->comp[oc]->df;
        rmv->comp[ic]->chol = rmv->comp[oc]->chol;
        rmv->comp[ic]->detL = rmv->comp[oc]->detL;
        mix_mvdens_free(&mmv);
        mmv = rmv;
      }
    }
#else
    testErrorRet(1==1,pmc_io,"Cannot do that without hdf5. Recompile with hdf5 !",*err,__LINE__,NULL);
#endif
  } else {
    ncomp = rc_get_integer(rc,"ncomp",err);
    forwardError(*err,__LINE__,NULL);

    ndim = rc_get_integer(rc,"ndim",err);
    forwardError(*err,__LINE__,NULL);

    mmv = mix_mvdens_alloc(ncomp,ndim,err);
    forwardError(*err,__LINE__,NULL);

    // deal with weight
    w = rc_get_weight(rc,ncomp,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(mmv->wght,w,sizeof(double)*ncomp);

    flg = rc_has_key(rc,"draw_from",err);
    forwardError(*err,__LINE__,NULL);

    if (flg==1) {
      // initialize by randomizeing under the first component
      gsl_rng* r;
      rca = rc_alias(rc,"draw_from",err);
      forwardError(*err,__LINE__,NULL);

      r = rc_gsl_rng(rca,"",err);
      forwardError(*err,__LINE__,NULL);

      rc_mvdens(rca,mmv->comp[0],err);
      forwardError(*err,__LINE__,NULL);
      for (ic=1;ic<ncomp;ic++) {
        mvdens_ran(mmv->comp[ic]->mean, mmv->comp[0], r,err);
        forwardError(*err,__LINE__,NULL);
        mmv->comp[ic]->chol = mmv->comp[0]->chol;
        mmv->comp[ic]->band_limit = mmv->comp[0]->band_limit;
        mmv->comp[ic]->df = mmv->comp[0]->df;
        mmv->comp[ic]->detL = mmv->comp[0]->detL;
        memcpy(mmv->comp[ic]->std,mmv->comp[0]->std,sizeof(double)*ndim*ndim);
      }
      rc_close(&rca);
      gsl_rng_free(r);
    } else {
      // I read everything from lua

      // deal with components
      ndum = rc_array_size(rc,"component",err);
      forwardError(*err,__LINE__,NULL);
      testErrorRetVA(ndum!=ncomp,pmc_dimension,"not the right number of elements in component (got %d expected %d)",*err,__LINE__,NULL,ndum,ncomp);

      for(ic=0;ic<ncomp;ic++) {
        char scomp[100];

        sprintf(scomp,"component[%d]",ic+1);
        rca = rc_alias(rc,scomp,err);
        forwardError(*err,__LINE__,NULL);

        rc_mvdens(rca,mmv->comp[ic],err);
        forwardError(*err,__LINE__,NULL);

        rc_close(&rca);
      }
    }    
  }
  flg = rc_has_key(rc,"save",err);
  forwardError(*err,__LINE__,NULL);
  if (flg==1) {
    char *file;
    FILE *f;
    file = rc_get_string(rc,"save",err);
    forwardError(*err,__LINE__,NULL);    
#ifdef _WITH_HDF5_
    mix_mvdens_hdfdump(mmv,file,err);
    forwardError(*err,__LINE__,NULL);    
#else
    f = fopen_err(file,"w",err);
    forwardError(*err,__LINE__,NULL);    
    mix_mvdens_dump(f,mmv);
    fclose(f);
#endif    
  }
  return mmv;
}

distribution* rcinit_mix_mvdens(confFile* rc,char* root, error **err) {
  mix_mvdens *mmv;
  distribution *trg;
  confFile *rca;  
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  mmv = rc_mix_mvdens(rca,err);
  forwardError(*err,__LINE__,NULL);
  
  trg = mix_mvdens_distribution(mmv->ndim, (void*) mmv, err);
  forwardError(*err,__LINE__,NULL);
  
  rc_close(&rca);
  
  return trg;
}

pmc_simu* init_importance_from_rc(confFile *rc,char *root,error **err) {
  pmc_simu *psim;
  long nsamples;
  confFile *rca;
  parabox *pb;
  int prop_print_step;
  distribution *proposal, *target;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);

  target = init_distribution_from_rc(rca,"target",err);
  forwardError(*err,__LINE__,NULL);
  
  proposal = init_distribution_from_rc(rca,"proposal",err);
  forwardError(*err,__LINE__,NULL);
  
  testErrorRetVA(proposal->ndim!=target->ndim,pmc_dimension,"incompatible dim for proposal (%d) and target (%d)",*err,__LINE__,NULL,
                 proposal->ndim,target->ndim);
  
  
  nsamples = rc_get_integer(rca,"nsample",err);
  forwardError(*err,__LINE__,NULL);
  
  
  psim = pmc_simu_init_plus_ded(nsamples,target->ndim,target->n_ded,err);
  forwardError(*err,__LINE__,NULL);

  pb = parabox_from_rc(rca,"pb",err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(pb->ndim!=target->ndim,pmc_dimension,"incompatible dim for pb (got %d expected %d)",*err,__LINE__,NULL,
                 pb->ndim,target->ndim);
  
  pmc_simu_init_target(psim,target,pb,err);
  forwardError(*err,__LINE__,NULL);
  
  prop_print_step = rc_safeget_integer(rca,"print_pc",0,err);
  forwardError(*err,__LINE__,NULL);
  
  pmc_simu_init_proposal(psim, proposal, prop_print_step, err);
  forwardError(*err,__LINE__,NULL);
  
  rc_close(&rca);
  
  return psim;
}

mcmc_run *init_mcmc_from_rc(confFile *rc,char *root,error **err) {
  mcmc_run *mcr;
  confFile *rca;
  distribution *target;
  mc_kernel *kernel;
  parabox *pb;
  int nsamples,nbatch,nd;
  double *pars0;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  target = init_distribution_from_rc(rca,"target",err);
  forwardError(*err,__LINE__,NULL);
  
  kernel = init_mc_kernel_from_rc(rca,"kernel",err);
  forwardError(*err,__LINE__,NULL);
  
  testErrorRetVA(kernel->ndim!=target->ndim,pmc_dimension,"incompatible dim for kernel (%d) and target (%d)",*err,__LINE__,NULL,
                 kernel->ndim,target->ndim);
  
  pb = parabox_from_rc(rca,"pb",err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(pb->ndim!=target->ndim,pmc_dimension,"incompatible dim for pb (got %d expected %d)",*err,__LINE__,NULL,
                pb->ndim,target->ndim);

  nsamples = rc_get_integer(rca,"nsample",err);
  forwardError(*err,__LINE__,NULL);

  nbatch = rc_safeget_integer(rca,"nbatch",nsamples,err);
  forwardError(*err,__LINE__,NULL);

  nd = rc_get_real_array(rca,"pars0",&pars0,err);
  forwardError(*err,__LINE__,NULL);
  
  testErrorRetVA(kernel->ndim!=nd,pmc_dimension,"incompatible dim for target (%d) and pars0 (%d)",*err,__LINE__,NULL,
                 kernel->ndim,nd);
  
  mcr = init_mcmc_run(nsamples,nbatch,kernel,target,pars0,pb,err);
  forwardError(*err,__LINE__,NULL);
  
  rc_close(&rca);
  return mcr;
}

pmc_simu* init_pmc_from_rc(confFile *rc,char *root,error **err) {
  pmc_simu *psim;
  
  psim = init_importance_from_rc(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  pmc_simu_init_pmc(psim, NULL, NULL, update_prop_rb_void, err);
  forwardError(*err,__LINE__,NULL);
  
  return psim;  
}

gsl_rng * rc_gsl_rng(confFile *rc, char *root,error **err) {
  gsl_rng *r;
  int flg;
  gsl_rng_type *rtype;
  confFile *rca;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);

  gsl_set_error_handler_off();
  
  rtype = gsl_rng_default;
  
  flg = rc_has_key(rca,"rng",err);
  forwardError(*err,__LINE__,NULL);

  if (flg==1) {
    const gsl_rng_type **rng_list;
    char* rngtype;
    int ii;
    rngtype = rc_get_string(rca,"rng",err);
    forwardError(*err,__LINE__,NULL);
    rng_list = gsl_rng_types_setup();
    for(ii=0;rng_list[ii]!=NULL;ii++) {
      if (strcasecmp((rng_list[ii])->name,rngtype)==0) {
        break;
      }
    }
    rtype = rng_list[ii];
    testErrorRetVA(rtype==NULL,pmc_undef,"Unknown generator '%s'",*err,__LINE__,NULL,rngtype);
  }
  
  r = gsl_rng_alloc(rtype);
  
  flg = rc_has_key(rca,"seed",err);
  forwardError(*err,__LINE__,NULL);
  if(flg==1) {
    long seed;
    seed = rc_get_integer(rca,"seed",err);
    forwardError(*err,__LINE__,NULL);
    gsl_rng_set(r,seed);
  }
  
  rc_close(&rca);
  return r;
}

char* get_itername(confFile* rc,char *key, int iter,error **err) {
  int flg;
  char *res,*oes;
  int ll;
  
  flg = rc_is_array(rc,key,err);
  forwardError(*err,__LINE__,NULL);
  
  if (flg) {
    char keyix[1000];
    sprintf(keyix,"%s[%d]",key,iter);
    res = rc_get_string(rc,keyix,err);
    forwardError(*err,__LINE__,NULL);
    return res;
  }
  res = rc_get_string(rc,key,err);
  forwardError(*err,__LINE__,NULL);
  
  ll = strlen(res);
  oes = rc_add_malloc(rc, ll+100,err);
  forwardError(*err,__LINE__,NULL);
  
  sprintf(oes,res,iter);
  return oes;
}

#ifdef _WITH_HDF5_
#include "hdf5.h"

void generic_hdfinit(char* fname,char *grpname, int ndim, int nded, int nsample, double *pars, double *deds, double *loglkl,hid_t *pfile_id, hid_t *pgroup_id,error **err) {
  hid_t       file_id, group_id;   
  herr_t      status;
  hsize_t     dims[2];
  
  /* Create a new file using default properties. */
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  testErrorRetVA(file_id<0,hdf5_base,"cannot create hdf5 file %s (got %d)",*err,__LINE__,,fname,file_id);
  
  group_id = H5Gcreate( file_id, grpname, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  testErrorRetVA(group_id<0,hdf5_base,"cannot create group /%s in file %s (got %d)",*err,__LINE__,,grpname,fname,group_id);

  status = H5LTset_attribute_int( group_id, ".", "ndim", &(ndim), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /%s/ndim in file %s (got %d)",*err,__LINE__,,grpname,fname,status);

  status = H5LTset_attribute_int( group_id, ".", "nded", &(nded), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /%s/nded in file %s (got %d)",*err,__LINE__,,grpname,fname,status);

  status = H5LTset_attribute_int( group_id, ".", "nsample", &(nsample), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /%s/nsample in file %s (got %d)",*err,__LINE__,,grpname,fname,status);

  dims[0] = nsample;
  dims[1] = ndim;
  status = H5LTmake_dataset(group_id,"pars",2,dims,H5T_NATIVE_DOUBLE,pars);
  testErrorRetVA(status<0,hdf5_base,"cannot save /%s/par in file %s (got %d)",*err,__LINE__,,grpname,fname,status);
  if (nded>0) {
    dims[1] = nded;
    status = H5LTmake_dataset(group_id,"deds",2,dims,H5T_NATIVE_DOUBLE,deds);
    testErrorRetVA(status<0,hdf5_base,"cannot save /%s/ded in file %s (got %d)",*err,__LINE__,,grpname,fname,status);
  }

  status = H5LTmake_dataset(group_id,"log_lkl",1,dims,H5T_NATIVE_DOUBLE,loglkl);
  testErrorRetVA(status<0,hdf5_base,"cannot save /%s/log_lkl in file %s (got %d)",*err,__LINE__,,grpname,fname,status);

  *pfile_id = file_id;
  *pgroup_id = group_id;
  
  return;
}
  
  
void mcmc_simu_hdfdump(mcmc_run* mcr,char* fname,error **err) {
  hid_t       file_id, group_id;   
  herr_t      status;

  generic_hdfinit(fname,"mcmc",mcr->ndim,mcr->nded,mcr->current,mcr->pars,mcr->deds,mcr->loglkl,&file_id,&group_id,err);
  forwardError(*err,__LINE__,);

  status = H5LTset_attribute_int( group_id, ".", "naccepted", &(mcr->naccepted), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mcmc/naccepted in file %s (got %d)",*err,__LINE__,,fname,status);
  //close
  status = H5Gclose(group_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save /mcmc in file %s (got %d)",*err,__LINE__,,fname,status);

  status = H5Fclose(file_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save file %s (got %d)",*err,__LINE__,,fname,status);

  return;
}
void pmc_simu_hdfdump(pmc_simu* psim,char* fname,error **err) {
  hid_t       file_id, group_id;   
  herr_t      status;
  hsize_t     dims[2];
  
  normalize_importance_weight(psim,err);
  forwardError(*err,__LINE__,);
  
  generic_hdfinit(fname,"pmc",psim->ndim,psim->n_ded,psim->nsamples,psim->X,psim->X_ded,psim->log_rho,&file_id,&group_id,err);
  forwardError(*err,__LINE__,);
  
  status = H5LTset_attribute_double( group_id, ".", "logSum", &(psim->logSum), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /pmc/logSum in file %s (got %d)",*err,__LINE__,,fname,status);
  
  dims[0] = psim->nsamples;
  status = H5LTmake_dataset(group_id,"flag",1,dims,H5T_NATIVE_SHORT,psim->flg);
  testErrorRetVA(status<0,hdf5_base,"cannot save /pmc/flag in file %s (got %d)",*err,__LINE__,,fname,status);
  status = H5LTmake_dataset(group_id,"weight",1,dims,H5T_NATIVE_DOUBLE,psim->weights);
  testErrorRetVA(status<0,hdf5_base,"cannot save /pmc/weight in file %s (got %d)",*err,__LINE__,,fname,status);

  //close
  status = H5Gclose(group_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save /pmc in file %s (got %d)",*err,__LINE__,,fname,status);
  
  status = H5Fclose(file_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save file %s (got %d)",*err,__LINE__,,fname,status);
  
  return;
}
#endif

