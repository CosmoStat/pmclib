/*
 *  mcmc.c
 *  likely
 *
 *  Created by Karim Benabed on 12/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef __PLANCK__
#include "HL2_likely/pmclib/mcmc.h"
#else
#include "mcmc.h"
#endif

mc_kernel *mc_kernel_init(size_t ndim, 
                          transition_func *sample,transition_ratio_func *log_pdf, 
                          void* data, posterior_log_free* freef, int nadapt, int adaptifaccepted, 
                          adapt_func *adapt, error ** err) {
                            
  mc_kernel* self;
  self = malloc_err(sizeof(mc_kernel),err);
  forwardError(*err,__LINE__,NULL);
  self->ndim = ndim;
  self->extra = data;
  self->sample = sample;
  self->transition_ratio = log_pdf;
  self->nadapt = nadapt;
  self->adaptifaccepted = adaptifaccepted;
  self->adapt = adapt;
  self->free = freef;
  self->dlhandle = NULL;
  return self;        
}

void free_mc_kernel(void** pmcl) {
  mc_kernel *mcl;
  mcl = *pmcl;
  mc_law_free((mc_law**)pmcl,mcl->free);
  if (mcl->dlhandle!=NULL) {
    dlclose(mcl->dlhandle);
  }
  
}

mcmc_run *init_mcmc_run(int nsample, int nbatch, mc_law* kernel, distribution *target, double* pars0, parabox * pb, error **err) {
  mcmc_run *mcr;
    
  testErrorRetVA(kernel->ndim!=target->ndim,mcmc_dimension,
                 "transition kernel and target have different number of dims (got %d and %d)",
                 *err,__LINE__,NULL,kernel->ndim, target->ndim);
   testErrorRetVA(kernel->ndim!=pb->ndim,mcmc_dimension,
                  "transition kernel and parabox have different number of dims (got %d and %d)",
                  *err,__LINE__,NULL,kernel->ndim, pb->ndim);

  mcr = malloc_err(sizeof(mcmc_run),err);
  forwardError(*err,__LINE__,0);
  
  mcr->ndim = target->ndim;
  mcr->nded = target->n_ded;
  mcr->ntot = target->ndim + target->n_ded;
  
  mcr->nsample = nsample;
  mcr->nbatch = nbatch;
  mcr->lastbatch = 0;
  
  mcr->nadapt = kernel->nadapt;
  mcr->adaptifaccepted = kernel->adaptifaccepted;
  mcr->lastadapt = 1;
  
  mcr->kernel = kernel;
  mcr->target = target;
  
  mcr->current = 1;
  mcr->naccepted = 0;
  mcr->cntadapt = 0;
  
  mcr->pars = malloc_err(sizeof(double)*((mcr->ntot+1)*mcr->nsample),err);
  forwardError(*err,__LINE__,0);
  mcr->deds = mcr->pars + (mcr->ndim)*mcr->nsample;
  mcr->loglkl = mcr->deds + (mcr->nded)*mcr->nsample;
  
  memcpy(mcr->pars,pars0,sizeof(double)*mcr->ndim);
  
  mcr->loglkl[0] = mcr->target->log_pdf(mcr->target, pars0,err);
  forwardError(*err,__LINE__,0);
  
  if(mcr->target->retrieve!=NULL) {
    mcr->target->retrieve(mcr->target,mcr->pars,err);
    forwardError(*err,__LINE__,0);
  }

  mcr->pb = pb;
  
  return mcr;
}

int mcmc_run_batch(mcmc_run *mcr, gsl_rng *r, error **err) {
  int flg,mx,i;
  
  
  if (mcr->current == mcr->nsample) {
    return 0;
  }
  
  mx = mcr->current + mcr->nbatch;
  if (mx>=mcr->nsample) {
    mx = mcr->nsample;
  }
  mx -= mcr->current;
  
  for(i=0;i<mx;i++) {
    mcr->loglkl[mcr->current] = mcr->loglkl[mcr->current-1];
    flg = mcmc_step(
              &mcr->pars[(mcr->current-1)*(mcr->ndim)],&mcr->pars[mcr->current*(mcr->ndim)],
              &mcr->loglkl[mcr->current],&distribution_lkl,mcr->target,mcr->kernel,mcr->pb,
              r,err);
    forwardError(*err,__LINE__,0);
    /*_DEBUGHERE_("%g,%g -> %g,%g :: %d",mcr->pars[(mcr->current-1)*(mcr->ndim)],mcr->pars[(mcr->current-1)*(mcr->ndim)+1],
                  mcr->pars[(mcr->current)*(mcr->ndim)],mcr->pars[(mcr->current)*(mcr->ndim)+1],flg); */ 
    if (flg != MC_AF_ACCEPT) { //rejected
      // copy n-1 pars and log into n
      memcpy(&mcr->pars[(mcr->current)*(mcr->ndim)],&mcr->pars[(mcr->current-1)*(mcr->ndim)],sizeof(double)*mcr->ndim);
      if(mcr->nded>0) {
        memcpy(&mcr->deds[(mcr->current)*(mcr->nded)],&mcr->deds[(mcr->current-1)*(mcr->nded)],sizeof(double)*mcr->nded);        
      }
      mcr->loglkl[mcr->current] = mcr->loglkl[mcr->current-1];
      if (mcr->adaptifaccepted==0) {
        mcr->cntadapt++;
      }
    } else { //accept
      if(mcr->target->retrieve!=NULL) {
        mcr->target->retrieve(mcr->target,&mcr->deds[mcr->current*(mcr->nded)],err);
        forwardError(*err,__LINE__,0);
      }
      mcr->cntadapt++;
      mcr->naccepted++;
    }
    mcr->current++;
    
    // deal with adaptation if needed
    if (mcr->nadapt>0 && mcr->cntadapt==mcr->nadapt) {
      _DEBUGHERE_("adapting at %d points",mcr->current);
      mcr->kernel->adapt(mcr->kernel,mcr,err);
      forwardError(*err,__LINE__,0);
      mcr->lastadapt = mcr->current;
      mcr->cntadapt = 0;
    }
  } 
  
  return 1;
}

void free_mcmc_run(void** pmcr) {
  mcmc_run *mcr;
  
  mcr = *pmcr;
  
  free(mcr->pars);
  free_distribution(&(mcr->target));
  free_mc_kernel((void**)(&(mcr->kernel)));
  free_parabox(&(mcr->pb));
  free(mcr);
  *pmcr = NULL;
}
 
mc_kernel * adaptMH_init(int ndim, double *sig0, int nadapt, int adaptifaccepted, double c_enlarge, double k_damp,error **err) {
  adaptMH *payload;
  mc_kernel *kernel;
  
  payload = malloc_err(sizeof(adaptMH),err);
  forwardError(*err,__LINE__,NULL);
  
  payload->gauss = mvdens_alloc(ndim,err);
  forwardError(*err,__LINE__,NULL);
  mvdens_from_meanvar(payload->gauss, NULL, sig0, 1);
  
  payload->c_enlarge = c_enlarge;
  payload->k_damp = k_damp;

  payload->buf = malloc_err(sizeof(double)*(ndim+2)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  payload->mu_t = payload->buf + ndim;
  payload->sig_t = payload->mu_t + ndim;
  memset(payload->mu_t,0,sizeof(double)*ndim);
  memset(payload->sig_t,0,sizeof(double)*ndim*ndim);
  
  kernel = mc_kernel_init(ndim, adaptMH_sample,NULL, payload, adaptMH_free, nadapt, adaptifaccepted, adaptMH_adapt, err);
  forwardError(*err,__LINE__,NULL);
  
  return kernel;
}

void adaptMH_free(void** pmh) {
  adaptMH *mh;
  
  mh = *pmh;
  mvdens_free(&mh->gauss);
  free(mh->buf);
  free(mh);
  *pmh = NULL;
}
 

void adaptMH_sample(double* state, void* self, gsl_rng* rng, error** err) {
  double* res;
  mc_kernel *kernel;
  adaptMH* mh;
  
  kernel = self;
  mh = kernel->extra;
  
  size_t i;
  res = mvdens_ran(mh->buf,mh->gauss,rng,err);
  forwardError(*err,__LINE__,);
  for(i=0;i<((mc_law*)self)->ndim;i++) {
    state[i]+=res[i];
  }
  return;
}


void adaptMH_adapt(mc_kernel *kernel, mcmc_run* mc, error **err) {
  int i,ip,jp,p;
  double psnk;
  double *mean,*var;
  adaptMH *mh;
  
  mh = kernel->extra;
  
  mean = malloc_err(sizeof(double)*mc->ndim*mc->ndim,err);
  forwardError(*err,__LINE__,);

  
  p = mc->current - mc->lastadapt;
  if (mh->k_damp!=0) {
    psnk = pow(p*1./mc->current,mh->k_damp);    
  } else {
    psnk = 1;
  }
  
  // compute mean

  memset(mean,0,sizeof(double)*mc->ndim);
  
  for(i=mc->lastadapt;i<mc->current;i++) {
    for(ip=0;ip<mc->ndim;ip++) {
      mean[ip] += mc->pars[i*mc->ndim+ip];
    }
  }
  
  for(ip=0;ip<mc->ndim;ip++) {
    mh->mu_t[ip] = (1-psnk) * mh->mu_t[ip] + psnk/p * mean[ip];
  }
  
  // compute var
  var = mean;
  memset(var,0,sizeof(double)*mc->ndim*mc->ndim);
  
  for(i=mc->lastadapt;i<mc->current;i++) {
    for(ip=0;ip<mc->ndim;ip++) {
      double vip;
      vip = mc->pars[i*mc->ndim+ip]-mh->mu_t[ip];
      var[ip*mc->ndim+ip] +=   vip * vip;
      for(jp=ip+1;jp<mc->ndim;jp++) {
        var[ip*mc->ndim+jp] +=  vip * (mc->pars[i*mc->ndim+jp] - mh->mu_t[jp]);
      }
    }
  }
  
  for(ip=0;ip<mc->ndim;ip++) {
    mh->sig_t[ip*mc->ndim+ip] = (1-psnk) * mh->sig_t[ip*mc->ndim+ip] + psnk/p * var[ip*mc->ndim+ip];
    for(jp=ip+1;jp<mc->ndim;jp++) {
      mh->sig_t[ip*mc->ndim+jp] = (1-psnk) * mh->sig_t[ip*mc->ndim+jp] + psnk/p * var[ip*mc->ndim+jp];
      mh->sig_t[jp*mc->ndim+ip] =mh->sig_t[ip*mc->ndim+jp];
    }
  }
  
  free(mean);
  
  mvdens_from_meanvar(mh->gauss, NULL, mh->sig_t, mh->c_enlarge*mh->c_enlarge/mc->ndim);
  
}
 
int mcmc_step(double *pstate,double* nstate, double *log_pdf_value, posterior_log_pdf_func *pdf_func, void *lkl_extra, mc_law* law, parabox *pb, gsl_rng* rng,error **err) {
  double current_lkl, previous_lkl,transition,proba,alea;

  
  previous_lkl = *log_pdf_value;
  
  /* draw new sample */
  memcpy(nstate,pstate,law->ndim*sizeof(double));
  law->sample(nstate,law,rng,err);
  forwardError(*err,__LINE__,0);
  
  /* test box */
  if (isinBox(pb,nstate,err)==0) {
    forwardError(*err,__LINE__,0);
    *log_pdf_value=-MC_INF;
    return MC_AF_OUTOFBOX;
  }
  
  /* compute transition ratio */
  transition=1;
  if (law->transition_ratio!=NULL) {
    transition = law->transition_ratio(pstate,nstate,law,err);
    forwardError(*err,__LINE__,0);      
  }
  
  /* compute lkl */
  current_lkl=pdf_func(lkl_extra,nstate,err);
  forwardError(*err,__LINE__,0);
  *log_pdf_value = current_lkl;
  /* compute proba */
  proba = exp(current_lkl-previous_lkl)*transition;
  //fprintf(stderr,"%g %g -> %g\n",previous_lkl,current_lkl,proba);
  /* accept/reject */
  //fprintf(stderr,"%g\n",proba);
  if (proba >=1) { /* inconditionnaly accept ! */
    //fprintf(stderr,"auto\n");
    //return 2; 
    return MC_AF_ACCEPT;
  }
  alea = gsl_ran_flat(rng,0.0,1.0);
  if (alea < proba) { /* accept */
    return MC_AF_ACCEPT;
  }
  /* reject */
  return MC_AF_REJECT;
}

mc_law * mc_law_init(size_t ndim, transition_func *sample,transition_ratio_func *log_pdf,void* extra,error ** err) {
  mc_law* self;
  self = mc_kernel_init(ndim,sample,log_pdf, 
                        extra,NULL, 0, 0, 
                        NULL, err);
  forwardError(*err,__LINE__,NULL);
  return self;  
}

void mc_law_free(mc_law** self, void (*extra_free)(void**)) {
  if ((*self)->extra !=NULL) {
    if (extra_free!=NULL) {
      extra_free(&((*self)->extra));
    }
    free((*self)->extra);
  }
  free(*self);
  self=NULL;
}

mc_law * mc_mvdens_init(mvdens* mself,error **err) {
  mc_law *self;
  self=mc_law_init(mself->ndim,&mc_mvdens_sample,NULL,mself,err);
  forwardError(*err,__LINE__,NULL);
  return self;
}

void mc_mvdens_free(mc_law** self) {
  mc_law_free(self, (void(*)(void**))(&mvdens_free));
}

void mc_mvdens_sample(double* state, void* self, gsl_rng* rng, error** err) {
  double* res;
  size_t i;
  res = mvdens_ran(NULL,((mc_law*)self)->extra,rng,err);
  forwardError(*err,__LINE__,);
  for(i=0;i<((mc_law*)self)->ndim;i++) {
    state[i]+=res[i];
  }
  free(res);
  return;
}

double trivial_transition_func(double* pstate, double* nstate, void* self, error** err)
{
   return 1.0;
}


/* ============================================================ *
 * Obsolete! Use estimate_param_covar_weight instead with       *
 * weigh=NULL.							*
 * ============================================================ */
void estimate_param_covar(size_t ndim, size_t nsamples, size_t psamples, const double* params, 
			  double *pmean, double *pvar, error **err)
{
   estimate_param_covar_weight(ndim, nsamples, psamples, params, NULL, pmean, pvar, err);
   forwardError(*err, __LINE__,);
}

/* ================================================================  *
 * Weighted mean and covariance. For unweighted estimates, call      *
 * with weight=NULL.						     *
 * See http://en.wikipedia.org/wiki/Sample_variance#Weighted_samples *
 * ================================================================= */
void estimate_param_covar_weight(size_t ndim, size_t nsamples, size_t psamples, const double* params, 
			  const double *weight, double *pmean, double *pvar, error **err)
{
  size_t i, j, isample, insample, indim;
  double *ppmean, w, sumw=0.0, sumww=0.0;
  
  ppmean = malloc_err(ndim*sizeof(double), err);   forwardError(*err, __LINE__,);
  
  /* Initialisation */
  for(i=0;i<ndim;i++) {
    indim=i*ndim;
    for(j=0;j<ndim;j++) {
      pvar[indim+j] *= psamples-1;
    }
    ppmean[i]  = pmean[i];
    pmean[i]  *= psamples;
  }

  /* Sum of weights */
  if (weight!=NULL) {
    for(isample=0,sumw=sumww=0.0; isample<nsamples; isample++) {
      sumw  += weight[isample];
      sumww += weight[isample]*weight[isample];
    }
    //fprintf(stderr, "Sum of weights = %g\n", sumw);
  }

  /* Mean */
  for(isample=0; isample<nsamples-psamples; isample++) {
    insample = isample*ndim;
    for(i=0;i<ndim;i++) {
      //fprintf(stderr,"%d %d %g\n",isample,i,params[insample+i]);
      if (weight==NULL) w = 1.0;
      else w = weight[isample];
      pmean[i] += w*params[insample+i];
    }
  }

  for(i=0;i<ndim;i++) {
    if (weight==NULL) pmean[i] /= nsamples;
    else pmean[i] /= sumw;
  }   
  
  /* Covariance */
  for(isample=0;isample<nsamples-psamples;isample++) {
    insample=isample*(ndim);
    for(i=0;i<ndim;i++) {
      indim=i*ndim;
      for(j=0;j<ndim;j++) {
	if (weight==NULL) w = 1;
	else w = weight[isample];
        pvar[indim+j] += w*(params[insample+i]-pmean[i])*(params[insample+j]-pmean[j]); 
      }
    }
  }
  for(i=0;i<ndim;i++) {
    indim=i*ndim;
    for(j=0;j<ndim;j++) {
      pvar[indim+j] += psamples*(ppmean[i]-pmean[i])*(ppmean[j]-pmean[j]);
      if (weight==NULL) pvar[indim+j] /= (nsamples-1);
      //if (i==0 && j==0) fprintf(stderr, "*** %g %g\n", pvar[indim+j], 1.0-sumww);
      else pvar[indim+j] /= 1.0-sumww; // ????
    }
  }   
}

