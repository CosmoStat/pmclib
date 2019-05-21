/*
 *  pmc.c
 *  likely
 *
 *  Created by Karim Benabed on 10/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef __PLANCK__
#include "HL2_likely/pmclib/pmc.h"
#else
#include "pmc.h"
#endif
#include <math.h>

/******************************************************************************/
/******************** high level interface ************************************/
/******************************************************************************/

pmc_simu* pmc_simu_init(long nsamples, int ndim, error **err)
{
   pmc_simu *self;
   self = pmc_simu_init_plus_ded(nsamples, ndim, 0, err);
   forwardError(*err, __LINE__, 0);
   return self;
}

pmc_simu* pmc_simu_init_plus_ded(long nsamples, int ndim, int n_ded, error **err) {
  pmc_simu* self;

  self = malloc_err(sizeof(pmc_simu), err);
  forwardError(*err, __LINE__, NULL);

  self->nsamples = nsamples;
  self->ndim     = ndim;
  self->n_ded    = n_ded;
  self->logSum   = 0;
  self->isLog    = 0;
  self->retry=0;
  
  self->data = malloc_err(nsamples*( 
                                (ndim+2)*sizeof(double)
				+n_ded*sizeof(double)
                                +sizeof(short) 
                                +sizeof(size_t)),err);
  forwardError(*err,__LINE__,NULL);


  
  self->X       = self->data;
  self->X_ded   = (void*)(((char*) self->X) + nsamples * ndim * sizeof(double));
  self->weights = (void*)(((char*) self->X_ded) + nsamples * n_ded * sizeof(double));
  self->log_rho = (void*)(((char*) self->weights) + nsamples * sizeof(double));
  self->indices = (void*)(((char*) self->log_rho) + nsamples * sizeof(double));
  self->flg     = (void*)(((char*) self->indices) + nsamples * sizeof(size_t));
  //fprintf(stderr,"pos : %p %p %d %d",self->X, self->weights, self->X-self->weights,nsamples * self->ndim * sizeof(double));

  self->prop_print_step = _PRSTP_;
  
  self->pb = NULL;

  self->proposal = NULL;
  self->target = NULL;
  
  self->filter_data = NULL;
  self->pmc_filter = NULL;
  self->pmc_update = NULL;

  self->mpi_rank = 0;
  self->mpi_size = 0;
  
  return self;
}

void pmc_simu_free(pmc_simu ** melf) {
  pmc_simu *self;
  self=*melf;
  free(self->data);
  if (self->target!=NULL)
    free_distribution(&self->target);
  if (self->proposal!=NULL)
    free_distribution(&self->proposal);
  
  free(self);
  *melf=NULL;
}

void pmc_simu_init_target(pmc_simu* psim, 
                          distribution *target, 
                          parabox *pb, error **err) {

  if (psim->target != NULL) {
    free_distribution(&psim->target);
  }
  testErrorRetVA(target->ndim!=psim->ndim,pmc_incompat,
               "Target invalid, got %d dims, while pmc initialized for %d",*err,__LINE__,,target->ndim,psim->ndim);
  testErrorRetVA(target->n_ded!=psim->n_ded,pmc_incompat,
                 "Target invalid, got %d deduced parameters, while pmc initialized for %d",
		 *err,__LINE__,,target->n_ded,psim->n_ded);
  testErrorRet(target->n_ded!=0 && target->retrieve==NULL,pmc_incompat,
               "Target invalid, expect deduced parameters, but nu function...",*err,__LINE__,);

  testErrorRet(target->log_pdf==NULL,pmc_incompat,
               "Target invalid, absent log_pdf function",*err,__LINE__,);
  
  psim->target = target;
  
  if (psim->pb!=NULL) {
    free_parabox(&psim->pb);
  }
  psim->pb = pb;
}

void pmc_simu_init_proposal(pmc_simu* psim, 
                            distribution *proposal, 
                            int prop_print_step, error **err) {
  
  if (psim->proposal != NULL) {
    free_distribution(&psim->proposal);
  }
  testErrorRetVA(proposal->ndim!=psim->ndim,pmc_incompat,
                 "Proposal invalid, got %d dims, while pmc initialized for %d",
                 *err,__LINE__,,proposal->ndim,psim->ndim);
  
  testErrorRet(proposal->log_pdf==NULL,pmc_incompat,
               "Proposal invalid, absent log_pdf function",*err,__LINE__,);

  testErrorRet(proposal->simulate==NULL,pmc_incompat,
               "Proposal invalid, absent simulate function",*err,__LINE__,);

  psim->proposal = proposal;
    
  psim->prop_print_step = prop_print_step;
}

void pmc_simu_init_pmc(pmc_simu* psim, 
                              void* filter_data,
                              filter_func* pmc_filter,
                              update_func* pmc_update) {
  
  psim->filter_data = filter_data;
  psim->pmc_filter  = pmc_filter;
  psim->pmc_update  = pmc_update;
}

void pmc_simu_init_classic_importance(pmc_simu* psim,
                                      distribution *target, parabox *pb,
                                      void *prop_data,int prop_print_step,error **err) {
  distribution *prop;
  // init target
  pmc_simu_init_target(psim,target,pb,err);
  forwardError(*err,__LINE__,);
  
  // init proposal
  prop = mix_mvdens_distribution(target->ndim, prop_data, err);
  forwardError(*err,__LINE__,);
  
  pmc_simu_init_proposal(psim, prop, prop_print_step, err);
  forwardError(*err,__LINE__,);
}

void pmc_simu_init_classic_pmc(pmc_simu* psim,
                               distribution *target, parabox *pb,
                               void *prop_data, int prop_print_step, void* filter_data,
                               filter_func* pmc_filter, error **err) {
  
  pmc_simu_init_classic_importance(psim,target,pb,
                                   prop_data,prop_print_step,err);
  forwardError(*err,__LINE__,);

  pmc_simu_init_pmc(psim, filter_data, pmc_filter, update_prop_rb_void);
}

size_t pmc_simu_importance(pmc_simu *psim, gsl_rng *r,error **err) {
  size_t nok;

  testErrorRet(psim->proposal->data==NULL,pmc_allocate,
               "proposal undefined",*err,__LINE__,0);
  testErrorRet(psim->target->data==NULL,pmc_allocate,
               "target undefined",*err,__LINE__,0);
  
  psim->proposal->simulate(psim,psim->proposal->data,r,psim->pb,err);
  forwardError(*err,__LINE__,0);
  nok = generic_get_importance_weight_and_deduced(psim,
                                                  psim->proposal,
                                                  &distribution_lkl,
                                                  &distribution_lkl,
                                                  &distribution_retrieve,
                                                  psim->target,err);
  forwardError(*err,__LINE__,0);
  return nok;
  
}

double pmc_simu_pmc_step(pmc_simu *psim, gsl_rng *r, error **err) {
  double perp;
  
  testErrorRet(psim->pmc_update==NULL,pmc_allocate,
               "update function undefined",*err,__LINE__,0);
  pmc_simu_importance(psim,r,err);
  forwardError(*err,__LINE__,0);
  
  if (psim->pmc_filter!=NULL) {
    psim->pmc_filter(psim,psim->filter_data,err);
    forwardError(*err,__LINE__,0);    
  }

  normalize_importance_weight(psim,err);
  forwardError(*err,__LINE__,0);
  
  
  psim->pmc_update(psim->proposal->data,psim,err);
  forwardError(*err,__LINE__,0);
  
  perp = perplexity(psim,MC_UNORM,err);
  forwardError(*err,__LINE__,0);
  return perp;  
}

/* Reallocates a PMC simulation to allow for a different sample size */ 
void pmc_simu_realloc(pmc_simu *psim,long newsamples,error **err) {
  if (psim->nsamples==newsamples) {
    return;
  }
  
  psim->nsamples = newsamples;
  free(psim->data);
  psim->data = malloc_err(newsamples*( 
                                      (psim->ndim+2)*sizeof(double) 
                                      +psim->n_ded*sizeof(double)
                                      +sizeof(short) 
                                      +sizeof(size_t)),err);
  forwardError(*err,__LINE__,);
  
  psim->X       = psim->data;
  psim->X_ded   = (void*)(((char*) psim->X) + newsamples * psim->ndim * sizeof(double));
  psim->weights = (void*)(((char*) psim->X_ded) + newsamples * psim->n_ded * sizeof(double));
  psim->log_rho = (void*)(((char*) psim->weights) + newsamples * sizeof(double));
  psim->indices = (void*)(((char*) psim->log_rho) + newsamples * sizeof(double));
  psim->flg     = (void*)(((char*) psim->indices) + newsamples * sizeof(size_t));
  psim->logSum  = 0;
  psim->isLog   = 0;
  
  return;
}

distribution* mix_mvdens_distribution(int ndim, void *mxmv, error **err) {
  distribution *dist;
  dist = init_distribution_full(ndim, mxmv, mix_mvdens_log_pdf_void, mix_mvdens_free_void,
				simulate_mix_mvdens_void, 0, NULL, err);
  forwardError(*err,__LINE__,NULL);
  return dist;
}

long simulate_mix_mvdens_void(struct _pmc_simu_struct_ *psim, void *proposal, gsl_rng * r,
			 parabox *pb, error **err)
{
   long i;
   i = simulate_mix_mvdens(psim, (mix_mvdens*)proposal, r, pb, err);
   forwardError(*err, __LINE__, 0);
   return i;
}

#define MNOK 5
long simulate_mix_mvdens(struct _pmc_simu_struct_ *psim, mix_mvdens *proposal, gsl_rng * r,
			 parabox *pb, error **err) {
  long i, inok;
  int ok;
  size_t ind;

  i=0;
  inok=0;
  while(i<psim->nsamples) {
    mix_mvdens_ran(psim->X+i*psim->ndim,&ind,proposal,r,err);
    forwardError(*err,__LINE__,0);
    ok=1;
    if (pb!=NULL) {

      if (isinBox(pb,psim->X+i*psim->ndim,err)==0) 
        ok=0;
      forwardError(*err,__LINE__,0);
    }

    if (ok) {	
      psim->indices[i] = ind;
      psim->flg[i] = 1;
      i++;
    } else {
      inok++;
      testErrorRetVA(inok>MNOK*psim->nsamples, pmc_tooManySteps,
                     "Too many points (%d) outside of box (inside: %d). Try to (1) decrease the variance of the initial\n"
                     "proposal (e.g. check the Fisher matrix); (2) increase nsamples",
                     *err, __LINE__, 0, inok, i);
    }
  }
  return i;
}
#undef MNOK


/******************************************************************************/
/******************** Simple IO ***********************************************/
/******************************************************************************/

void out_pmc_simu(const char *name, const pmc_simu *psim, error **err)
{
  FILE *OUT;
  int i, j;
  
  OUT = fopen_err(name, "w", err);
  forwardError(*err,__LINE__,);
  
  for (i=0; i<psim->nsamples; i++) {
    if (psim->flg[i]==0) continue;
    for (j=0; j<psim->ndim; j++) {
      fprintf(OUT, "% .5e ", psim->X[j+i*psim->ndim]);
    }
    for (j=0; j<psim->n_ded; j++) {
      fprintf(OUT, "% .5e ", psim->X_ded[j+i*psim->n_ded]);
    }
    fprintf(OUT, "% .5e\n", psim->weights[i]);
  }
  fclose(OUT);
}

void pmc_simu_print(FILE *where, pmc_simu *psim,double nrm) {
  long i;
  for(i=0;i<psim->nsamples;i++)
    print_step(where,psim->flg[i],nrm*(psim->weights[i]),psim->ndim,psim->X+i*psim->ndim);
}

#define packToFile(file,data,size,cursize,line) {\
size_t sz,tl; \
sz=size; \
testErrorRet(fwrite(data,sizeof(char),sz,file)!=size,pmc_io,"Cannot write to file",*err,line,); \
cursize += sz; \
tl=ftell(file); \
testErrorRetVA(tl!=cursize,pmc_io,"size and file size are not compatible (%ld %ld)",*err,line,,tl,cursize); \
}

void pmc_simu_dump(FILE *where,pmc_simu *psim,error **err) {
  size_t cursize;

  cursize = 0;
  
  // write header
  packToFile(where,"PMC_SIMD",sizeof(char)*8,cursize,__LINE__);
  packToFile(where,&(psim->nsamples),sizeof(long),cursize,__LINE__);
  packToFile(where,&(psim->ndim),sizeof(int),cursize,__LINE__);
  packToFile(where,&(psim->n_ded),sizeof(int),cursize,__LINE__);
  packToFile(where,&(psim->isLog),sizeof(int),cursize,__LINE__);
  packToFile(where,&(psim->logSum),sizeof(double),cursize,__LINE__);
 
  // write flag
  packToFile(where,(psim->flg),sizeof(short)*psim->nsamples,cursize,__LINE__);
  
  // write weights
  packToFile(where,(psim->weights),sizeof(double)*psim->nsamples,cursize,__LINE__);

  // write log_prop
  packToFile(where,(psim->log_rho),sizeof(double)*psim->nsamples,cursize,__LINE__);

  // write data
  packToFile(where,(psim->X),sizeof(double)*psim->nsamples*psim->ndim,cursize,__LINE__);

  // write deduced
  packToFile(where,(psim->X_ded),sizeof(double)*psim->nsamples*psim->n_ded,cursize,__LINE__);
  
  return;
}

pmc_simu *pmc_simu_from_file(FILE *PMCSIM, int nsamples, int npar, int n_ded,
			     mix_mvdens *proposal, int nclipw, error **err)
{
   pmc_simu *psim;
   fpos_t headpos;

   /* Here: No call to MPI init function, also used in non-MPI programs */
   psim = pmc_simu_init_plus_ded(nsamples, npar, n_ded, err);
   forwardError(*err, __LINE__, NULL);
   psim->X = read_mkmc_chain(npar, n_ded, nsamples, PMCSIM, &(psim->nsamples), &headpos,
			     psim->X_ded, psim->weights, psim->indices, err);
   forwardError(*err, __LINE__, NULL);

   fill_read_pmc_simu(psim, proposal, err);
   forwardError(*err, __LINE__, NULL);

   /* NEW: clip weights moved here */
   if (nclipw>0) {
      clip_weights(psim, nclipw, stderr, err);
      forwardError(*err, __LINE__, NULL);
   }

   return psim;
}


/* Necessary after a pmc simulation has been read from a file. */
double fill_read_pmc_simu(pmc_simu *psim, mix_mvdens *proposal, error **err)
{
  int i;
  double MR, MW, weight_sum, logwi;
  
  MR = MW = -1.0e30;
  weight_sum = 0.0;
  for (i=0; i<psim->nsamples; i++) {
    
    /*if (psim->weights[i]==0) {
     psim->flg[i] = 0;
     continue;
     }*/
    
    psim->flg[i] = 1;
    //weight_sum += psim->weights[i];
    
    if (proposal!=NULL) {
      psim->log_rho[i] = mix_mvdens_log_pdf(proposal, &(psim->X[i*psim->ndim]), err);
      forwardError(*err, __LINE__, 0.0);
      if (psim->log_rho[i]>MR) MR = psim->log_rho[i];
    }

    //logwi = log(psim->weights[i]);
    /* New (cosmo_pmc) v1.2: In read_mkmc_chain no exp of weights */
    logwi = psim->weights[i];

    if (logwi>MW) MW = logwi;
    /* weights have to be in log for normalize_importance_weight */
    psim->weights[i] = logwi;
  }

  psim->maxW   = MW;              /* As in get_importance_weight */
  psim->maxR   = MR;              /* As in get_importance_weight */
  psim->isLog  = 1;		  /* So that weights are normalized in next line */
  weight_sum = normalize_importance_weight(psim, err);
  forwardError(*err, __LINE__, 0.0);

  return weight_sum;
}


/******************************************************************************/
/******************** mvdens specific functions *******************************/
/******************************************************************************/


size_t get_importance_weight(pmc_simu *psim,mix_mvdens *m,
                             posterior_log_pdf_func *posterior_log_pdf, 
                             void *extra, error **err) {
  
  /* Compute normalized weights from the set of multivariate samples local_X
   local_X should be of size [nsamples_per_proc,ndim]
   local_weights should be [nsamples_per_proc]
   weights should be [nsamples]; only set for master proc
   */
  size_t rr;
  
  rr=generic_get_importance_weight_and_deduced(psim, m, &mix_mvdens_log_pdf_void, posterior_log_pdf, NULL,extra, err);
  forwardError(*err,__LINE__,0);
  return rr;
}


/******************************************************************************/
/******************** importance functions ************************************/
/******************************************************************************/


size_t generic_get_importance_weight_and_deduced(pmc_simu *psim, void *proposal_data,
                             posterior_log_pdf_func *proposal_log_pdf,
                             posterior_log_pdf_func *posterior_log_pdf,
                             retrieve_ded_func *retrieve_ded, void *target_data, error **err) {
   size_t res;

   res = generic_get_importance_weight_and_deduced_verb(psim, proposal_data, proposal_log_pdf, posterior_log_pdf,
							retrieve_ded, target_data, -1, 0, err);
   forwardError(*err, __LINE__, 0.0);
   return res;
}

/* 
 * this is the generic function
 */
size_t generic_get_importance_weight_and_deduced_verb(pmc_simu *psim, void *proposal_data,
						      posterior_log_pdf_func *proposal_log_pdf,
						      posterior_log_pdf_func *posterior_log_pdf,
						      retrieve_ded_func *retrieve_ded, void *target_data, 
						      double t_iter_tempered, int quiet, error **err)
{
  double *x;
  size_t i;
  size_t ndim,n_ded;
  
  size_t cnt;
  double *local_weights;
  double *local_X;
  double *local_R;
  double *local_X_ded;
  short *flg;
  size_t nsamples_per_proc;
  double print_step;
  int last_print_step;
  double MW,MR;
  double rloc, post;
  char preamble[200];

  /* MKDEBUG: To be removed */
  testErrorRetVA(t_iter_tempered<0, pmc_infinite, "t_iter_tempered=%g should be > 0", *err, __LINE__, -1, t_iter_tempered);
  if (t_iter_tempered < 0) {
    t_iter_tempered = 1;  /* No tempering */
  }

  printf("MKDEBUG t_iter_tempered = %g\n", t_iter_tempered);
  
  local_weights     = psim->weights;
  local_X           = psim->X;
  local_X_ded       = psim->X_ded;
  local_R           = psim->log_rho;
  flg               = psim->flg;
  nsamples_per_proc = psim->nsamples;
  n_ded             = psim->n_ded;
  
  if (psim->prop_print_step>0)
    print_step = 100.0/psim->prop_print_step/nsamples_per_proc;
  else
    print_step = 1000;

  last_print_step   = 1;
  if (psim->mpi_size<=0) {
    preamble[0]='\0';
  } else{
    sprintf(preamble,"%d/%d : ",psim->mpi_rank,psim->mpi_size);
  }
  
  ndim  = psim->ndim;
  cnt   = 0;
  MR = MW = -1.0e30;

  for (i=0;i<nsamples_per_proc;i++) {
    //init and set flag to 0
    x=&(local_X[i*ndim]);
    flg[i]=0;
    local_weights[i]=0;
    
    /* Compute log(density) according to proposal */
    rloc = proposal_log_pdf(proposal_data, x, err);
    forwardErrorNoReturn(*err,__LINE__);
    ParameterErrorVerb(*err,x,quiet,ndim); // test if there was an error, print it out and continue

    testErrorVA((isnan(rloc)!=0 || isinf(rloc)!=0), pmc_infinite,
                "Invalid (%g) return value for proposal density", *err, __LINE__, rloc);
    //ParameterErrorVerb(*err,x,quiet,ndim); // test if there was an error, print it out and continue   
    
    if ((i==0) || rloc>MR)
      MR=rloc;
    local_R[i] = rloc;
    
    
    /* Compute log(weight) = log(posterior) - log(proposal) */
    post = posterior_log_pdf(target_data, x, err);
    if (getErrorValue(*err)==tls_cosmo_par) {
       forwardError(*err, __LINE__, 0.0);
    } else {
       forwardErrorNoReturn(*err,__LINE__);
       ParameterErrorVerb(*err,x,quiet,ndim);
    }
    testErrorVA((isnan(rloc)!=0 || isinf(rloc)!=0), pmc_infinite,
                "Invalid (%g) return value for target density", *err, __LINE__, rloc);
    //ParameterError(*err,x,ndim);    

    // MKDEBUG: New, tempering exponent, posterior^t_iter_tempered
    rloc = t_iter_tempered * post - rloc;

    if (retrieve_ded!=NULL && n_ded>0) {
       retrieve_ded(target_data, &(local_X_ded[i*n_ded]), err);
       forwardErrorNoReturn(*err,__LINE__);
       ParameterErrorVerb(*err,x,quiet,ndim); 
       // if the computation of deduced parameter fail, I just want to flag the data and keep going 
    }

    if ((i==0) || rloc>MW)
      MW = rloc;
    local_weights[i] = rloc;

    /* Everything fine, set flag to one and continue */
    flg[i]  = 1;
    cnt    += 1;

    if (i*print_step>=last_print_step) {
       if (! quiet) fprintf(stderr,"%sDone %d%%\n",preamble,last_print_step*psim->prop_print_step);
      last_print_step++;

    }
  }

  //remember that we only computed log 
  psim->maxW = MW;
  psim->maxR = MR;
  
  psim->isLog = 1;
  return cnt;
} 

/* 
 * Some specialization to mvdens
 */
size_t generic_get_importance_weight(pmc_simu *psim, void *m,
                                     posterior_log_pdf_func *proposal_log_pdf,
                                     posterior_log_pdf_func *posterior_log_pdf,
                                     void *extra, error **err) {
  
  size_t rr;
  
  rr=generic_get_importance_weight_and_deduced(psim, m, proposal_log_pdf, posterior_log_pdf, NULL,extra,err);
  forwardError(*err,__LINE__,0);
  return rr;
}


size_t get_importance_weight_plus_ded(pmc_simu *psim, mix_mvdens *m,
                                      posterior_log_pdf_func *posterior_log_pdf, 
                                      retrieve_ded_func *retrieve_ded, void *extra, error **err) {
  size_t rr;
  
  rr=generic_get_importance_weight_and_deduced(psim, m, &mix_mvdens_log_pdf_void, posterior_log_pdf,
                                               retrieve_ded, extra, err);
  forwardError(*err,__LINE__,0);
  return rr;
}

/******************************************************************************/
/******************** processing of importance weights ************************/
/******************************************************************************/

/* ============================================================ *
 * Normalizes the importance weights, original values are       *
 * replaced. Returns the sum of weights. If psim->isLog=0, the  *
 * weights are assumed to be already normalized, and the sum of *
 * (unnormalized) weights stored in psim->logSum.		*
 * After calling this function, weights are *real* importance   *
 * weights, not the log.					*
 * ============================================================ */
double normalize_importance_weight(pmc_simu *psim,error **err) {
  double weight_max, weight_sum;
  size_t i;
  double *weights, wtmp;
  short *flg;
  size_t nsamples;
  int count;
  
  if(psim->isLog==0) {
    fprintf(stderr, "normalize_importance_weight: isLog=0, returning %g\n", psim->logSum);
    return psim->logSum;
  }
  weights  = psim->weights;
  flg      = psim->flg;
  nsamples = psim->nsamples;
  
  /* Get max of weights over samples */
  weight_max = psim->maxW;
  
  /* Compute sum of exp(log(weights)) and store exp(log(weights))/renorm in psim */
  weight_sum = 0.0;
  count = 0;
  for (i=0;i<nsamples;i++) {
    if (flg[i]==1) {
      wtmp = weights[i];
      weights[i] = exp(wtmp-weight_max+DYNMAX);
      if (!finite(weights[i])) {
	 fprintf(stderr, "Got infinite weight w=%g (log w=%g, max=%g) at %ld, flagging it out\n",
		 weights[i], wtmp, weight_max, i);
        flg[i] = 0;
        continue;
      }

      weight_sum += weights[i];
      count++;
     }
  }


  testErrorRet(count==0, pmc_nosamplep,
	       "Not a single point in sample with flag=ok.\n"
	       "Check your likelihood function.", *err, __LINE__, 0);
  testErrorRetVA(weight_sum<=0, pmc_negWeight, "Sum of weights <= 0 (got %g)", *err, __LINE__, 0,weight_sum);
  
  /* Renormalize weights */
  for (i=0;i<nsamples;i++) {
    weights[i] /= weight_sum;
  }
  
  /* Store logSum */
  psim->logSum=log(weight_sum)+weight_max-DYNMAX;

  /* We now have the real weights, up to a renormalization */
  psim->isLog = 0;

  return exp(psim->logSum);
}

int* filter_histogram_init(int nbins,int clipleft,int clipright, error **err) {
  int *self;
  self = malloc_err(sizeof(int)*3,err);
  forwardError(*err,__LINE__,NULL);
  self[0] = nbins;
  self[1] = clipleft;
  self[2] = clipright;
  return self;
}

void filter_histogram(pmc_simu *psim ,void* filter_data, error **err) {
  int *fil,nbin,cl,cr;
  int i;
  double minW, maxW;
  double step,w;
  double bl,br;
  int nfilt;
  
  _DEBUGHERE_("%s","");
  fil = filter_data;
  nbin = fil[0];
  cl = fil[1];
  cr = fil[2];
  _DEBUGHERE_("%d %d %d",nbin,cl,cr);

  testErrorRet(psim->isLog==0,pmc_isLog,"can only filter un-normalized data",*err,__LINE__,);
  // find min and max
  _DEBUGHERE_("%s","");
  minW = 0;
  maxW = 0;
  nfilt = 0;
  _DEBUGHERE_("%s","");
  for(i=0;i<psim->nsamples;i++) {
    if (psim->flg[i]==1) {
      minW = psim->weights[0];
      maxW = psim->weights[0];
      break;
    }
  }
  _DEBUGHERE_("%d",i);
  
  for(;i<psim->nsamples;i++) {
    if (psim->flg[i]==1) {
      w = psim->weights[i];
      if (w>maxW) {
        maxW = w;
      }
      if (w<minW) {
        minW = w;
      }
    }
  }
  _DEBUGHERE_("%s","");
  step = (maxW - minW)*1./nbin;
  bl = minW + step*cl;
  br = maxW - step*cr;
  _DEBUGHERE_("::: %g %g | %g %g %g",minW,maxW,step,bl,br);

  for(i=0;i<psim->nsamples;i++) {
    if (psim->flg[i]==1) {
      w = psim->weights[i];
      if (w<bl || br<w) {
        _DEBUGHERE_("::filt: %g %g | %g %g | %g",minW,maxW,bl,br,w);
        psim->flg[i]=0;
        nfilt++;
      }
    }
  }  
  _DEBUGHERE_("Filtering out %d points",nfilt);
  
  
}

void filter_importance_weight_sig(pmc_simu *psim,double nsig,error **err) {
  size_t i;
  double std,mean,max=0.0;
  size_t cnt;
  
  //fprintf(stderr,"IN FILTERING 1\n");
  if(psim->isLog==1) {
    normalize_importance_weight(psim, err);
    forwardError(*err,__LINE__,);
  }
  //fprintf(stderr,"IN FILTERING 2\n");
  
  // compute std
  std=0;
  mean=0;
  cnt=0;
  for(i=0;i<psim->nsamples;i++) {
    if (psim->flg[i]==0)
      continue;
    mean+=psim->weights[i];
    cnt++;
  }
  //fprintf(stderr,"\n");
  mean*=1./cnt;
  //fprintf(stderr,"mean %g %d\n",mean,cnt);
  
  for(i=0;i<psim->nsamples;i++) {
    if (psim->flg[i]==0)
      continue;
    std+=(psim->weights[i]-mean)*(psim->weights[i]-mean);
  }
  std*=1./(cnt-1);
  //fprintf(stderr,"std %g",std);
  std=sqrt(std);
  //fprintf(stderr," %g",std);
  std*=nsig;
  //fprintf(stderr," %g\n",std);
  
  // flag big values
  for(i=0;i<psim->nsamples;i++) {
    //fprintf(stderr,"%d |",i);
    if (psim->flg[i]==0)
      continue;
    //fprintf(stderr,"go ");
    if (fabs(psim->weights[i]-mean)>nsig) {
      psim->flg[i]=0;
      fprintf(stderr,"Flagging %zu %g %g\n ",i,fabs(psim->weights[i]-mean),nsig);
      print_step(stderr, 0, psim->weights[i],psim->ndim, &(psim->X[i*psim->ndim]));
      continue;
    }
    psim->weights[i]=log(psim->weights[i]);  
    if ((i==0) || (psim->weights[i]>max))
      max=psim->weights[i];
  }
  psim->isLog=1;
  psim->maxW=max;
}

void filter_importance_weight(pmc_simu *psim,double nsig,error **err) {
  size_t i,ip;
  double *weights;
  size_t *sind;
  //int deall;
  //double median;
  
  if(psim->isLog!=1) {
    fprintf(stderr,"Filtering normalized importance weights... It is probably too late...\n");
    weights = malloc_err(sizeof(double)*psim->nsamples,err);
    forwardError(*err,__LINE__,);
    //deall=1;
    for(i=0;i<psim->nsamples;i++) {
      if (psim->flg[i]==0) {
        continue;
      }
      weights[i]=log(psim->weights[i]);
    }
  } else {
    weights=psim->weights;
    //deall=0;
  }
  
  // sort table
  sind = malloc_err(sizeof(size_t)*psim->nsamples,err);
  forwardError(*err,__LINE__,)
  ip=0;
  for(i=0;i<psim->nsamples;i++) {
    if (psim->flg[i]==0) {
      continue;
    }
    sind[ip]=ip;
    ip++;
  }
  
  sort_indices(weights, ip, sind, NULL, 1, err);
  forwardError(*err,__LINE__,);
  
  // get median
  /*
  if (ip%2==0) {
    median = .5 * (weights[sind[ip/2]]+weights[sind[ip/2-1]]);
  } else {
    median = weights[sind[ip/2]];
  }
  */

  // var of 95% of points
  
  
}
size_t* sort_indices(double *list, size_t size, size_t *in_ind,size_t *buf,int order,error **err) {
  double pivot,U,L,M;
  size_t ir,il,i,ipivot,Ui,Li,Mi;
  size_t *buff,*in_indices;
  int del_buf;
  
  // init in_indices if needed
  if (in_ind==NULL) {
    in_indices=malloc_err(sizeof(size_t)*size,err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<size;i++) {
      in_indices[i]=i;
    }
  } else {
    in_indices=in_ind;
  }
  
  /*fprintf(stderr,"--> (%d)> ",size);
  for(i=0;i<size;i++) {
    fprintf(stderr,"%d, ",(int)in_indices[i]);
  }
  fprintf(stderr,"\n");
  */

  // special cases
  if (size<=1) {
    return in_indices;
  }

  if (size==2) {
    if (order*(list[in_indices[0]]-list[in_indices[1]])>0) {
      ir=in_indices[1];
      in_indices[1]=in_indices[0];
      in_indices[0]=ir;
    }
    return in_indices;
  }

  if (size==3) {
    Ui=in_indices[0];
    U=list[Ui];
    Li=in_indices[1];
    L=list[Li];
    Mi=in_indices[2];
    M=list[Mi];
    if (order*(L-U)>0) {
      pivot=L;
      ir=Li;
      L=U;
      Li=Ui;
      U=pivot;
      Ui=ir;
    }
    if (order*(M-L)>0) {
      if (order*(M-U)>0) {
        ir=Ui;
        Ui=Mi;
        Mi=ir;
      }
    } else {
      ir=Li;
      Li=Mi;
      Mi=ir;
    }
    in_indices[0]=Li;
    in_indices[1]=Mi;
    in_indices[2]=Ui;
    return in_indices;
  }
  
  del_buf=0;
  if (buf==NULL) {
    del_buf=1;
    buff=malloc_err(sizeof(size_t)*size,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    buff=buf;
  }
  
  // look for pivot value
  Li = in_indices[0];
  L = list[Li];
  Ui = in_indices[1];
  U = list[Ui];
  il = 2;
  
  while ((order*(L-U)<0) && (il<size)) {
    Li = Ui;
    L = U;
    Ui = in_indices[il];
    U = list[Ui];
    il++;
  }
  
  if ((order*(L-U)<0) && (il==size)) {
    return in_indices;
  }
  
  pivot=L;
  ipivot=Li;
  in_indices[il-2]=Ui;
  
  il=il-1;
  ir=0;
  for(i=il+1;i<size;i++) {
  
    Mi=in_indices[i];
    if (order*(list[Mi]-pivot)<0) {
      in_indices[il]=Mi;
      il++;
    } else {
      buff[ir]=Mi;
      ir++;
    }
    testErrorRet((il+ir+1>size),pmc_sort,"Fatal error in sort",*err,__LINE__,NULL);
  }
  testErrorRet((il+ir+1!=size),pmc_sort,"Fatal error in sort",*err,__LINE__,NULL);
  
  in_indices[il]=ipivot;
  /*fprintf(stderr,"pivot %d %g\n",ipivot,pivot);
  fprintf(stderr,"out -> ");*/
  for(i=0;i<ir;i++) {
    in_indices[il+i+1]=buff[i];
    //fprintf(stderr,"%d %g, ",buff[i],list[buff[i]]);
  }
  /*
  fprintf(stderr,"\n");
  fprintf(stderr,"in -> ");
  for(i=0;i<il;i++) {
    fprintf(stderr,"%d %g, ",in_indices[i],list[in_indices[i]]);
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"tot -> ");
  for(i=0;i<size;i++) {
    fprintf(stderr,"%d %g, ",in_indices[i],list[in_indices[i]]);
  }
  fprintf(stderr,"\n");
  */
  
  if (il>1) {
    sort_indices(list,il,in_indices,buff,order,err);
    forwardError(*err,__LINE__,NULL);
  }
  if (ir>1) {
    sort_indices(list,ir,&in_indices[size-ir],buff,order,err);
    forwardError(*err,__LINE__,NULL);
  }
  
  if (del_buf==1) {
    free(buff);
  }
  
  return in_indices;
}

/* ============================================================ *
 * Returns the perplexity and calculates the effective sample   *
 * size (ESS).							*
 * ============================================================ */
#define perp_zero 1e-20
double perplexity_and_ess(pmc_simu *psim, int normalize, double *ess, error **err) {
  size_t i,nexit;
  double perp,wi;
  double *weights;
  short* flg;
  size_t nsamples;
  
  weights  = psim->weights;
  flg      = psim->flg;
  nsamples = psim->nsamples;
  
  if (normalize==MC_NORM) {
    normalize_importance_weight(psim, err);
    forwardError(*err,__LINE__,0);
  }
  perp  = 0;
  nexit = 0;
  *ess  = 0;
  for(i=0;i<nsamples;i++) {
    if (flg[i]==0) {
      nexit++;
      continue;
    }
    wi = weights[i];
    if (wi<perp_zero)
      continue;
    perp += wi*log(wi);
    *ess += wi*wi;
    //fprintf(stderr,"%g %g %g %g\n",nrm,wi,log(wi),perp);
  }
  testErrorRet(nexit==nsamples,pmc_undef,"all points are flagged",*err,__LINE__,0);
  
  testErrorRet(ess==0,pmc_infinite,"ess is infinite",*err,__LINE__,0);
  *ess = 1/(*ess);
  //fprintf(stderr,"%g\n",perp);
  return exp(-perp)/(nsamples-nexit);
}
#undef perp_zero

double perplexity(pmc_simu *psim, int normalize,error **err) {
  double ess,perp;
  perp = perplexity_and_ess(psim,normalize,&ess, err);
  forwardError(*err,__LINE__,0);
  return perp;
}

/* Returns Bayesian evidence E. In ln_evi the logarithm of E is stored */
double evidence(pmc_simu *psim, double *ln_evi,error **err) {
  int i, n;
  double lne;
  
  for (i=0,n=0; i<psim->nsamples; i++) {
    if (psim->flg[i]==0) continue;
    n++;
  }

  testErrorRet(n==0,pmc_undef,"All points are flagged",*err,__LINE__,0);

  lne = psim->logSum - log((double)n);    
  
  if (ln_evi!=NULL) {
    *ln_evi =  lne;
  }
  return exp(lne);
}

/* Returns mean of a-th parameter */
double mean_from_psim(double *X, double *w, short *flg, int nsamples, int ndim, int a)
{
  double m;
  int i;
  
  /* TODO(?): check isLog */
  for (i=0,m=0.0; i<nsamples; i++) {
    if (flg[i]==0) continue;
    m += w[i]*X[i*ndim+a];
  }
  
  return m;
}

int double_cmp(const void *av, const void *bv)
{
  const double *a, *b;
  
  a = (const double*)av;
  b = (const double*)bv;
  if (*a<*b) return -1;
  if (*a>*b) return +1;
  return 0;
}

/* Sets the highest n weights to zero. Does not normalizes to one. */
void clip_weights(pmc_simu *psim, int n, FILE *OUT, error **err)
{
  double *max_weights, sum_w;
  int i, j;
  
  if (n==0) return;
  
  max_weights = calloc_err(n, sizeof(double), err);  forwardError(*err, __LINE__,);
  
  /* Initialise temp array */
  i = j = 0;
  while (i<n && j<psim->nsamples) {
    if (psim->flg[i]==0) {
      j++;
      continue;
    }
    max_weights[i] = psim->weights[j];
    i++;
  }

  qsort(max_weights, n, sizeof(double), double_cmp);
  
  /* Store n largest weights */
  for (j=n; j<psim->nsamples; j++) {
    if (psim->flg[j]==0) continue;
    if (max_weights[0]<psim->weights[j]) {
      max_weights[0] = psim->weights[j];
      qsort(max_weights, n, sizeof(double), double_cmp);
    }
  }
  
  /* Set to zero. TODO: Replace with Karim's method */
  if (OUT) fprintf(OUT, "clip_weights: setting %d weights larger than %g to zero:\n", n, max_weights[0]);
  for (j=0,sum_w=0.0; j<psim->nsamples; j++) {
    if (psim->flg[j]==0) continue;
    if (psim->weights[j]>=max_weights[0]) {
      if (OUT) fprintf(OUT, "%g ", psim->weights[j]);
      psim->weights[j] = 0.0;
      psim->flg[j] = 0;
    } else {
      sum_w += psim->weights[j];
    }
  }
  if (OUT) fprintf(OUT, "\n");
  
  testErrorRetVA(sum_w<=0, pmc_infnan, "Sum of weights (%g) has to be larger than zero",
		 *err, __LINE__,, sum_w);
  /* Normalize */
  for (j=0; j<psim->nsamples; j++) {
    if (psim->flg[j]==0) continue;
    psim->weights[j] /= sum_w;
  }
  
  free(max_weights);
}

/******************************************************************************/
/******************** updates (mvdens) ****************************************/
/******************************************************************************/

void update_prop(mix_mvdens *m, pmc_simu *psim, error **err)
{
  /* Recompute MOG parameters (weights, means and covariances) from samples
   and associated normalized IS weights */
  size_t ndim, ncomp;
  size_t *count;
  double w8,gw8,gamma_d;
  size_t i,j,k,l,index,insamples;
  
  double minwgt;
  double *X,*gamma_rho_d,*weight_gamma_d;
  size_t *indices; 
  double *weights;
  short *flg;
  size_t nsamples;
  double rloc,vloc;
  double *buf;
  
  X=psim->X;
  indices=psim->indices; 
  weights=psim->weights;
  flg=psim->flg;
  nsamples=psim->nsamples;
  
  ncomp=m->ncomp;
  ndim=m->ndim;
  minwgt=1./nsamples;
    
  gsl_set_error_handler_off();
  
  count=calloc_err(sizeof(long),ncomp,err);
  forwardError(*err,__LINE__,);

  gamma_rho_d = malloc_err(sizeof(double)*nsamples,err);
  forwardError(*err,__LINE__,);
  
  weight_gamma_d = calloc_err(sizeof(double), ncomp,err);
  forwardError(*err,__LINE__,);

  prepare_update(m, count,psim,&buf,err);
  forwardError(*err,__LINE__,);

  /* Compute the updates, eq 8 of Cappe & Robert */
  /* First count occurences of indices in samples  (compute indicatrices)*/
  gsl_vector_set_zero(m->wght_view);
  for (i=0;i<ncomp;i++) {
    gsl_vector_set_zero(m->comp[i]->mean_view);
    gsl_matrix_set_zero(m->comp[i]->std_view);
  }
  
  // update means
  for (i=0;i<nsamples;i++) {
    if (flg[i]==0)
      continue;
    insamples=i*ndim;
    index=indices[i];
    if (count[index] > MINCOUNT) {
      w8 = weights[i];
      if(m->comp[index]->df!=-1) {
        gamma_d= (m->comp[index]->df + ndim) /
          (m->comp[index]->df + 
            scalar_product(m->comp[index], X+insamples, err));
        forwardError(*err,__LINE__,);
        gw8=gamma_d*w8;
      } else {
        gw8=w8;
      }
      gamma_rho_d[i]=gw8;
      m->wght[index]+=w8;
      weight_gamma_d[index] += gw8;
      for(j=0;j<ndim;j++) {
        m->comp[index]->mean[j] += X[insamples+j]*gw8;
      }
    }
  }
  
  for (i=0;i<ncomp;i++) {
    if (count[i] > MINCOUNT) {
      gsl_vector_scale(m->comp[i]->mean_view,1.0/weight_gamma_d[i]);
    }
  }
  
  // Now update the covariances, once means have been updated
  for (i=0;i<nsamples;i++) {
    if (flg[i]==0)
      continue;
    index=indices[i];
    insamples=i*ndim;
    if (count[index] > MINCOUNT) {
      gw8=gamma_rho_d[i];
      for (k=0;k<ndim;k++) {
        vloc=(X[insamples+k] - m->comp[index]->mean[k]);
        m->comp[index]->std[k*ndim+k]+=vloc*vloc*gw8;
        for (l=k+1;l<ndim;l++) {
          if(abs(k-l)>=m->comp[index]->band_limit)
            break;
          rloc=vloc*(X[insamples+l] - m->comp[index]->mean[l]) * gw8;
          m->comp[index]->std[k*ndim+l] += rloc;
          m->comp[index]->std[l*ndim+k] += rloc;
        }
      }
    }
  }
  
  for (i=0;i<ncomp;i++) {
    if (count[i] > MINCOUNT) {
      gsl_matrix_scale(m->comp[i]->std_view,1.0/m->wght[i]);
    }
  }
  
  cleanup_after_update(m,count,minwgt,psim,buf,err);
  forwardError(*err,__LINE__,);
  
  free(count);
  free(weight_gamma_d);
  free(gamma_rho_d);
}

/* ============================================================ *
 * Rao-Blackwell version of proposal update.			*
 * Recomputes mixture parameters (weights, means and covs) from *
 * samples and associated normalized IS weights.		*
 * ============================================================ */
void update_prop_rb_void(void *m, pmc_simu *psim, error **err)
{
   update_prop_rb((mix_mvdens*)m, psim, err);
   forwardError(*err, __LINE__,);
}

void update_prop_rb(mix_mvdens *m, pmc_simu *psim, error **err)
{
  size_t ndim, ncomp;
  size_t *count;
  double w8, gw8;
  size_t i,j,k,l,index,insamples;
  
  double* gamma_rho_d,*mean_d,*weight_d;
  double *X;
  //size_t *indices; 
  double *weights;
  short *flg;
  size_t nsamples;
  double minwgt,gamma_d,rho_d,rho_d0;
  double *weight_gamma_d;
  double vloc,rloc;
  double * log_rho;
  double *buf;
  
  X        = psim->X;
  log_rho  = psim->log_rho;
  //indices  = psim->indices; 
  weights  = psim->weights;
  flg      = psim->flg;
  nsamples = psim->nsamples;
  minwgt   = 1.0/nsamples;

  ncomp    = m->ncomp;
  ndim     = m->ndim;

  gsl_set_error_handler_off();

  count=calloc_err(sizeof(long),ncomp,err);
  forwardError(*err,__LINE__,);

  gamma_rho_d=calloc_err(sizeof(double), ncomp*nsamples,err);
  forwardError(*err,__LINE__,);

  weight_gamma_d = calloc_err(sizeof(double), ncomp,err);
  forwardError(*err,__LINE__,);
  
  weight_d = calloc_err(sizeof(double), ncomp,err);
  forwardError(*err,__LINE__,);
  
  mean_d = calloc_err(sizeof(double), ncomp*ndim,err);
  forwardError(*err,__LINE__,);
  
  
  prepare_update(m, count,psim,&buf,err);
  forwardError(*err,__LINE__,);
  
  
  // update means
  for (i=0;i<nsamples;i++) {
    if (flg[i]==0)
      continue;
    insamples=i*ndim;
    rho_d0 = log_rho[i]-psim->maxR+DYNMAX;
    for(index=0;index<ncomp;index++) {
      if (count[index] > MINCOUNT) {
	 /* Compute rho_d */
        rho_d = exp(mvdens_log_pdf(m->comp[index], X+insamples, err)-rho_d0);
        forwardError(*err,__LINE__,);
        rho_d*=m->wght[index];
        w8 = weights[i]*rho_d;
        if(m->comp[index]->df!=-1) {
          gamma_d= (m->comp[index]->df + ndim) /
            (m->comp[index]->df + scalar_product(m->comp[index], X+insamples, err));
          forwardError(*err,__LINE__,);
          gw8=gamma_d*w8;
        } else {
          gw8=w8;
        }
        gamma_rho_d[i*ncomp+index]=gw8;
        weight_d[index]+=w8;
	// if (index==0) fprintf(stderr, "### %d %f %f %f\n", i, weight_d[index], weights[i], rho_d);
        weight_gamma_d[index] += gw8;
        for(j=0;j<ndim;j++) {
          mean_d[index*ndim+j] += X[insamples+j]*gw8;
        }
      }	
    }
  }
  
  // clean up std
  // copy weights, rescale means;
  for (i=0;i<ncomp;i++) {
    m->wght[i]=weight_d[i];
    /* This should be useless 
    testErrorRet(!finite(m->wght[i]), pmc_infnan,
		 "Proposal weight not finite. Check whether your likelihood function catches\n"
		 "all nans and infs\n", *err, __LINE__,); */
    //fprintf(stderr, "%d -> %g\n",weight_d[i]);
    gsl_matrix_set_zero(m->comp[i]->std_view);
    if (count[i] > MINCOUNT) {
      for(j=0;j<ndim;j++) {
        m->comp[i]->mean[j]=mean_d[i*ndim+j]/weight_gamma_d[i];
      }
    } else {
      for(j=0;j<ndim;j++) {
        m->comp[i]->mean[j]=0;
      }
    }
  }
  free(mean_d);
  free(weight_d);
  free(weight_gamma_d);
  
  // Now update the covariances, once means have been updated
  for (i=0;i<nsamples;i++) {
    if (flg[i]==0)
      continue;
    insamples=i*ndim;
    for(index=0;index<ncomp;index++) {
      if (count[index] > MINCOUNT) {
        gw8=gamma_rho_d[index+i*ncomp];
        for (k=0;k<ndim;k++) {
          vloc=(X[insamples+k] - m->comp[index]->mean[k]);
          //fprintf(stderr,"%d %d %d-> %g %g %g -> ",index,i,k,m->comp[index]->mean[k],vloc,gw8);
          m->comp[index]->std[k*ndim+k]+=vloc*vloc*gw8;
          for (l=k+1;l<ndim;l++) {
            if(abs(k-l)>=m->comp[index]->band_limit)
              break;
            rloc=vloc*(X[insamples+l] - m->comp[index]->mean[l]) * gw8;
            m->comp[index]->std[k*ndim+l] += rloc;
            m->comp[index]->std[l*ndim+k] += rloc;
          }
        }
      }
    }
  }
  
  for (i=0;i<ncomp;i++) {
    if (count[i] > MINCOUNT) {
      // fprintf(stderr, "*** %d %d %f\n", i, count[i], m->wght[i]);
      // testErrorRet(1 || m->wght[i]<=0, pmc_infnan, "Division by zero", *err, __LINE__,);
      // BEWARE THIS IS WITH m->wght and not weight_gamma !
      gsl_matrix_scale(m->comp[i]->std_view,1.0/m->wght[i]);
    }
  }
  
  cleanup_after_update(m,count,minwgt,psim,buf,err);
  forwardError(*err,__LINE__,);
  
  free(count);
  free(gamma_rho_d);
}

void cleanup_after_update(mix_mvdens *m,size_t *count,double minwgt,pmc_simu *psim,double *buf,error **err) {
  double wght_sum;
  size_t i,ncomp;
  int nok;
  
  nok=0;
  ncomp=m->ncomp;
  
  // Make sure mixture weights sum to one (first time)
  wght_sum=0.0;
  for (i=0;i<ncomp;i++) {
    if (count[i] > MINCOUNT) {
      wght_sum += m->wght[i];
    }
  }
  gsl_vector_scale(m->wght_view,1.0/wght_sum);
  
  // Now cholesky decompose all valid covariances
  for (i=0;i<ncomp;i++) {
    if (m->wght[i] > minwgt) {
      m->comp[i]->chol=0;
      //fprintf(stdout,"component %d (%g)\n",i,m->wght[i]);
      //mvdens_print(stdout,m->comp[i]);
      mvdens_cholesky_decomp(m->comp[i],err);
      if isError(*err) {
	fprintf(stderr, "Bogus component %zu\n", i); fflush(stderr);
        mvdens_dump(stderr, m->comp[i]);
        fprintf(stdout, "Component %zu:\n", i);
        printError(stdout, *err);
        purgeError(err);
        if (psim->retry==0) {
          memcpy(m->comp[i]->std,buf+i*m->ndim*m->ndim,sizeof(double)*m->ndim*m->ndim); 
          nok=1;
        } else {
          m->wght[i]=0;
        }
        m->comp[i]->chol=1;
      }
    }
  }
  psim->retry=nok;
  if (psim->retry!=0) {
    fprintf(stderr,"RETRY state -> %d\n",psim->retry);
  }
  // Get rid of unwanted components
  for (i=0;i<ncomp;i++) {
    if (count[i] < MINCOUNT || m->wght[i] < minwgt) {
      m->wght[i]=0.0;
      gsl_vector_set_zero(m->comp[i]->mean_view);
      gsl_matrix_set_zero(m->comp[i]->std_view);
    }
  }
  
  // Renormalize after cleanup
  wght_sum=0.0;
  for (i=0;i<ncomp;i++) {
    wght_sum += m->wght[i];
  }

  testErrorRetVA(wght_sum<=0.0, pmc_negWeight, "Sum of weight is %g", *err, __LINE__,, wght_sum);

  gsl_vector_scale(m->wght_view,1.0/wght_sum);
}

void prepare_update(mix_mvdens *m, size_t *count,pmc_simu *psim,double** pbuf,error **err) {
  size_t i,index,nsamples;
  size_t *indices;
  short *flg;
  double* buf = NULL;
  
  /* Check psim status */
  if (psim->isLog==1) {
    normalize_importance_weight(psim,err);
    forwardError(*err,__LINE__,);
  }
  
  //ncomp=m->ncomp;
  indices = psim->indices;
  flg=psim->flg;
  nsamples=psim->nsamples;
  
  /* Reinitialize all MOG parameters to zero */
  m->init_cwght=0.0;
  
  /* Compute the updates, eq 12 of Cappe & Robert */
  /* First count occurences of indices in samples */
  
  // compute count
  for (i=0;i<nsamples;i++) {
    if (flg[i]==1) {
      index        = indices[i];
      count[index] = count[index]+1;
    }
  }    
  
  mix_mvdens_cholesky_decomp(m,err);
  forwardError(*err,__LINE__,);
  
  if (psim->retry==0) {
    buf=malloc_err(sizeof(double)*m->ncomp*m->ndim*m->ndim,err);
    forwardError(*err,__LINE__,);
    
    for(i=0;i<m->ncomp;i++) {
      memcpy(buf + i*m->ndim*m->ndim,m->comp[i]->std,sizeof(double)*m->ndim*m->ndim);
    }
  }
  (*pbuf)=buf;
  
}


/******************************************************************************/
/******************** ressampling functions ***********************************/
/******************************************************************************/
size_t *sample_from_pmc_simu(const pmc_simu *psim, const gsl_rng * rng, error **err) {
  long i;
  gsl_ran_discrete_t *lut;
  size_t *isample;
  
  /* generate lookup table */
  lut     = gsl_ran_discrete_preproc(psim->nsamples, psim->weights);
  isample = malloc_err(sizeof(size_t)*psim->nsamples, err);  
  forwardError(*err, __LINE__, NULL);
  
  for (i=0; i<psim->nsamples; i++) {
    isample[i] = gsl_ran_discrete(rng, lut);
  }
  
  gsl_ran_discrete_free(lut);
  return isample;
}

/* Resamples a pmc simulation (residual resampling) and returns an array of the multiplicity for each point. *
 * See Douc, Cappe & Moulines 2008 */
unsigned int *resample_residual(const pmc_simu *psim, const gsl_rng *rng, error **err)
{
   double *weights, R;
   unsigned int i, Ntrials, Nsamples, *n;

   Ntrials  = psim->nsamples; /* This could be a different number */
   Nsamples = psim->nsamples; 
   n        = malloc_err(sizeof(unsigned int)*Ntrials, err);            forwardError(*err, __LINE__, NULL);
   weights  = malloc_err(sizeof(double)*Nsamples, err);                 forwardError(*err, __LINE__, NULL);

   R = exp(psim->logSum);

   for (i=0; i<Nsamples; i++) {
      weights[i] = (Nsamples*psim->weights[i] - (int)(Nsamples*psim->weights[i]))/((double)Nsamples - R);
   }

   gsl_ran_multinomial(rng, psim->nsamples, Ntrials, weights, n);

   free(weights);

   return n;
}

