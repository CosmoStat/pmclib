/*
 *  pmc_mpi.c
 *  likely
 *
 *  Created by Karim Benabed on 13/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef __PLANCK__
#include "HL2_likely/pmclib/pmc_mpi.h"
#else
#include "pmc_mpi.h"
#endif

/******************************************************************************/
/******************** high level interface ************************************/
/******************************************************************************/

pmc_simu* pmc_simu_init_mpi(long nsamples, int ndim, int nded, error **err) {
  int mpi_rank,mpi_size;
  pmc_simu *psim;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); 

  if (mpi_rank==0) {
    psim=pmc_simu_init_plus_ded(nsamples, ndim,nded, err);
  } else {
    psim=pmc_simu_init_plus_ded(1, ndim,nded, err);
  }
  psim->mpi_rank = mpi_rank;
  psim->mpi_size = mpi_size;
  return psim;  
}

void pmc_simu_realloc_mpi(pmc_simu *psim,long newsamples,error **err) {
  if (psim->mpi_rank==0) {
    pmc_simu_realloc(psim,newsamples,err);
    forwardError(*err,__LINE__,);
  }
}


size_t pmc_simu_importance_mpi(pmc_simu *psim, gsl_rng *r,error **err) {
  size_t nok,mysamples,master_samples=0;

  testErrorRet(psim->proposal->data==NULL,pmc_allocate,
               "proposal undefined",*err,__LINE__,0);
  testErrorRet(psim->target->data==NULL,pmc_allocate,
               "target undefined",*err,__LINE__,0);
  
  if (psim->mpi_rank==0) { //master
    
    psim->proposal->simulate(psim,psim->proposal->data,r,psim->pb,err);
    forwardError(*err,__LINE__,0);
    
    mysamples=send_simulation(psim,psim->mpi_size,err);
    forwardError(*err,__LINE__,0);
    
    master_samples=psim->nsamples;
    psim->nsamples=mysamples;
    
  } else {//client
    
    //only receives simu
    receive_simulation(psim,psim->mpi_size,psim->mpi_rank,err);
    forwardError(*err,__LINE__,0);      
    
  }
  
  /**** Compute importance weights ********************************************/
  
  nok = generic_get_importance_weight_and_deduced(psim,
                                                        psim->proposal->data,
                                                        psim->proposal->log_pdf,
                                                        psim->target->log_pdf,
                                                        psim->target->retrieve,
                                                        psim->target->data,err);
  forwardError(*err,__LINE__,0);
  
  /**** cleanup ***************************************************************/  
  if (psim->mpi_rank==0) { //master
    //receive computation from clients
    mysamples=psim->nsamples;
    psim->nsamples=master_samples;      
    nok=receive_importance_weight(psim,psim->mpi_size,nok,mysamples,err);
    forwardError(*err,__LINE__,0);
  } else {
    //send results    
    send_importance_weight(psim->mpi_rank,psim->mpi_size,psim,nok);
  }
  
  return nok;
}

double pmc_simu_pmc_step_mpi(pmc_simu *psim, gsl_rng *r,error **err) {
  double perp;
  
  testErrorRet(psim->pmc_update==NULL,pmc_allocate,
               "update function undefined",*err,__LINE__,0);
  pmc_simu_importance_mpi(psim,r,err);
  forwardError(*err,__LINE__,0);
  perp = 0;
  if (psim->mpi_rank==0) {
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
  }
  
  distribution_broadcast(psim->proposal,err);
  forwardError(*err,__LINE__,0);
  
  return perp;  
}

void pmc_simu_init_classic_importance_mpi(pmc_simu* psim,
                                      distribution *target, parabox *pb,
                                      void *prop_data,int prop_print_step,error **err) {
  distribution *prop;
  // init target
  pmc_simu_init_target(psim,target,pb,err);
  forwardError(*err,__LINE__,);
  
  // init proposal
  prop = mix_mvdens_distribution_mpi(target->ndim, prop_data, err);
  forwardError(*err,__LINE__,);
  
  pmc_simu_init_proposal(psim, prop, prop_print_step, err);
  forwardError(*err,__LINE__,);
}

void pmc_simu_init_classic_pmc_mpi(pmc_simu* psim,
                               distribution *target, parabox *pb,
                               void *prop_data,int prop_print_step,void* filter_data,
                               filter_func* pmc_filter,error **err) {
  
  pmc_simu_init_classic_importance_mpi(psim,target,pb,
				       prop_data,prop_print_step,err);
  forwardError(*err,__LINE__,);
  
  pmc_simu_init_pmc(psim, filter_data, pmc_filter, &update_prop_rb_void);
}

void* mvdens_broadcast_mpi(void* data,error **err) {
  int mpi_rank,mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); 
  void *rdata;
  
  if (mpi_rank==0) {
    send_mix_mvdens(data,mpi_size,err);
    forwardError(*err,__LINE__,NULL);
    return data;
  } else {
    rdata = receive_mix_mvdens(mpi_rank,mpi_size,err);
    forwardError(*err,__LINE__,NULL);
    return rdata;
  }
}

distribution* init_distribution_full_mpi(int ndim,
                                       void* data, 
                                       posterior_log_pdf_func* log_pdf,
                                       posterior_log_free* freef,
                                       simulate_func *simulate,
                                       int nded,
                                       retrieve_ded_func* retrieve,
                                       mpi_exchange_func *broadcast_mpi, error **err) {

  void *rdata;
  distribution* dist;
  
  if (broadcast_mpi!=NULL) {
    rdata = broadcast_mpi(data,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    rdata = data;
  }
  dist = init_distribution_full(ndim,rdata,log_pdf,freef,simulate,nded,retrieve,err);
  forwardError(*err,__LINE__,NULL);
  dist->broadcast_mpi = broadcast_mpi;
  return dist;  
}

distribution* init_distribution_mpi(int ndim,
                                    void* data, 
                                    posterior_log_pdf_func* log_pdf,
                                    posterior_log_free* freef,
                                    simulate_func *simulate,
                                    mpi_exchange_func *broadcast_mpi, error **err) {
  
  distribution* dist;
  
  dist = init_distribution_full_mpi(ndim,data,log_pdf,freef,simulate,0,NULL,broadcast_mpi,err);
  forwardError(*err,__LINE__,NULL);
  return dist;  
}

void distribution_broadcast(distribution *dist,error **err) {
  void *rdata;
  testErrorRet(dist->broadcast_mpi==NULL,pmc_incompat,
               "Cannot broadcast the distribution",*err,__LINE__,);
  
  rdata = dist->broadcast_mpi(dist->data,err);
  forwardError(*err,__LINE__,);
  if (rdata!=dist->data) {
    if (dist->free!=NULL) {
      dist->free(&dist->data);
    }
    dist->data = rdata;
  }
  return;
}

distribution *mix_mvdens_distribution_mpi(int ndim,void *mxmv, error **err) {
  distribution *dist;
  dist = init_distribution_full_mpi(ndim, mxmv, mix_mvdens_log_pdf_void, mix_mvdens_free_void,
                                    simulate_mix_mvdens_void, 0, NULL, mvdens_broadcast_mpi, 
                                    err);
  forwardError(*err,__LINE__,NULL);
  return dist;
}

#ifdef HAS_LUA

mix_mvdens* rc_mix_mvdens_mpi(confFile* rc, error **err) {
  int mpi_rank,mpi_size;
  mix_mvdens* mmv;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); 
  
  if (mpi_rank==0) {
    mmv = rc_mix_mvdens(rc,err);
    forwardError(*err,__LINE__,NULL);
  }
  mmv = mvdens_broadcast_mpi(mmv,err);
  forwardError(*err,__LINE__,NULL);
  
  return mmv;
}

distribution* rcinit_mix_mvdens_mpi(confFile* rc, error **err) {
  int mpi_rank,mpi_size;
  mix_mvdens* mmv;
  distribution *dist;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); 

  mmv = NULL;
  if (mpi_rank==0) {
    mmv = rc_mix_mvdens(rc,err);
    forwardError(*err,__LINE__,NULL);
  }
  mmv = mvdens_broadcast_mpi(mmv,err);
  forwardError(*err,__LINE__,NULL);
  
  dist = init_distribution_full(mmv->ndim, mmv, mix_mvdens_log_pdf_void, mix_mvdens_free_void,
                                    simulate_mix_mvdens_void, 0, NULL, err);
  forwardError(*err,__LINE__,NULL);
  dist->broadcast_mpi = mvdens_broadcast_mpi;
  return dist;
}

pmc_simu *init_importance_mpi_from_rc(confFile* rc,char* root,error **err) {
  pmc_simu *psim;
  long nsamples;
  confFile *rca;
  parabox *pb;
  int prop_print_step;
  distribution* proposal, *target;
  
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
  
  
  psim = pmc_simu_init_mpi(nsamples,target->ndim,target->n_ded,err);
  forwardError(*err,__LINE__,NULL);

  pb = parabox_from_rc(rca,"pb",err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(pb->ndim!=target->ndim,pmc_dimension,"incompatible dim for pb (got %d expected %d)",*err,__LINE__,NULL,
                 pb->ndim,target->ndim);
  
  pmc_simu_init_target(psim,target,pb,err);
  forwardError(*err,__LINE__,);
  
  prop_print_step = rc_safeget_integer(rca,"print_pc",0,err);
  forwardError(*err,__LINE__,);
  
  pmc_simu_init_proposal(psim, proposal, prop_print_step, err);
  forwardError(*err,__LINE__,);
  
  rc_close(&rca);
  
  return psim;
}

pmc_simu* init_pmc_mpi_from_rc(confFile *rc,char *root,error **err) {
  pmc_simu *psim;
    
  psim = init_importance_mpi_from_rc(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  pmc_simu_init_pmc(psim, NULL, NULL, update_prop_rb_void);
  
  return psim;  
}

#endif

/******************************************************************************/
/******************** exchanging pmc data *************************************/
/******************************************************************************/

/* Master sends to clients */
size_t send_simulation(pmc_simu* psim, int nproc, error **err) {
  size_t nsamples_per_proc,nsamples_master,current_sample,buffsize,remaining;
  int tag,i;
  
  tag = nproc;
  nsamples_per_proc = (psim->nsamples/nproc);
  nsamples_master = nsamples_per_proc;
  current_sample = nsamples_master;
  buffsize = nsamples_per_proc;
  remaining=psim->nsamples-nsamples_per_proc*nproc;
  if (remaining>0) {
    buffsize=nsamples_per_proc+1;
    remaining--;
  }
  fprintfDEBUG(stderr,"Master %ld Clients %ld\n",nsamples_master,nsamples_per_proc);
  for (i=1;i<nproc;i++) {
    fprintfDEBUG(stderr,"I am the master, I am sending buffsize to %d on %d (%zu)\n",i,tag*mc_tag_simulation_size+i,buffsize);
    MPI_Send(&buffsize,1,MPI_LONG,i,tag*mc_tag_simulation_size+i,MPI_COMM_WORLD);
    
    fprintfDEBUG(stderr,"I am the master, I am sending sample points to %d on %d \n",i,tag*mc_tag_simulation+i);
    fprintfDEBUG(stderr,"  ndim, current_sample, nsamples = %d %zu %zu\n",psim->ndim,current_sample,psim->nsamples);
    MPI_Send(&(psim->X[current_sample*psim->ndim]),buffsize*psim->ndim,MPI_DOUBLE,i,tag*mc_tag_simulation+i,MPI_COMM_WORLD);
    
    fprintfDEBUG(stderr,"I am the master, I am sending X_ded to %d on %d\n",i,tag*mc_tag_flg+i);
    MPI_Send(&(psim->X_ded[current_sample*psim->n_ded]),buffsize*psim->n_ded,MPI_DOUBLE,i,tag*mc_tag_simulation+i,
             MPI_COMM_WORLD);
    
    fprintfDEBUG(stderr,"I am the master, I am sending flags to %d on %d \n",i,tag*mc_tag_flg+i);
    MPI_Send(&(psim->flg[current_sample]),buffsize,MPI_SHORT,i,tag*mc_tag_flg+i,MPI_COMM_WORLD);
    
    current_sample += buffsize;
    buffsize=nsamples_per_proc;
    if (remaining>0) {
      buffsize=nsamples_per_proc+1;
      remaining--;
    }  
  }
  return nsamples_master;
}

/* Clients receive from master */
void receive_simulation(pmc_simu* psim, int basetag,int myid, error **err)
{
  size_t buffsize;
  int tag;
  
  tag = basetag;
  
  fprintfDEBUG(stderr,"%d:: I am waiting for buffsize on %d\n",myid,tag*mc_tag_simulation_size+myid);
  MPI_Recv(&buffsize,1,MPI_LONG,0,tag*mc_tag_simulation_size+myid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
  fprintfDEBUG(stderr,"%d:: %zu\n",myid,buffsize);
  if (buffsize!=psim->nsamples) {
    pmc_simu_realloc(psim,buffsize,err);
    forwardError(*err,__LINE__,);
  }
  
  fprintfDEBUG(stderr,"%d:: I am waiting for sample points on %d (%zu)\n",myid,tag*mc_tag_simulation+myid,buffsize);
  MPI_Recv(psim->X,buffsize*psim->ndim,MPI_DOUBLE,0,tag*mc_tag_simulation+myid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
  fprintfDEBUG(stderr,"%d:: I am waiting for X_ded on %d (%zu)\n",myid,tag*mc_tag_simulation+myid,buffsize);
  MPI_Recv(psim->X_ded,buffsize*psim->n_ded,MPI_DOUBLE,0,tag*mc_tag_simulation+myid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  fprintfDEBUG(stderr,"%d:: I am waiting for flags on %d (%zu)\n",myid,tag*mc_tag_flg+myid,buffsize);
  MPI_Recv(psim->flg,buffsize,MPI_SHORT,0,tag*mc_tag_flg+myid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  fprintfDEBUG(stderr,"%d:: received simu on %d\n",myid,tag*mc_tag_simulation+myid);
  
  return;
}

/* Clients send to master */
void send_importance_weight(int myid, int basetag, pmc_simu *psim, size_t nok)
{
  // Send unnormalized local weights to master
  size_t buffsize;

  fprintfDEBUG(stderr, "%d start send_importance_weight\n", myid);
  buffsize=psim->nsamples;
  MPI_Send(&nok,1,MPI_LONG,0,basetag*mc_tag_nok+myid,MPI_COMM_WORLD);
  MPI_Send(&buffsize,1,MPI_LONG,0,basetag*mc_tag_data_size+myid,MPI_COMM_WORLD);
  if (nok==0) {
    fprintfDEBUG(stderr,"%d Argl\n",myid);
    return;
  }
  
  fprintfDEBUG(stderr,"%d Sending buffsize (%zu)\n",myid,buffsize);
  MPI_Send(psim->weights,buffsize,MPI_DOUBLE,0,basetag*mc_tag_data+myid,MPI_COMM_WORLD);
  
  fprintfDEBUG(stderr,"%d Sending X_ded\n",myid);
  MPI_Send(psim->X_ded,buffsize*psim->n_ded,MPI_DOUBLE,0,basetag*mc_tag_data+myid,MPI_COMM_WORLD);
  
  fprintfDEBUG(stderr,"%d Sending log_rho\n",myid);
  MPI_Send(psim->log_rho,buffsize,MPI_DOUBLE,0,basetag*mc_tag_logrho+myid,MPI_COMM_WORLD);
  
  fprintfDEBUG(stderr,"%d Sending flags\n",myid);
  MPI_Send(psim->flg,buffsize,MPI_SHORT,0,basetag*mc_tag_flag+myid,MPI_COMM_WORLD);
  
  fprintfDEBUG(stderr,"%d maxW \n",myid);
  MPI_Send(&psim->maxW, 1, MPI_DOUBLE, 0,basetag*mc_tag_maxw+myid,MPI_COMM_WORLD);
  
  fprintfDEBUG(stderr,"%d sending maxR=%g &maxR=%p\n",myid,psim->maxR,(void*)(&psim->maxR));
  MPI_Send(&psim->maxR, 1, MPI_DOUBLE, 0,basetag*mc_tag_maxr+myid,MPI_COMM_WORLD);
  
  fprintfDEBUG(stderr,"%d end send\n",myid);
  
  return;
}

/* Master receives from clients */
size_t receive_importance_weight(pmc_simu *psim,int nproc,int master_cnt,int master_sample,error **err)
{
  size_t current_sample;
  size_t i,j,cnt,inok,buffsize;
  int tag;
  double max;
  
  cnt=master_cnt;  
  current_sample=master_sample;
  tag = nproc;
  for (i=1;i<nproc;i++) { 
    
    fprintfDEBUG(stderr,"I am waiting for inok from %zu on %zu\n",i,tag*mc_tag_nok+i);
    MPI_Recv(&inok,1,MPI_LONG,i,tag*mc_tag_nok+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    fprintfDEBUG(stderr,"master(%zu):: inok=%zu\n",i,inok);
    
    fprintfDEBUG(stderr,"I am the master, I am receiving buffsize from %zu on %zu\n",i,tag*mc_tag_data_size+i);
    MPI_Recv(&buffsize,1,MPI_LONG,i,tag*mc_tag_data_size+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    fprintfDEBUG(stderr,"master(%zu):: buffsize=%zu\n",i,buffsize);
    
    cnt+=inok;
    if (inok==0) {
      for(j=0;j<buffsize;j++) {
        psim->weights[current_sample+j]=0;
        psim->log_rho[current_sample+j]=0;
        psim->flg[current_sample+j]=0;
      }
      current_sample+=buffsize;
      fprintfDEBUG(stderr,"node %d is dead ?\n", (int)i);
      continue;
    }
    
    testErrorRetVA(current_sample+buffsize>psim->nsamples,pmc_badComm,
                   "Not the right number of elements (Expected %d got %d)",
                   *err,__LINE__,0,psim->nsamples-current_sample,buffsize);
    
    fprintfDEBUG(stderr,"I am waiting for weights from %zu on %zu\n",i,tag*mc_tag_data+i);
    MPI_Recv(&(psim->weights[current_sample]),buffsize,
             MPI_DOUBLE,i,tag*mc_tag_data+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    fprintfDEBUG(stderr,"I am waiting for X_ded from %zu on %zu\n",i,tag*mc_tag_data+i);
    MPI_Recv(&(psim->X_ded[current_sample*psim->n_ded]),buffsize*psim->n_ded,
             MPI_DOUBLE,i,tag*mc_tag_data+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    fprintfDEBUG(stderr,"I am waiting for log_rho from %zu on %zu\n",i,tag*mc_tag_data+i);
    MPI_Recv(&(psim->log_rho[current_sample]),buffsize,
             MPI_DOUBLE,i,tag*mc_tag_logrho+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    fprintfDEBUG(stderr,"I am waiting for flg from %zu on %zu\n",i,tag*mc_tag_data+i);
    MPI_Recv(&(psim->flg[current_sample]),buffsize,
             MPI_DOUBLE,i,tag*mc_tag_flag+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    fprintfDEBUG(stderr,"I am waiting for max from %zu on %zu\n",i,tag*mc_tag_flag+i);
    MPI_Recv(&max,1,MPI_DOUBLE,i,tag*mc_tag_maxw+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    if (max>psim->maxW)
      psim->maxW=max;
    
    fprintfDEBUG(stderr, "received maxW=%g\n", psim->maxW);
    MPI_Recv(&max,1,MPI_DOUBLE,i,tag*mc_tag_maxr+i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    if (max>psim->maxR)
      psim->maxR=max;
    current_sample+=buffsize;
  }
  
  return cnt;
}


/******************************************************************************/
/******************** exchanging mvdens proposal ******************************/
/******************************************************************************/

// This routine is called by master to send update MOG to clients
void send_mix_mvdens(mix_mvdens *m, int nproc,error **err) {
  int tag,i;
  size_t buffsz;
  void* buff;
  
  tag=nproc;
  //mix_mvdens_print(stderr,m);
  buff = serialize_mix_mvdens(m,&buffsz,err);
  forwardError(*err,__LINE__,);
  for (i=1;i<nproc;i++) {
    fprintfDEBUG(stderr,"I am the master, I am sending buffsize to %d on %d (%zu)\n",i,tag*mc_tag_mix_mvdens_size+i,buffsz);
    MPI_Send(&buffsz,1,MPI_LONG,i,tag*mc_tag_mix_mvdens_size+i,MPI_COMM_WORLD);
    
    fprintfDEBUG(stderr,"I am the master, I am sending data to %d on %d\n",i,tag*mc_tag_mix_mvdens_buffer+i);
    MPI_Send(buff,buffsz,MPI_CHAR,i,tag*mc_tag_mix_mvdens_buffer+i,MPI_COMM_WORLD);
  }
}

/* Clients receive from master */
mix_mvdens * receive_mix_mvdens(int myid, int nproc,error **err) {
  int tag;
  size_t buffsz;
  void* buff;
  mix_mvdens *self;
  
  tag = nproc;

  fprintfDEBUG(stderr,"%d:: I am waiting for guess on %d\n",myid,tag*mc_tag_mix_mvdens_size+myid);
  MPI_Recv(&buffsz,1,MPI_LONG,0,tag*mc_tag_mix_mvdens_size+myid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  fprintfDEBUG(stderr,"%d:: buffsz=%zu\n",myid,buffsz);
  buff = malloc_err(buffsz,err);
  forwardError(*err,__LINE__,NULL);

  fprintfDEBUG(stderr,"%d:: I am waiting for guess on %d\n",myid,tag*mc_tag_mix_mvdens_buffer+myid);
  MPI_Recv(buff,buffsz,MPI_CHAR,0,tag*mc_tag_mix_mvdens_buffer+myid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  fprintfDEBUG(stderr,"%d:: received guess on %d\n",myid,tag*mc_tag_mix_mvdens_buffer+myid);
  
  self=deserialize_mix_mvdens(buff,buffsz,err);
  forwardError(*err,__LINE__,NULL);
  //mix_mvdens_print(stderr,self);
  free(buff);
  return self;
}


