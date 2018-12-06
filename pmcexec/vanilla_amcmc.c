/*
 *  test_mcmc.c
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


// some includes
#include "pmc.h"

#include "readConf.h"
#include "pmc_rc.h"


#define __FUNC__ "main"
int main(int argc, char** argv) {
  error *_err, **err;   // error management

  parabox *pb;          // This is a parameter box
  mcmc_run *mcr;        // the mcmc object
  
  confFile *rc;   // parameter file
  
  gsl_rng *r;
  
  FILE *f,*pfile;
  char *mcmcFile;
  double evi,ess,perp;
  int niter,iter,flg;
  int more;
  
  /**** inits  ****************************************************************/ 
  // initialize error management 
  _err = initError();
  err = &_err;
  
  // and nuke gsl_error_handler
  gsl_set_error_handler_off();
  
  printf("----- starting %s -----\n",argv[0]);
  
  // read the parameter file
  rc = rc_init_from_args(argc,argv,err);
  quitOnError(*err,__LINE__,stderr);

  //init gsl_rng from parameter file
  r = rc_gsl_rng(rc,"rng",err);
  quitOnError(*err,__LINE__,stderr);

  // init mcmc from parfile
  mcr = init_mcmc_from_rc(rc,"mcmc",err);
  quitOnError(*err,__LINE__,stderr);
  
  // read name of file to save in
  mcmcFile = rc_get_string(rc,"mcmc.file",err);
  quitOnError(*err,__LINE__,stderr);
  
  for(more=1;more==1;) {
    // run batches of mcmc until we have enough points
    more = mcmc_run_batch(mcr, r, err);
    quitOnError(*err,__LINE__,stderr);
    
    printf("simulated %d points (acceptance ration %g%%)\n",mcr->current,mcr->naccepted*100./mcr->current);
#ifdef _WITH_HDF5_    
    // save in hdf5
    mcmc_simu_hdfdump(mcr,mcmcFile,err);
    quitOnError(*err,__LINE__,stderr); 
    
# else
    fprintf(stderr,"IMPLEMENT YOUR OWN MCMC DUMP HERE\n");
#endif

  }
  
  //cleanup
  free_mcmc_run(&mcr);
  gsl_rng_free(r);
  rc_close(&rc);
  endError(err);

  printf("----- end %s -----\n",argv[0]);

  exit(0);
}
  
    
