/*
 *  test_pmc.c
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
  pmc_simu *psim;       // the pmc object
  
  confFile *rc;   // parameter file
  
  gsl_rng *r;
  
  FILE *f,*pfile;
  char *pmcFile;
  double evi,ess,perp;
  int niter,iter,flg;
  
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

  // init importance run from parameter file
  psim = init_pmc_from_rc(rc,"pmc",err);
  quitOnError(*err,__LINE__,stderr);  

  // get number of iteration
  niter = rc_get_integer(rc,"pmc.niter",err);
  quitOnError(*err,__LINE__,stderr);  
  
  flg = rc_has_key(rc,"pmc.perpfile",err);
  quitOnError(*err,__LINE__,stderr);  
  if (flg==1) {
    char *cpfile;
    cpfile = rc_get_string(rc,"pmc.perpfile",err);
    quitOnError(*err,__LINE__,stderr);  

    pfile = fopen_err(cpfile,"w",err);
    quitOnError(*err,__LINE__,stderr);  
  } else {
    pfile = NULL;
  }

  //run loop
  for (iter=0;iter<niter;iter++) {
    char *rname;
    fprintf(stdout,"Iteration %d\n",iter);
    mix_mvdens_print(stderr,psim->proposal->data);

    perp = pmc_simu_pmc_step(psim,r,err);
    quitOnError(*err,__LINE__,stderr);
    
    fprintf(stdout,"Iteration %d : save to disk\n",iter);
      
    // save sample to disk 
    rname = get_itername(rc,"pmc.resfile",iter,err);
    quitOnError(*err,__LINE__,stderr);
    
    fprintf(stdout,"data -> %s\n",rname);
#ifdef _WITH_HDF5_
    pmc_simu_hdfdump(psim,rname, err);
    quitOnError(*err,__LINE__,stderr);
#else
    f=fopen_err(rname,"w",err);
    quitOnError(*err,__LINE__,stderr);
    pmc_simu_dump(f, psim, err);
    quitOnError(*err,__LINE__,stderr);
    fclose(f);
#endif

    rname = get_itername(rc,"pmc.propfile",iter,err);
    quitOnError(*err,__LINE__,stderr);

    fprintf(stdout,"updated proposal -> %s\n",rname);
#ifdef _WITH_HDF5_
    mix_mvdens_hdfdump(psim->proposal->data,rname,err);
    quitOnError(*err,__LINE__,stderr);
#else
    f=fopen_err(rname,"w",err);
    quitOnError(*err,__LINE__,stderr);
    mix_mvdens_dump(f,psim->proposal->data);
    fclose(f);
#endif
    // printout proposal
    mix_mvdens_print(stderr,psim->proposal->data);
    
    perp = perplexity_and_ess(psim, MC_NORM, &ess,err);
    quitOnError(*err,__LINE__,stderr);

    evi = evidence(psim, NULL,err);
    quitOnError(*err,__LINE__,stderr);

    printf("perplexity %g effectiveSampleSize %g evidence %g\n",perp,ess,evi);
    if (pfile != NULL) {
      fprintf(pfile,"%g %g %g\n",perp,ess,evi);      
    }
    
  }
  if (pfile!=NULL) {
    fclose(pfile);
  }
  
  //cleanup
  pmc_simu_free(&psim);
  gsl_rng_free(r);
  rc_close(&rc);
  endError(err);
  printf("----- end %s -----\n",argv[0]);
    
  exit(0);
}
