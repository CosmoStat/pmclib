/*
 *  test_importance_mpi.c
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


// some includes
#include "pmc_mpi.h"

#include "readConf.h"


#define __FUNC__ "main"
int main(int argc, char** argv) {
  error *_err, **err;   // error management

  parabox *pb;          // This is a parameter box
  pmc_simu *psim;       // the pmc object
  
  confFile *rc;   // parameter file
  
  gsl_rng *r;
  
  FILE* f;
  char *pmcFile;
  double evi,ess,perp;
  
  int mpi_rank;
  
  /**** inits  ****************************************************************/ 
  // the main has to init mpi
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  // initialize error management 
  _err = initError();
  err = &_err;
  
  // and nuke gsl_error_handler
  gsl_set_error_handler_off();
  
  if (mpi_rank==0) {
    printf("----- starting %s -----\n",argv[0]);
  }
  
  // read the parameter file
  rc = rc_init_from_args(argc,argv,err);
  quitOnError(*err,__LINE__,stderr);

  //init gsl_rng from parameter file (only master)
  if (mpi_rank==0) {
    r = rc_gsl_rng(rc,"rng",err);
    quitOnError(*err,__LINE__,stderr);    
  }
  
  // init importance run from parameter file
  psim = init_importance_mpi_from_rc(rc,"importance",err);
  quitOnError(*err,__LINE__,stderr);  
  
  //run loop
  pmc_simu_importance_mpi(psim,r,err);
  quitOnError(*err,__LINE__,stderr);  
  
  // save result (only master)
  
  if (mpi_rank==0) {
      pmcFile = rc_get_string(rc,"importance.file",err);
      
      quitOnError(*err,__LINE__,stderr);

      fprintf(stdout,"Saving in '%s'\n",pmcFile);
    #ifdef _WITH_HDF5_
      pmc_simu_hdfdump(psim,pmcFile, err);
      quitOnError(*err,__LINE__,stderr);
    #else
      f=fopen_err(pmcFile,"w",err);
      quitOnError(*err,__LINE__,stderr);
      pmc_simu_dump(f, psim, err);
      quitOnError(*err,__LINE__,stderr);
      fclose(f);
    #endif    
  
    perp = perplexity_and_ess(psim, MC_NORM, &ess,err);
    quitOnError(*err,__LINE__,stderr);
  
    evi = evidence(psim, NULL,err);
    quitOnError(*err,__LINE__,stderr);
  
    printf("perplexity %g effectiveSampleSize %g evidence %g\n",perp,ess,evi);
  }
  
  //cleanup
  pmc_simu_free(&psim);
  if (mpi_rank==0) {
    gsl_rng_free(r);    
  }
  rc_close(&rc);
  endError(err);

  if (mpi_rank==0) {
    printf("----- end %s -----\n",argv[0]);
  }
  
  MPI_Finalize();
  
  exit(0);
}
