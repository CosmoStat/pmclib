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
#include "optimize.h"
#include "readConf.h"
#include "pmc_rc.h"


#define __FUNC__ "main"
int main(int argc, char** argv) {
  error *_err, **err;   // error management
  confFile *rc;   // parameter file
  optimize_struct *opt;
  double *pars,*var;
  int i;
  char *fname;
  
  /**** inits  ****************************************************************/ 
  // initialize error management 
  
  _err = initError();
  err = &_err;
    
  printf("----- starting %s -----\n",argv[0]);
  
  // read the parameter file
  rc = rc_init_from_args(argc,argv,err);
  quitOnError(*err,__LINE__,stderr);
  
  // init optimize
  opt = init_optimize_from_rc(rc,"optimize",err);
  quitOnError(*err,__LINE__,stderr);  

  optimize(opt,err);
  quitOnError(*err,__LINE__,stderr);  
  
  pars = optimize_get_best_pars(opt,err);
  quitOnError(*err,__LINE__,stderr);  
  
  fprintf(stdout,"Best parameters : \n");
  for(i=0;i<opt->ndim;i++) {
    fprintf(stdout,"%g ",pars[i]);
  }
  fprintf(stdout,"\n");
  
  fname = rc_safeget_string(rc,"parfile","",err);
  quitOnError(*err,__LINE__,stderr);  
  
  if (fname[0]!='\0') {
    FILE* f;
    f=fopen_err(fname,"w",err);
    quitOnError(*err,__LINE__,stderr);  
    for(i=0;i<opt->ndim;i++) {
      fprintf(f,"%40g ",pars[i]);
    }
    fprintf(f,"\n");
    fclose(f);
  }

  free(pars);
  
  fname = rc_safeget_string(rc,"varfile","",err);
  quitOnError(*err,__LINE__,stderr);  
  
  if (fname[0]!='\0') {
    FILE* f;
    
    var = optimize_get_variance(opt,err);
    quitOnError(*err,__LINE__,stderr);  
    
    f=fopen_err(fname,"w",err);
    quitOnError(*err,__LINE__,stderr);  
    for(i=0;i<opt->ndim*opt->ndim;i++) {
      fprintf(f,"%40g ",var[i]);
    }
    fprintf(f,"\n");
    fclose(f);
    free(var);
  }
  
  
  //cleanup
  optimize_free(&opt);
  rc_close(&rc);
  endError(err);
  printf("----- end %s -----\n",argv[0]);
    
  exit(0);
}
