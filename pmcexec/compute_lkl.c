#include "readConf.h"
#include "pmc_rc.h"

#define __FUNC__ "main"
int main(int argc, char **argv) {
  confFile* rc;
  error *_err;
  error **err;
  distribution *dist;
  double *test_pars,res,dres;
  size_t n_test_pars;
  int flg;
  double rel;
  
  // init error management
  _err = initError();
  err = &_err;
  
  printf("----- starting %s -----\n",argv[0]);
  
  // read the parameter file
  rc = rc_init_from_args(argc,argv,err);
  quitOnError(*err,__LINE__,stderr);
  
  // initialize the distribution
  dist = init_distribution_from_rc(rc,"testbed.target",err);
  quitOnError(*err,__LINE__,stderr);
  
  // read testpars et result
  n_test_pars = rc_get_real_array(rc,"testbed.pars",&test_pars,err);
  quitOnError(*err,__LINE__,stderr);
  
  testErrorExitVA(n_test_pars!=dist->ndim,-10,"not the right number of dims ! (got %d expected %d)",*err,__LINE__,n_test_pars,dist->ndim);
  
  // compute the log likelihood
  dres = distribution_lkl(dist,test_pars,err);
  quitOnError(*err,__LINE__,stderr);
  
  printf("log lkl = %g\n",dres);
  
  // if result is given in the parfile, compare the result
  flg = rc_has_key(rc,"testbed.result",err);
  quitOnError(*err,__LINE__,stderr);
  
  if (flg==1) {
    res = rc_get_real(rc,"testbed.result",err);
    quitOnError(*err,__LINE__,stderr);
    if (res==0) {
      rel = res-dres;
    } else {
      rel = (res-dres)/res;
    }
    printf("expected log lkl = %g -- diff = %g relative error = %g\n",res,res-dres,rel);
  }

  // cleanup
  free_distribution(&dist);
  rc_close(&rc);
  endError(err);
  
  printf("----- end %s -----\n",argv[0]);
  return 0; 
}