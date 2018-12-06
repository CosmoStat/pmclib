#include "sampleTarget.h"

#define __FUNC__ "init_bentGauss"
bentGauss * init_bentGauss(int ndim, double b, double sig1sq, error **err) {
  bentGauss *bg;
  
  testErrorRetVA(ndim<2,-1000,"Number of dimensions too small (got %d, expected at least 2)",*err,__LINE__,NULL,ndim);
  
  bg = malloc_err(sizeof(bentGauss),err);
  forwardError(*err,__LINE__,NULL);
  
  bg->ndim = ndim;
  bg->b = b;
  bg->sig1 = sig1sq;
  
  return bg;
}
#undef __FUNC__

#define __FUNC__ "init_bentGauss_distribution"
distribution* init_bentGauss_distribution(int ndim, double b, double sig1, error **err) {
  distribution *dst;
  bentGauss *bg;
  char **chname;
  int i;
  
  bg = init_bentGauss(ndim,b,sig1,err);
  forwardError(*err,__LINE__,NULL);
  
  dst = init_distribution(ndim, bg, bentGauss_lkl, bentGauss_free, NULL, err);
  forwardError(*err,__LINE__,NULL);
  
  chname = malloc_err(sizeof(char*)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  chname[0] = malloc_err(sizeof(char)*100,err);
  sprintf(chname[0],"bent");

  chname[1] = malloc_err(sizeof(char)*100,err);
  sprintf(chname[1],"sig1_square");

  for(i=2;i<ndim;i++) {
    chname[i] = malloc_err(sizeof(char)*100,err);
    sprintf(chname[i],"extra_%d",i-1);
    
  }
  
  distribution_set_names(dst,chname,err);
  forwardError(*err,__LINE__,NULL);
  
  for(i=0;i<ndim;i++) {
    free(chname[i]);
  }
  free(chname);
  
  return dst;  
}
#undef __FUNC__

#define __FUNC__ "rcinit_bentGauss"
void* rcinit_bentGauss(confFile* rc, char* root, error **err) {
  long ndim;
  double b,sig1;
  distribution *dst;
  confFile *rca;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  ndim = rc_get_integer(rca,"ndim",err);
  forwardError(*err,__LINE__,NULL);
  
  b = rc_get_real(rca,"bent",err);
  forwardError(*err,__LINE__,NULL);

  sig1 = rc_get_real(rca,"sig1_square",err);
  forwardError(*err,__LINE__,NULL);
    
  dst =  init_bentGauss_distribution(ndim,b,sig1,err);
  forwardError(*err,__LINE__,NULL);
  rc_close(&rca);

  return (void*) dst;
}
#undef __FUNC__


#define __FUNC__ "bentGauss_lkl"
double bentGauss_lkl(void* vbg, double* pars, error **err) {
  double res,x2;
  int i;
  bentGauss *bg;
  
  bg = vbg;
  /*for(i=0;i<bg->ndim;i++) {
    _DEBUGHERE_("%d->%g",i,pars[i]);
  }*/
  res = 0;
  x2 = pars[0]*pars[0];
  res += x2/bg->sig1;
  x2 -= bg->sig1;
  x2 *= bg->b;
  x2 += pars[1];
  res += x2*x2;
  for(i=2;i<bg->ndim;i++) {
    res += pars[i]*pars[i];
  }
  //_DEBUGHERE_("r %g",res);
  return -res/2.;
}
#undef __FUNC__


#define __FUNC__ "bentGauss_free"
void bentGauss_free(void** pbg) {
  free(*pbg);
  *pbg = NULL;
}
#undef __FUNC__

