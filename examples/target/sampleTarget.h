#ifndef __SMP_TARGET__
#define __SMP_TARGET__

#include "distribution.h"
#include "readConf.h"

typedef struct { 
  int ndim;
  double b, sig1;
  } bentGauss;

  
bentGauss * init_bentGauss(int ndim, double b, double sig1, error **err);
distribution *init_bentGauss_distribution(int ndim, double b, double sig1, error **err);
double bentGauss_lkl(void* bg, double* pars, error **err);
void bentGauss_free(void** bg);

void* rcinit_bentGauss(confFile* rc, char* root,error **err);

#endif
