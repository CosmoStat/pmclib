/*
 *  gridMe.h
 *  likely
 *
 *  Created by Karim Benabed on 08/07/08.
 *  Copyright 2008 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifndef __GRID_H
#define __GRID_H
#ifdef __PLANCK__

#include "HL2_likely/tools/errorlist.h"
#include "HL2_likely/pmclib/parabox.h"
#include "HL2_likely/pmclib/distribution.h"
#include "HL2_likely/pmclib/pmc_mpi.h"

#else

#include "errorlist.h"
#include "parabox.h"
#include "distribution.h"
#include "pmc_mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct {
  int ndim;
  long totbins;
  void *data;
  double *limits;
  double *bins;
  double *res;
  short *flag;
  int *nbin;
  int isLog;
  double **pos;
  
  parabox *pb;
  distribution *target;
  
} grid_simu;

// high level
grid_simu *init_grid(distribution *target, parabox* pb, int *nbins, double *mM, int isLog, double** pos,error **err);
void doGrid(grid_simu* grid, error **err);
void free_grid_simu(grid_simu **grid);
void grid_simu_dump(FILE *where,grid_simu *grid,error **err);

//low_level
grid_simu *init_grid_simu(int ndim,double *mM,int *nbins,error **err);
void loc_grid(grid_simu* grid,int* locint, double *pars,error **err);
void int2loc(grid_simu *grid,long pos,int* mloc,error **err);
void gridMe(grid_simu *grid, void *target, posterior_log_pdf_func *target_func, parabox *pb,error **err);
void setLog_grid(grid_simu* grid,error **err);
void compute_pos(grid_simu* grid,error **err) ;
grid_simu *alloc_grid_simu(int ndim,int *nbins,error **err);

// hdf
#ifdef _WITH_HDF5_
#include "hdf5.h"
void grid_simu_hdfdump(grid_simu *grid, char* fname, error **err);
#define hdf5_base  -1000
#endif
// rc
#ifdef _WITH_RC_
#include "readConf.h"
#include "pmc_rc.h"
grid_simu *init_grid_from_rc(confFile *rc, char * root, error **err);
#endif
#endif
