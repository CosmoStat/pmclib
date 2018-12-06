/*
 *  test_grid.c
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 30/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


// some includes
#include "gridMe.h"

#include "readConf.h"
#include "pmc_rc.h"


#define __FUNC__ "main"
int main(int argc, char** argv) {
  error *_err, **err;   // error management

  distribution *target; // And the target there
  
  parabox *pb;          // This is a parameter box
  grid_simu *grid;      // the grid object
  
  confFile *rc;   // parameter file
  
  int ndim;
  int mpi_rank;
  
  char* gridFile;
  FILE* f;
  
  /**** inits  ****************************************************************/ 
  // the main has to init mpi
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  // initialize error management 
  _err = initError();
  err = &_err;
  
  if(mpi_rank == 0) {
    printf("----- starting %s -----\n",argv[0]);
  }

  // read the parameter file
  rc = rc_init_from_args(argc,argv,err);
  quitOnError(*err,__LINE__,stderr);

  
  // init grid from parameter file
  grid = init_grid_from_rc(rc,"grid",err);
  quitOnError(*err,__LINE__,stderr);  
  
  // run loop
  doGrid(grid,err);
  quitOnError(*err,__LINE__,stderr);
  
  // save result
  if (mpi_rank==0) {
    gridFile = rc_get_string(rc,"grid.file",err);
    quitOnError(*err,__LINE__,stderr);
    
    fprintf(stdout,"Saving in '%s'\n",gridFile);
    #ifdef _WITH_HDF5_
    grid_simu_hdfdump(grid, gridFile, err);
    quitOnError(*err,__LINE__,stderr);
    #else
    f=fopen_err(gridFile,"w",err);
    quitOnError(*err,__LINE__,stderr);
    grid_simu_dump(f, grid, err);
    quitOnError(*err,__LINE__,stderr);
    fclose(f);
    #endif
  }
  
  //cleanup
  free_grid_simu(&grid);
  rc_close(&rc);
  endError(err);

  if(mpi_rank == 0) {
    printf("----- end %s -----\n",argv[0]);
  }
  
  MPI_Finalize();
  
  exit(0);
}
  