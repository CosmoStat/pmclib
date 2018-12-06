/*
 *  gridMe.c
 *  likely
 *
 *  Created by Karim Benabed on 08/07/08.
 *  Copyright 2008 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */
#ifdef __PLANCK__
#include "HL2_likely/pmclib/gridMe.h"
#else
#include "gridMe.h"
#endif

// high level interface

grid_simu *init_grid(distribution *target, parabox* pb, int *nbins, double *mM,  int isLog, double** pos,error **err) {
  int ndim;
  grid_simu *grid;
  
  ndim = target->ndim;
  if (pb!=NULL) {
    testErrorRetVA(ndim != pb->ndim,-10,"target and parabox have different dimensions ! (found %d and %d)",*err,__LINE__,NULL,ndim,pb->ndim);
  }
  
  if (mM!=NULL) {
    grid = init_grid_simu(ndim, mM, nbins, err);
    forwardError(*err,__LINE__,NULL);

    if (isLog) {
      setLog_grid(grid,err);
      forwardError(*err,__LINE__,NULL);
    }
  } else {
    int i;
    grid = alloc_grid_simu(ndim,nbins,err);
    grid->pos = malloc_err(sizeof(double*)*grid->ndim,err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<grid->ndim;i++) {
      grid->pos[i] = malloc_err(sizeof(double)*grid->nbin[i],err);
      forwardError(*err,__LINE__,NULL);
      memcpy(grid->pos[i],pos[i],sizeof(double)*grid->nbin[i]);
    }
  }
  
  grid->target = target;
  grid-> pb = pb;
  
  return grid;
}

void doGrid(grid_simu *grid,error **err) {
  gridMe(grid, grid->target, distribution_lkl, grid->pb, err);
  forwardError(*err,__LINE__,);
}

void free_grid_simu(grid_simu **grid) {
  int i;
  
  if((*grid)->data!=NULL) {
    free((*grid)->data);
  }
  free((*grid)->nbin);
  free((*grid)->res);
  free((*grid)->flag);
  if ((*grid)->pos !=NULL) {
    for(i=0;i<(*grid)->ndim;i++) {
      free((*grid)->pos[i]);
    }
    free((*grid)->pos);
  }
  if ((*grid)->pb != NULL) {
    free_parabox(&((*grid)->pb));
  }
  if ((*grid)->target != NULL) {
    free_distribution(&((*grid)->target));
  }
  free(*grid);
  *grid = NULL;
}

// low level interface & io
grid_simu *alloc_grid_simu(int ndim,int *nbins,error **err) {
  grid_simu *grid;
  int i;
  
  grid = malloc_err(sizeof(grid_simu),err);
  forwardError(*err,__LINE__,NULL);
  
  grid->ndim=ndim;
  
  grid->nbin = malloc_err(sizeof(int)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  grid->totbins = 1;
  
  for(i=0;i<ndim;i++) {
    grid->nbin[i]=nbins[i];
    grid->totbins *= nbins[i];
  }
  
  grid->res = malloc_err(sizeof(double)*grid->totbins,err);
  forwardError(*err,__LINE__,NULL);
  
  grid->flag = malloc_err(sizeof(short)*grid->totbins,err);
  forwardError(*err,__LINE__,NULL);
  
  grid->isLog = 0;
  
  grid->target = NULL;
  grid->pb = NULL;
  grid->data = NULL;
  grid->limits = NULL;
  grid->bins = NULL;
  grid->pos = NULL;
  return grid;
  
}

grid_simu *init_grid_simu(int ndim,double *mM,int *nbins,error **err) {
  grid_simu *grid;
  int i;
  
  grid = alloc_grid_simu(ndim,nbins,err);
  forwardError(*err,__LINE__,NULL);
  
  grid->data = malloc_err(sizeof(double)*ndim*3,err);
  forwardError(*err,__LINE__,NULL);
  
  grid->limits = grid->data;
  grid->bins = grid->limits + 2*ndim;
  
  for(i=0;i<ndim;i++) {
    grid->limits[i*2] = mM[i*2];
    grid->limits[i*2+1] = mM[i*2+1];
    grid->bins[i] = (mM[i*2+1]-mM[i*2])/grid->nbin[i];
  }
  
  grid->isLog = 0;
  compute_pos(grid,err);
  forwardError(*err,__LINE__,NULL);

  return grid;
}

void setLog_grid(grid_simu* grid,error **err) {
  int i;
  
  testErrorRet(grid->limits==NULL,pmc_undef,"NO !",*err,__LINE__,);
  
  grid->isLog=1;
  for(i=0;i<grid->ndim;i++) {
    grid->bins[i] = (log(grid->limits[i*2+1])-log(grid->limits[i*2]))/grid->nbin[i];
    //fprintf(stderr,"ppppp %g %g %d %g\n",log(grid->limits[i*2+1]),log(grid->limits[i*2]), grid->bins[i]);
  }
  compute_pos(grid,err);
  forwardError(*err,__LINE__,);
  
  return ;
}

void compute_pos(grid_simu* grid,error **err) {
  double **pos;
  int i,j;
  double *pars;
  int *locint;
  
  if (grid->pos !=NULL) {
    for(i=0;i<grid->ndim;i++) {
      free(grid->pos[i]);
    }
    free(grid->pos);
    grid->pos = NULL;
  }
  
  pars = malloc_err(sizeof(double)*grid->ndim,err);
  forwardError(*err,__LINE__,);
  locint = malloc_err(sizeof(int)*grid->ndim,err);
  forwardError(*err,__LINE__,);
  memset(locint,0,sizeof(int)*grid->ndim);
  
  pos = malloc_err(sizeof(double*)*grid->ndim,err);
  forwardError(*err,__LINE__,);
  for(i=0;i<grid->ndim;i++) { 
    pos[i] = malloc_err(sizeof(double)*grid->nbin[i],err);
    forwardError(*err,__LINE__,);
    for(j=0;j<grid->nbin[i];j++) {
      locint[i]=j;
      loc_grid(grid,locint, pars,err);
      forwardError(*err,__LINE__,);
      pos[i][j] = pars[i];
    }
    locint[i]=0;
  }
  grid->pos = pos;
  free(pars);
  free(locint);
}

void loc_grid(grid_simu* grid,int* locint, double *pars,error **err) {
  int i;

  if (grid->pos!=NULL) {
    for(i=0;i<grid->ndim;i++) {
      testErrorRet((locint[i]>=grid->nbin[i]) || (locint[i]<0),pmc_outOfBound,"Indice out of bounds",*err,__LINE__,);
      pars[i] = grid->pos[i][locint[i]];
    }
    return;
  }
  for(i=0;i<grid->ndim;i++) {
    //fprintf(stderr,"%d %d \n",locint[i],grid->nbin[i]);
    testErrorRet((locint[i]>=grid->nbin[i]) || (locint[i]<0),pmc_outOfBound,"Indice out of bounds",*err,__LINE__,);
    if (grid->isLog==1) {
      pars[i] = grid->limits[i*2]*exp( (locint[i]+0.5)*grid->bins[i]);
      //fprintf(stderr,"-> %g %d %g %g\n",grid->limits[i*2],locint[i],grid->bins[i],pars[i]);
    } else {
      pars[i] = grid->limits[i*2] + grid->bins[i] * (locint[i] + 0.5);
    }
  }
  
  return ;
}

void int2loc(grid_simu *grid,long pos,int* mloc,error **err) {
  int i;
  long cur;
  
  testErrorRet((pos>=grid->totbins) || (pos<0),pmc_outOfBound,"Indice out of bounds",*err,__LINE__,);
  
  cur=pos;
  for(i=0;i<grid->ndim;i++) {
    mloc[grid->ndim-1-i] = cur % grid->nbin[grid->ndim-1-i];
    cur /= grid->nbin[grid->ndim-1-i];
  }
}

void gridMe(grid_simu *grid, void *target, posterior_log_pdf_func *target_func, parabox *pb,error **err) {
  long pos,rmn,step,begin,end,iproc;
  int *loc;
  double *pars;
  double res;
  int size,rank;
  int ppc;
  
  loc = malloc_err(sizeof(int)*grid->ndim,err);
  forwardError(*err,__LINE__,);
  
  pars = malloc_err(sizeof(double)*grid->ndim,err);
  forwardError(*err,__LINE__,);
  
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  //fprintf(stderr, ">> %d 1\n",rank);
  
  step = grid->totbins / size;
  rmn = grid->totbins % size;
  if (rank==0) {
    begin = 0;
    end = step;
  } else if (rank <= rmn) {
    begin = step * rank + rank - 1;
    end = begin + step + 1;
  } else {
    begin = step * rank + rmn;
    end = begin + step;
  }
  
  //fprintf(stderr,"%d:: %d -> %d\n\n\n",rank,begin,end);
  
  //fprintf(stderr, ">> %d 2 (%d %d %d)\n",rank,begin,end,end-begin);
  ppc=0;
  
  for(pos=begin;pos<end;pos++) {
    int pc;
    pc=((pos-begin)*1./(end-begin)*100.)/5.;
    
    if (pc>ppc) {
      ppc=pc;
      //fprintf(stderr,"%d: %d done\n",rank,ppc*5);
    }
    //fprintf(stderr, ">> %d %d %d %d\n",rank,pos,begin,end);
    // compute position
    int2loc(grid, pos, loc, err);
    forwardError(*err,__LINE__,);
    loc_grid(grid,loc,pars,err);
    forwardError(*err,__LINE__,);
    //fprintf(stderr,">> %d (%d %d) (%g %g)\n",rank,loc[0],loc[1],pars[0],pars[1]);
    // compute res
    grid->flag[pos] = 0;
    grid->res[pos]=0;
    if (pb!=NULL) {
      if (isinBox(pb,pars,err)==0) {
        forwardError(*err,__LINE__,);
        continue;
      }
      forwardError(*err,__LINE__,);
    }
    //fprintf(stderr, ">> %d before\n",rank);    
    res = target_func(target,pars,err);
    //fprintf(stderr,"%d :: %g %g -> %g\n",rank,pars[0],pars[1],res);
    /*if((pars[0]==1) && (pars[1]==0))
     fprintf(stderr,"%g %g : %g (%d)\n",pars[0],pars[1],res,pos);*/
    forwardErrorNoReturn(*err,__LINE__);
    if (isError(*err)) {
      // if error, print, let the flag to 0 and continue
      printError(stderr,*err);
      purgeError(err);
      //fprintf(stderr,">> %d ERROR (%d %d) (%g %g)\n",rank,loc[0],loc[1],pars[0],pars[1]);
      continue;
    }
    grid->flag[pos]=1;
    //fprintf(stderr, ">> %d after\n",rank);    
    grid->res[pos]=res;
    
  }
  
  if (rank == 0) {
    pos=end;
    for(iproc=1;iproc<size;iproc++) {
      //fprintf(stderr,">> 0 rcv to %d\n",iproc);
      MPI_Recv(&step,1,MPI_LONG,iproc,size*mc_tag_simulation_size+iproc,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      //fprintf(stderr,">> 0 rcv flag to %d\n",iproc);
      MPI_Recv(grid->flag+pos,step,MPI_SHORT,iproc,size*mc_tag_flag+iproc,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      //fprintf(stderr,">> 0 rcv data to %d\n",iproc);
      MPI_Recv(grid->res+pos,step,MPI_DOUBLE,iproc,size*mc_tag_data+iproc,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      //fprintf(stderr,"rcv %d %g %g %g %g %g \n",iproc,grid->res[pos+0],grid->res[pos+1],grid->res[pos+2],grid->res[pos+3],grid->res[pos+4]);
      pos+=step;
    }
  } else {
    step = end-begin;
    //fprintf(stderr,">> %d send\n",rank);
    MPI_Send(&step,1,MPI_LONG,0,size*mc_tag_simulation_size+rank,MPI_COMM_WORLD);
    //fprintf(stderr,">> %d send flag\n",rank);
    MPI_Send(grid->flag+begin,step,MPI_SHORT,0,size*mc_tag_flag+rank,MPI_COMM_WORLD);
    //fprintf(stderr,">> %d send data\n",rank);
    MPI_Send(grid->res+begin,step,MPI_DOUBLE,0,size*mc_tag_data+rank,MPI_COMM_WORLD);
    //fprintf(stderr,"send %d %g %g %g %g %g \n",rank,grid->res[begin+0],grid->res[begin+1],grid->res[begin+2],grid->res[begin+3],grid->res[begin+4]);
  }
  
  free(loc);
  free(pars);
  return;
}

char* pack_data(char* data, char* dest, size_t size, size_t *cursize,size_t checksize, error **err) {
  char * res;
  //fprintf(stderr,">>%d %d %d\n",size,*cursize,checksize);
  memcpy(dest,data,size);
  (*cursize) += size;
  res = dest + size;
  //fprintf(stderr,"%d %d %d<<\n",size,*cursize,checksize);
  testErrorRet((*cursize)>checksize,pmc_io,"Overflowing",*err,__LINE__,NULL);
  //res[0]='\0';
  //fprintf(stderr,"--->>> %s\n",dest);
  return res;
}

void grid_simu_dump(FILE *where,grid_simu *grid,error **err) {
  size_t cursize,maxsize;
  char *buffer,*curbuf,*pbuf;
  
  cursize = 0;
  maxsize =  sizeof(char)*4;
  maxsize += sizeof(int)*(2+grid->ndim);
  maxsize += sizeof(long);
  maxsize += sizeof(short)*(grid->totbins);
  maxsize += sizeof(double)*(grid->ndim*2+grid->totbins);
  //fprintf(stderr,">>>> %d %d\n",cursize,maxsize);
  
  buffer = malloc_err(maxsize*2,err);
  forwardError(*err,__LINE__,);
  
  curbuf = buffer;
  // write header
  //fprintf(stderr,"H\n");
  curbuf = pack_data("GRI1",curbuf,sizeof(char)*4,&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  
  curbuf = pack_data((char*)&(grid->isLog),curbuf,sizeof(int),&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  
  //fprintf(stderr,"ndim %d\n",grid->ndim);  
  curbuf = pack_data((char*)&(grid->ndim),curbuf,sizeof(int),&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  
  //fprintf(stderr,"totb %d\n",grid->totbins);  
  curbuf = pack_data((char*)&(grid->totbins),curbuf,sizeof(long),&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  
  //fprintf(stderr,"limits\n");  
  // write grid description
  curbuf = pack_data((char*)(grid->limits),curbuf,sizeof(double)*grid->ndim*2,&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  
  //fprintf(stderr,"nbin\n");  
  curbuf = pack_data((char*)(grid->nbin),curbuf,sizeof(int)*grid->ndim,&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  
  // write flag
  //fprintf(stderr,"flags\n"); 
  curbuf = pack_data((char*)(grid->flag),curbuf,sizeof(short)*grid->totbins,&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  // write res
  //fprintf(stderr,"res\n"); 
  pbuf=curbuf;
  curbuf = pack_data((char*)(grid->res),curbuf,sizeof(double)*grid->totbins,&cursize,maxsize,err);
  forwardError(*err,__LINE__,);
  //fprintf(stderr,"--->> %g %g\n\n",grid->res[12650],*((double *)(pbuf+12650*sizeof(double))));
  testErrorRet(fwrite(buffer,sizeof(char),cursize,where)!=cursize,pmc_io,"Cannot write to file",*err,__LINE__,);
  
  return;
}

#ifdef _WITH_HDF5_
void grid_simu_hdfdump(grid_simu *grid, char* fname, error **err) {
  hid_t       file_id, group_id;   
  herr_t      status;
  hsize_t     *dims;
  int i;
  
  dims = malloc_err(sizeof(hsize_t)*grid->ndim,err);
  forwardError(*err,__LINE__,);
  
  for(i = 0;i<grid->ndim;i++) {
    dims[i] = grid->nbin[i];
  }
  /* Create a new file using default properties. */
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  testErrorRetVA(file_id<0,hdf5_base,"cannot create hdf5 file %s (got %d)",*err,__LINE__,,fname,file_id);

  // create group grid
  group_id = H5Gcreate( file_id, "grid", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  testErrorRetVA(group_id<0,hdf5_base,"cannot create group /grid in file %s (got %d)",*err,__LINE__,,fname,group_id);
  
  
  // save flags & res
  status = H5LTmake_dataset(group_id,"flag",grid->ndim,dims,H5T_NATIVE_SHORT,grid->flag);
  testErrorRetVA(status<0,hdf5_base,"cannot save /grid/flag in file %s (got %d)",*err,__LINE__,,fname,status);
  status = H5LTmake_dataset(group_id,"log_lkl",grid->ndim,dims,H5T_NATIVE_DOUBLE,grid->res);
  testErrorRetVA(status<0,hdf5_base,"cannot save /grid/data in file %s (got %d)",*err,__LINE__,,fname,status);
  
  
  // set attributes
  status = H5LTset_attribute_int( group_id, ".", "ndim", &(grid->ndim), 1);
  testErrorRetVA(status<0,hdf5_base,"cannot save /grid/ndim in file %s (got %d)",*err,__LINE__,,fname,status);
  status = H5LTset_attribute_int( group_id, ".", "nbin", grid->nbin, grid->ndim);
  testErrorRetVA(status<0,hdf5_base,"cannot save /grid/nbin in file %s (got %d)",*err,__LINE__,,fname,status);
  for(i=0;i<grid->ndim;i++) {
    char cpos[20];
    sprintf(cpos,"pos_%d",i+1);
    status = H5LTset_attribute_double( group_id, ".", cpos, grid->pos[i], grid->nbin[i]);
    testErrorRetVA(status<0,hdf5_base,"cannot save /grid/%s in file %s (got %d)",*err,__LINE__,,cpos,fname,status);
  }

  //close
  status = H5Gclose(group_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save /grid in file %s (got %d)",*err,__LINE__,,fname,status);
  
  status = H5Fclose(file_id);
  testErrorRetVA(status<0,hdf5_base,"cannot save file %s (got %d)",*err,__LINE__,,fname,status);
  
  free(dims);
  
  return;
}
#endif
  
#ifdef _WITH_RC_
grid_simu *init_grid_from_rc(confFile *rc, char * root, error **err) {
  confFile *rca;
  parabox *pb;
  int ndim,i,nbins_size;
  double* mM;
  grid_simu *grid;
  long *lbins;
  int *nbins;
  int islog,isar;
  int usepb,allbin;
  int cas1;
  distribution *target;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  target = init_distribution_from_rc(rca,"target",err);
  forwardError(*err,__LINE__,NULL);
  
  pb = parabox_from_rc(rca,"pb",err);
  forwardError(*err,__LINE__,NULL);

  ndim = pb->ndim;

  // init nbins
  isar = rc_is_array(rca,"nbin",err);
  forwardError(*err,__LINE__,NULL);
  
  nbins = malloc_err(sizeof(int)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  if (isar==1) {
    nbins_size = rc_get_integer_array(rca,"nbin",&lbins,err);
    forwardError(*err,__LINE__,NULL);
    testErrorVA(nbins_size!=ndim,-1,"nbins size (%d) different from ndim (%d) aborting",*err,__LINE__,nbins_size,ndim);
  } else {
    allbin = rc_get_integer(rca,"nbin",err);
    forwardError(*err,__LINE__,NULL);
  }
  for(i=0;i<ndim;i++) {
    if (isar ==1) {
      nbins[i] = lbins[i];
    } else {
      nbins[i] = allbin;
    }
  }  
  
  // cas 1 : pos
  
  cas1 = rc_has_key(rca,"pos",err);
  forwardError(*err,__LINE__,NULL);
  
  if(cas1==1) {
    
    double **pos;
    pos = malloc_err(sizeof(double*)*ndim,err);
    forwardError(*err,__LINE__,NULL);


    isar = rc_is_array(rca,"pos",err);
    forwardError(*err,__LINE__,NULL);
    
    if (isar == 1) {
      for(i=0;i<ndim;i++) {
        int sz;
        char cpos[20];
        sprintf(cpos,"pos[%d]",i+1);
        sz = rc_get_real_array(rca,cpos,&(pos[i]),err);
        forwardError(*err,__LINE__,NULL);
        testErrorRetVA(sz!=nbins[i],pmc_dimension,"too many elements for direction %d (got %d expected %d)",*err,__LINE__,NULL,i,sz,nbins[i]);
      }
    } else {
#ifdef _WITH_HDF5_
      hid_t file_id;
      herr_t hstat;
      int ldim;
      
      char *hdfname;
      hdfname = rc_get_string(rca,"pos",err);
      forwardError(*err,__LINE__,NULL);
      
      file_id = H5Fopen( hdfname, H5F_ACC_RDONLY, H5P_DEFAULT);
      testErrorRetVA(file_id<0,hdf5_base,"cannot open  file %s (got %d)",*err,__LINE__,NULL,hdfname,file_id);
      hstat = H5LTget_attribute_int( file_id, "/grid", "ndim",  &ldim);
      testErrorRetVA(hstat<0,hdf5_base,"cannot read /grid/ndim in file %s (got %d)",*err,__LINE__,NULL,hdfname,hstat);
      testErrorRetVA(ldim != ndim,pmc_dimension,"Bad size for file %s (got %d expected %d)",*err, _LINE__,NULL,hdfname,ldim,ndim);
      for(i=0;i<ndim;i++) {
        hsize_t sz;
        H5T_class_t dumclass;
        size_t dumsz;
        char cpos[20];
        sprintf(cpos,"pos_%d",i+1);
        hstat = H5LTget_attribute_info( file_id, "/grid", cpos, &sz, &dumclass, &dumsz );        
        testErrorRetVA(hstat<0,hdf5_base,"cannot read /grid/%s in file %s (got %d)",*err,__LINE__,NULL,cpos,hdfname,hstat);
        testErrorRetVA(sz!=nbins[i],pmc_dimension,"too many elements for direction %d (got %d expected %d)",*err,__LINE__,NULL,i,sz,nbins[i]);
        pos[i] = malloc_err(sizeof(double)*nbins[i],err);
        forwardError(*err,__LINE__,NULL);
        hstat = H5LTget_attribute_double( file_id, "/grid", cpos, (pos[i]));
        testErrorRetVA(hstat<0,hdf5_base,"cannot read /grid/%s in file %s (got %d)",*err,__LINE__,NULL,cpos,hdfname,hstat); 
      }
      hstat = H5Fclose(file_id);
      testErrorRetVA(hstat<0,hdf5_base,"cannot save file %s (got %d)",*err,__LINE__,,hdfname,hstat);
#else
      testErrorRet(1==1,pmc_io,"Cannot do that without hdf5. Recompile with hdf5 !",*err,__LINE__,NULL);
#endif
    }
    
    grid = init_grid(target, pb,nbins,NULL, 0,pos,err);
    forwardError(*err,__LINE__,NULL);
    if(isar==0) {
      for(i=0;i<ndim;i++) {
        
        free(pos[i]);
      }
    }
    free(pos);
    
  } else {
    mM = malloc_err(sizeof(double)*ndim*2,err);
    forwardError(*err,__LINE__,NULL);

    usepb = rc_safeget_integer(rca,"usepb",0,err);
    forwardError(*err,__LINE__,NULL);

    islog = rc_safeget_integer(rca,"islog",0,err);
    forwardError(*err,__LINE__,NULL);


    for(i=0;i<ndim;i++) {
      double rmin,rmax;
      char par[100];
      if (usepb == 1) {
      sprintf(par,"pb.min[%d]",i+1);
      } else {
      sprintf(par,"min[%d]",i+1);
      }
      rmin = rc_get_real(rca,par,err);
      forwardError(*err,__LINE__,NULL);
      if (usepb == 1) {
      sprintf(par,"pb.max[%d]",i+1);
      } else {
      sprintf(par,"max[%d]",i+1);
      }
      rmax = rc_get_real(rca,par,err);
      forwardError(*err,__LINE__,NULL);
      mM[i*2] = rmin;
      mM[i*2+1] = rmax;
      if (isar ==1) {
        nbins[i] = lbins[i];
      } else {
        nbins[i] = allbin;
      }
    }

    grid = init_grid(target, pb,nbins,mM, islog,NULL,err);
    forwardError(*err,__LINE__,NULL);

    free(nbins);
    free(mM);    
  }
  
  rc_close(&rca);

  return grid;
}
#endif
