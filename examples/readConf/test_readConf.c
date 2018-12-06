/*
 *  test_readConf.c
 *  readConf
 *
 *  Created by Karim Benabed on 19/03/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "readConf.h"

#define __FUNC__ "main"
int main(int argc, char **argv) {
  confFile* rc, *rca;
  error *_err;
  error **err;
  double dres;
  char *sres;
  long ires;
  long *lar;
  size_t n_lar,i;
  char **sar;
  size_t n_sar;
  int d1,d2,d3,d4,d5;
  int go;
  int ii;
  
  _err = initError();
  err = &_err;

  printf("----- starting %s -----\n",argv[0]);
  printf("reading parfile from command line (assuming it's the last parameter of the command line)\n");
  testErrorExit(argc==1,-1,"not enough parameters on the command line",*err,__LINE__);  
  
  /* old way
  
  // initialize parameter parsing with a file
  rc = rc_open(argv[argc-1],err);
  quitOnError(*err,__LINE__,stderr);
  
  // tell readConf to be verbose on stderr  
  rc_toggle_display(rc, stderr, err);
  quitOnError(*err,__LINE__,stderr);
  
  
  // look for modifiers for the parameter file on the command line, prepended with -e or --extra
  rc_do_args(rc,argc,argv,'e',"extra",err);
  quitOnError(*err,__LINE__,stderr);
*/

  // quick init ! Take care of command line and verbosity also
  rc = rc_init_from_args(argc,argv,err);
  quitOnError(*err,__LINE__,stderr);
  
  // read a few parameters
  dres = rc_get_real(rc, "test_double", err);
  quitOnError(*err,__LINE__,stderr);
  fprintf(stderr,"read test_double -> %g\n",dres);
  
  
  ires = rc_get_integer(rc, "test_int", err);
  quitOnError(*err,__LINE__,stderr);
  fprintf(stderr,"read test_int -> %ld\n",ires);

 sres = rc_get_string(rc, "test_string", err);
 quitOnError(*err,__LINE__,stderr);
 fprintf(stderr,"read test_string -> %s\n",sres);
  
  dres = rc_get_real(rc, "test_object.test_double", err);
  quitOnError(*err,__LINE__,stderr);
  fprintf(stderr,"read test_object.test_double -> %g\n",dres);
    
  
  
  n_lar = rc_get_integer_array(rc, "test_object.test_array_integer", &lar, err);
  quitOnError(*err,__LINE__,stderr);
  for(i=0;i<n_lar;i++) {
    fprintf(stderr,"test_object.test_array_integer[%ld] = %ld\n",i,lar[i]);
  }
  
  n_sar = rc_get_string_array(rc, "test.string_array", &sar, err);
  quitOnError(*err,__LINE__,stderr);
  for(i=0;i<n_sar;i++) {
    fprintf(stderr,"test_object.test_string_array[%ld] = '%s'\n",i,sar[i]);
  }
  
  rca = rc_alias(rc,"test",err);
  quitOnError(*err,__LINE__,stderr);
  n_sar = rc_get_string_array(rca, "string_array", &sar, err);
  quitOnError(*err,__LINE__,stderr);
  for(i=0;i<n_sar;i++) {
    fprintf(stderr,"test_object.test_string_array[%ld] = '%s'\n",i,sar[i]);
  }
  rc_close(&rca);
  
  fprintf(stderr,"now for somethong completly different...\n");
  
  // read lot's of parameters in one function !
  go = rc_get_list(rc, 5, "", -1, err,
              "oest_double","g",&dres,&d1,
              "test_string","s",&sres,&d2,
              "test_int","d",&ires,&d3,
              "test_object.test_double","g",&dres,&d4,
              "test_object.test_array_integer","d*",&lar,&n_lar,&d5);
  quitOnError(*err,__LINE__,stderr);
  fprintf(stderr,"got %d params %d %d %d %d %d\n",go,d1,d2,d3,d4,d5);
  
  go = rc_get_list(rc, 2, "test_object", 0, err,
              "test_double","g",&dres,
              "test_array_integer","d*",&lar,&n_lar);
  quitOnError(*err,__LINE__,stderr);
  //fprintf(stderr,"got %d params %d %d %d %d %d\n",go,d1,d2,d3,d4,d5);
  
  // I do not need parameters anymore, close rc
  rc_close(&rc);
  endError(err);
  
  return 0;
}
#undef __FUNC__