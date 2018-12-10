/*
 *  tools.c
 *  likely
 *
 *  Created by Karim Benabed on 13/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "tools.h"



int print_if_err(error **err, FILE *F1, FILE *F2)
{
   if (isError(*err)) {
      printError(F1, *err);
      fprintf(F1, "continuing...\n");
      fflush(F1);

      printError(F2, *err);
      fprintf(F2, "continuing...\n");
      fflush(F2);

      purgeError(err);
      return 1;
   }

   return 0;
}

void out_err_log(FILE *FLOG, const char* str, ...)
{
  va_list v;
 
  va_start(v, str);
  //vfprintf(stderr, str, v);
  vfprintf(FLOG, str, v);
  va_end(v);
}

/* On my Mac there is no getline... */

#ifndef SIZE_MAX
# define SIZE_MAX ((size_t) -1)
#endif
#ifndef SSIZE_MAX
# define SSIZE_MAX ((ssize_t) (SIZE_MAX / 2))
#endif
#if !HAVE_FLOCKFILE
# undef flockfile
# define flockfile(x) ((void) 0)
#endif
#if !HAVE_FUNLOCKFILE
# undef funlockfile
# define funlockfile(x) ((void) 0)
#endif

ssize_t getdelim(char **lineptr, size_t *n, int delimiter, FILE *fp)
{
   ssize_t result=0;
   size_t cur_len = 0;

   if (lineptr == NULL || n == NULL || fp == NULL)
   {
      errno = EINVAL;
      return -1;
   }

   flockfile (fp);

   if (*lineptr == NULL || *n == 0)
   {
      *n = 120;
      *lineptr = (char *) malloc (*n);
      if (*lineptr == NULL)
      {
	 result = -1;
	 goto unlock_return;
      }
   }

   for (;;)
   {
      int i;

      i = getc (fp);
      if (i == EOF)
      {
	 result = -1;
	 break;
      }

      /* Make enough space for len+1 (for final NUL) bytes.  */
      if (cur_len + 1 >= *n)
      {
	 size_t needed_max =
	   SSIZE_MAX < SIZE_MAX ? (size_t) SSIZE_MAX + 1 : SIZE_MAX;
	 size_t needed = 2 * *n + 1;   /* Be generous. */
	 char *new_lineptr;

	 if (needed_max < needed)
	   needed = needed_max;
	 if (cur_len + 1 >= needed)
	 {
	    result = -1;
	    goto unlock_return;
	 }

	 new_lineptr = (char *) realloc (*lineptr, needed);
	 if (new_lineptr == NULL)
	 {
	    result = -1;
	    goto unlock_return;
	 }

	 *lineptr = new_lineptr;
	 *n = needed;
      }

      (*lineptr)[cur_len] = i;
      cur_len++;

      if (i == delimiter)
	break;
   }
   (*lineptr)[cur_len] = '\0';
   result = cur_len ? cur_len : result;

 unlock_return:
   funlockfile (fp);
   return result;
}

ssize_t getline(char **lineptr, size_t *n, FILE *stream)
{
   return getdelim(lineptr, n, '\n', stream);
}

void read_header_lc(const char *line, int npar, int n_ded, error **err)
{
   int NPAR, N_DED, nread;

   nread = sscanf(line, "# npar = %d, n_ded = %d\n", &NPAR, &N_DED);

   if (nread==2) {
      testErrorRetVA(NPAR!=npar, tls_Nparam, "Chain/sample header (npar=%d) not compatible with config file (npar=%d)",
		     *err, __LINE__,, NPAR, npar);
      testErrorRetVA(N_DED!=n_ded, tls_Nparam, "Chain/sample header (n_ded=%d) not compatible with config file (n_ded=%d)",
		     *err, __LINE__,, N_DED, n_ded);
   } else {
      /* Different header line */
   }

   /*
   nread = sscanf(line, "# number 1 %d %s\n", &nparam, dummy);
   if (n_ded>0 && nread==2 && strcmp(dummy, "param_ded")==0) {
      testErrorRet(nparam!=n_ded, tls_Nparam, "Chain/sample not compatible with config file "
		   "(different number of deduced parameters)", *err, __LINE__,);
   } else if (nread==2 && strcmp(dummy, "param")==0) {
      testErrorRetVA(nparam!=npar, tls_Nparam, "Chain/sample (%d parameters) not compatible with config file (npar=%d)",
		     *err, __LINE__,, nparam, npar);
   }
   */
}

/* CHPRE has to be positioned after header. Used in mkmc. Returns non-zero value if end of file reached. */
int read_next_mkmc_step(FILE *CHPRE, int npar, int n_ded, double *pstate, double *pstate_ded, 
			double *logL, double *accept, error **err)
{
   char *line, *line0;
   size_t nread, length;
   int eof;

   do {
      line = NULL;
      nread = getline(&line, &length, CHPRE);
      line0 = line;          /* remember starting point for free() */
      /* nread=-1 also for eof */
      /*testErrorRet(nread==-1, mcmc_file,
	"Error while reading mcmc chain file. Maybe input/output error.",
	*err, __LINE__, -1);*/
      if (nread==-1) {
	 fprintf(stderr, "Eof in previous chain/sample, ending reading\n");
	 return feof(CHPRE);
      }
   } while (line[0]=='#');

   testErrorRet(line[0]=='#', mcmc_file, "Header line not expected that late in chain/sample file.",
		*err, __LINE__, -1);

      read_step_lc(&line, npar, n_ded, pstate, pstate_ded, logL, accept, err); forwardError(*err, __LINE__, -1);
      if (line0) free(line0);

   eof = feof(CHPRE);
   //fprintf(stderr, "eof = %d\n", eof);
   return eof;
}

/* Reads mcmc parameter from line, using functions from sn1a.c */
int read_step_lc(char **line, int npar, int n_ded, double *param, double *param_ded, double *logL, double *accept,
		 error **err)
{
   int ok, i;

   ok = read_double(line, accept);
   testErrorRet(ok==0, tls_file, "Bad line: First entry (mcmc:accept, pmc:weight) not a double.",
		*err, __LINE__, ok);
   ok = read_double(line, logL);
   *logL = -(*logL);
   testErrorRet(ok==0, tls_file, "Bad line: Second entry (mcmc:logL, pmc:icomp) not a double",
		*err, __LINE__, ok);
   for (i=0; i<npar; i++) {
      ok = read_double(line, param+i);
      testErrorRet(ok==0, tls_file, "Bad line: Parameter entry not a double", *err, __LINE__, ok);
   }
   if (n_ded>0) testErrorRet(param_ded==NULL, mcmc_ded, "param_ded is NULL but n_ded>0", *err, __LINE__, -1);
   for (i=0; i<n_ded; i++) {
      ok = read_double(line, param_ded+i);
      testErrorRet(ok==0, tls_file, "Bad line: Deduced parameter not a double", *err, __LINE__, ok);
   }
   return ok;
}

/* ============================================================= *
 * Reads a MCMC chain or PMC simulation (lc format). If weight   *
 * is not NULL, it is filled with the 'weight' field. If indices *
 * is not NULL, it is filled with (size_t)(-logL) = mixture      *
 * component.							 *
 * ============================================================= */
double *read_mkmc_chain(int npar, int n_ded, int nchainmax, FILE *CHAIN, long *naccepted,
			fpos_t *headpos, double *param_ded, double *weight, size_t *indices, error **err)
{
   size_t length = 0;
   ssize_t nread;
   char *line, *line0;
   double *accstates, logL, accept;

   //fprintf(FLOG, "read_mkmc_chain:\n"); fprintf(FLOG, "================\n");

   accstates = calloc_err(npar*nchainmax, sizeof(double), err);
   forwardError(*err, __LINE__, NULL);
   *naccepted = 0;

   do {

      line = NULL;
      nread = getline(&line, &length, CHAIN);
      line0 = line;          /* remember starting point for free() */
      if (feof(CHAIN)) break;
      testErrorRet(nread==-1, mcmc_file,
		   "Error while reading chain/sample file. Maybe input/output error",
		   *err, __LINE__, NULL);

      if (line[0]=='#') {

	 read_header_lc(line, npar, n_ded, err);
	 forwardError(*err, __LINE__, NULL);
	 fgetpos(CHAIN, headpos);

      } else {

	 /* data */
	 testErrorRetVA(*naccepted>=nchainmax, tls_overflow, "Chain/sample length is larger than expected (n=%d). "
			"For mcmc: check nchain (in config file). "
			"For pmc: check nsamples, or nsamples*fsfinal in case of last iteration.",
			*err, __LINE__, NULL, nchainmax);

	 read_step_lc(&line, npar, n_ded, accstates+npar*(*naccepted),
		      param_ded+n_ded*(*naccepted), &logL, &accept, err);
	 forwardError(*err, __LINE__, NULL);
	 if (weight!=NULL) {
	    /* Weights are stored in file as log(weights) */
	    //weight[*naccepted]   = exp(accept);
	    /* New (cosmo_pmc) v1.2: Avoid inf due to very large log-posteriors */
	    weight[*naccepted] = accept;
	 }
	 if (indices!=NULL) indices[*naccepted] = (size_t)(logL);
	 (*naccepted)++;

      }

      if (line0) free(line0);

   } while (1);

   //fprintf(FLOG, "naccepted = %d\n\n", *naccepted);

   return accstates;
}

void print_parameter(FILE *where, size_t npar, const double *params)
{
   int i;
   for (i=0;i<npar;i++) {
      fprintf(where, "%g ", params[i]);
   }
   fprintf(where, "\n");
}

void print_step(FILE* where, int accept, double loglkl, size_t ndim, double *params) {
  FILE *rf;
  size_t i;
  
  rf=where;
  if (where==NULL)
    rf=stdout;
  fprintf(rf,"%1d %10g",accept,loglkl);
  for (i=0;i<ndim;i++) {
    fprintf(rf," %10g",params[i]);
  }
  fprintf(rf,"\n");
}

size_t read_chain(size_t *accepted,int accept_flag, double *loglk, double* params, int nsamples, int ndim, FILE* where, error **err) {
  FILE *rf;
  size_t i;
  double _loglk;
  int _accept,nread;
  double *_params;
  size_t lread;
  
  rf=where;
  if (where==NULL)
    rf=stdin;
  lread=0;
  
  while (lread<nsamples) {
    nread=fscanf(rf,"%d %lg ",&_accept,&_loglk);
    if (nread!=2) {
      /* assume I have finished */
      break;
    }
    _params=&(params[lread*ndim]);
    for(i=0;i<ndim-1;i++) {
      testErrorRetVA(fscanf(rf,"%lg",&(_params[i]))!=1,tls_io,"Bad format for file at index %d dim %d",*err,__LINE__,0,lread,i);
      //fprintf(stderr,"%g ",_params[i]);
    }
    testErrorRetVA(fscanf(rf,"%lg\n",&(_params[i]))!=1,tls_io,"Bad format for file at index %d dim %d",*err,__LINE__,0,lread,i);
    //fprintf(stderr,"%g\n",_params[i]);
    if ((accept_flag==MC_AF_ALL) || ((accept_flag==_accept))) {
      if (accepted!=NULL) {
        accepted[lread]=_accept;
      }
      if (loglk!=NULL) {
        loglk[lread]=_loglk;
      }
      lread++;
    }
  }
  return lread;
}

repar* repar_init(int nin,int nout,void* payload, posterior_log_pdf_func *payload_pdf, posterior_log_free* payload_free,void* repar_func_data, repar_func_type* rfnt, posterior_log_free* repar_func_free, error **err) {
  repar* self;
  
  self = malloc_err(sizeof(repar),err);
  forwardError(*err, __LINE__,NULL);
  
  self->nin=nin;
  self->nout=nout;
  self->rpars = malloc_err(sizeof(double)*nout,err);
  forwardError(*err, __LINE__,NULL);
  self->repar_func = rfnt;
  self->repar_func_data = repar_func_data;
  self->repar_func_free = repar_func_free;
  self->payload = payload;
  self->payload_pdf = payload_pdf;
  self->payload_free = payload_free;
  return self;
}

double repar_lkl(void* pelf,double* pars, error **err) {
  repar* self;
  double* rpars;
  double res;
  
  self=pelf;
  rpars=self->rpars;
  self->repar_func(self->repar_func_data,pars,rpars,err);
  forwardError(*err, __LINE__,0);
  
  res = self->payload_pdf(self->payload,rpars,err);
  forwardError(*err, __LINE__,0);
  
  return res;  
}

void repar_free(repar** pelf) {
  repar* self;
  self=*pelf;
  free(self->rpars);
  
  if(self->payload_free!=NULL)
    self->payload_free(&(self->payload));
  if(self->repar_func_free!=NULL)
    self->repar_func_free(&(self->repar_func_data));
  free(self);
  *pelf=NULL;
}

repar_select* repar_select_init(int nin, int nout, int* select,double* base,error **err) {
  repar_select* self;
  int i,ns,no;
  
  self = malloc_err(sizeof(repar_select),err);
  forwardError(*err, __LINE__,NULL);
  self->nin=nin;
  self->nout=nout;
  self->mode = REPAR_SELECT_MODE_EXPAND;
  ns=nout;
  no=nin;
  if (nin<nout) {
    self->mode = REPAR_SELECT_MODE_REDUCE;
    ns=nin;
    no=nout;
    testErrorRet(base!=NULL,tls_undef,"Base has no meaning in reduce mode",*err,__LINE__,NULL);
  }
  self->select = malloc_err(sizeof(int)*ns,err);
  forwardError(*err, __LINE__,NULL);
  
  
  for(i=0;i<ns;i++) {
    int ml;
    ml = select[i];
    testErrorRet(ml>=no,tls_outOfBound,"Out of bounds !",*err,__LINE__,NULL);
    self->select[i] = ml; 
  }
  
  self->base =NULL;
  if (base!=NULL) {
    self->base = malloc_err(sizeof(int)*nout,err);
    forwardError(*err, __LINE__,NULL);
    memcpy(self->base, base, sizeof(double)*no);
  }
  
  return self;
}

void repar_select_func(void* pelf, double* parin, double* parout, error **err) {
  int i;
  repar_select *self;
  
  self=pelf;
  
  if (self->mode == REPAR_SELECT_MODE_REDUCE) {
    for(i=0;i<self->nout;i++) {
      parout[i]=parin[self->select[i]];
    }
  } else {
    memcpy(parout,self->base,sizeof(double)*self->nout);
    for(i=0;i<self->nin;i++) {
      parout[self->select[i]]=parin[i];
    }
  }
  return;
}

void repar_select_free(void **pelf) {
  repar_select *self;
  self=*pelf;
  free(self->select);
  if (self->base!=NULL) {
    free(self->base);
  }
  free(self);
  *pelf=NULL;
}
