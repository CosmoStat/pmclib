/*
 *  optimize.c
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 09/10/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifdef __PLANCK__
#include "HL2_likely/pmclib/optimize.h"
#else
#include "optimize.h"
#endif


optimize_struct * optimize_init(distribution *target,void *data, optimize_func *optf, optimize_get_func *get_p, optimize_get_func *get_v,posterior_log_free *fropt, error **err) {
  optimize_struct *opt;
  
  testErrorRet(optf==NULL,-10001000,"optimize slot must be filled",*err,__LINE__,NULL);
  testErrorRet(get_p==NULL,-10001000,"get_result slot must be filled",*err,__LINE__,NULL);
  opt = malloc_err(sizeof(optimize_struct),err);
  forwardError(*err,__LINE__,NULL);
  opt->ndim = target->ndim;
  opt->target = target;
  opt->data  = data;
  opt->optimize = optf;
  opt->get_result = get_p;
  opt->get_variance = get_v;
  opt->free = fropt;
  opt->cvg = 0;
  return opt;
}

void optimize(optimize_struct *opt,error **err) {
  opt->optimize(opt,err);
  forwardError(*err,__LINE__,);
  if (isError(*err)) {
    return;
  }
  opt->cvg=1;
}

double* optimize_get_best_pars(optimize_struct *opt, error **err) {
  double* res;
  
  testErrorRet(opt->cvg==0,-10001000,"did not converge",*err,__LINE__,NULL);
  res = opt->get_result(opt,err);
  forwardError(*err,__LINE__, NULL);
  
  return res;
}

double* optimize_get_variance_default(optimize_struct *opt, error **err) {
  double *res,*pos;
  int ndim;
  int i,j;
  char uplo;
  int info;
  
#ifndef _WITH_LAPACK_
  testErrorRet(1==1,-1,"Cannot do without lapack !",*err,__LINE__,NULL);
#endif
  ndim = opt->ndim;

  pos = optimize_get_best_pars(opt,err);
  forwardError(*err,__LINE__,NULL);
  
  res = distribution_second_derivative_matrix(opt->target,pos,err);
  forwardError(*err,__LINE__,NULL);
  free(pos);
  
  for(i=0;i<ndim;i++) {
    for(j=0;j<ndim;j++) {
      res[i+j*ndim] = - res[i+j*ndim];
    }
  }
  
  uplo = 'L';
#ifdef _WITH_LAPACK_
  dpotrf(&uplo,&ndim,res,&ndim,&info);
  testErrorRetVA(info!=0,-101001010,"Could not cholesky decompose using dpotrf (status %d)",*err,__LINE__,,info);

  dpotri(&uplo, &ndim, res, &ndim, &info );
#endif
  
  for(i=0;i<ndim;i++) {
    for(j=i+1;j<ndim;j++) {
      res[i+j*ndim] = res[j+i*ndim];
    }
  }
  
  return res;
}

double* optimize_get_variance(optimize_struct *opt, error **err) {
  double* res;
  

  testErrorRet(opt->cvg==0,-10001000,"did not converge",*err,__LINE__,NULL);
  //testErrorRet(opt->get_variance==NULL,-10001000,"can't get variance",*err,__LINE__,NULL);
  if (opt->get_variance==NULL) {
    res = optimize_get_variance_default(opt,err);
    forwardError(*err,__LINE__,NULL);
    return res;
  }
  res = opt->get_variance(opt,err);
  forwardError(*err,__LINE__,NULL);
  
  return res;
}

void optimize_free(void **popt) {
  optimize_struct *opt;
  
  opt = *popt;
  if (opt->free) {
    opt->free(&(opt->data));
  } 
  
  free(opt);
  *popt = NULL;
}

#ifdef _WITH_RC_

void *init_optimize_from_rc(confFile *rc, char * root, error **err) {
  optimize_struct *target;
  confFile *rca;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);

  INIT_FROM_RC(rca,target);
  
  rc_close(&rca);
   
  return target;
}

void *init_linesearch_from_rc(confFile *rc, char *root, error **err) {
  linesearch_struct *target;
  confFile *rca;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);

  INIT_FROM_RC(rca,target);

  rc_close(&rca);

  return target;
}

void* rcinit_bydir(confFile *rc, char *root, error **err) {
  confFile *rca;
  optimize_struct *opt;
  distribution *target;
  int ndum;
  double tol;
  int maxstep;
  linesearch_struct *ls;
  double *scales,*guess;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  target = init_distribution_from_rc(rca,"target",err);
  forwardError(*err,__LINE__,NULL);

  ndum = rc_get_real_array(rca,"guess",&guess,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(ndum!=target->ndim,-10001000,"distribution has %d dimension while guess has %d",*err,__LINE__,NULL,target->ndim,ndum);
  
  scales= NULL;
  ndum = rc_safeget_real_array(rca,"scales",&scales,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(ndum!=target->ndim && scales!=NULL,-10001000,"distribution has %d dimension while scales has %d",*err,__LINE__,NULL,target->ndim,ndum);
  
  tol = rc_safeget_real(rca,"tol",1e-6,err);
  forwardError(*err,__LINE__,NULL);
  
  maxstep = rc_safeget_integer(rca,"maxstep",100,err);
  forwardError(*err,__LINE__,NULL);
  
  ls = init_linesearch_from_rc(rca,"linesearch",err);
  forwardError(*err,__LINE__,NULL);

  opt = bydir_init(target, guess, scales,ls, tol, maxstep,err);
  forwardError(*err,__LINE__,NULL);

  rc_close(&rca);

  return opt;
}

#ifdef _WITH_LAPACK_
void* rcinit_bfgs(confFile *rc, char *root, error **err) {
  confFile *rca;
  optimize_struct *opt;
  distribution *target;
  int ndum;
  double tol;
  int maxstep;
  linesearch_struct *ls;
  double *scales,*guess;
  
  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  target = init_distribution_from_rc(rca,"target",err);
  forwardError(*err,__LINE__,NULL);
  
  ndum = rc_get_real_array(rca,"guess",&guess,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(ndum!=target->ndim,-10001000,"distribution has %d dimension while guess has %d",*err,__LINE__,NULL,target->ndim,ndum);
  
  scales= NULL;
  ndum = rc_safeget_real_array(rca,"scales",&scales,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(ndum!=target->ndim && scales!=NULL,-10001000,"distribution has %d dimension while scales has %d",*err,__LINE__,NULL,target->ndim,ndum);
  
  tol = rc_safeget_real(rca,"tol",1e-6,err);
  forwardError(*err,__LINE__,NULL);

  
  maxstep = rc_safeget_integer(rca,"maxstep",100,err);
  forwardError(*err,__LINE__,NULL);
  
  ls = init_linesearch_from_rc(rca,"linesearch",err);
  forwardError(*err,__LINE__,NULL);
  
  opt = bfgs_init(target, guess, scales,ls, tol, maxstep,err);
  forwardError(*err,__LINE__,NULL);
  
  rc_close(&rca);
  
  return opt;
}
#endif

void* rcinit_secant_linesearch(confFile *rc, char *root, error **err) {
  confFile *rca;
  linesearch_struct *ls;
  int ndim;
  double tol,step;
  int maxstep;

  rca = rc_alias(rc,root,err);
  forwardError(*err,__LINE__,NULL);
  
  tol = rc_safeget_real(rca,"tol",1e-6,err);
  forwardError(*err,__LINE__,NULL);

  step = rc_safeget_real(rca,"step",.01,err);
  forwardError(*err,__LINE__,NULL);
  
  maxstep = rc_safeget_integer(rca,"maxstep",100,err);
  forwardError(*err,__LINE__,NULL);
  
  ndim = rc_get_integer(rca,"ndim",err);
  forwardError(*err,__LINE__,NULL);
  
  ls = secant_linesearch_init( ndim, step, maxstep, tol, err);
  forwardError(*err,__LINE__,NULL);
  
  rc_close(&rca);
  
  return ls;
}
#endif

linesearch_struct* linesearch_init(int ndim, void* data, ls_func* fnc, posterior_log_free * frr, error **err) {
  linesearch_struct *ls;
  
  ls = malloc_err(sizeof(linesearch_struct),err);
  forwardError(*err,__LINE__,NULL);
  
  ls->ndim = ndim;
  ls->data = data;
  ls->ls_fnc = fnc;
  ls->ls_free = frr;
  ls->PIM = malloc_err(sizeof(double)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  return ls;
}

double linesearch(linesearch_struct *ls, distribution *dist, double *PIM, double *pk, double *gk,error **err) {
  double res;
  
  res = ls->ls_fnc(ls, dist, PIM, pk, gk,err);
  forwardError(*err,__LINE__,-1);
  return res;
}

void free_linesearch(void **pls) {
  linesearch_struct *ls;
  
  ls = *pls;
  ls->ls_free(&ls->data);
  free(ls->PIM);
  free(ls);
  *pls = NULL;
}

linesearch_struct* secant_linesearch_init(int ndim, double step, int maxstep, double tol, error **err) {
  linesearch_secant_struct *lst;
  linesearch_struct *ls;
  
  lst = malloc_err(sizeof(linesearch_secant_struct),err);
  forwardError(*err,__LINE__,NULL);
  
  lst->step = step;
  lst->maxstep = maxstep;
  lst->tol = tol;
  
  lst->golden = 1;
  
  ls = linesearch_init(ndim,lst, secant_linesearch, secant_free,err);
  forwardError(*err,__LINE__,NULL);
  
  return ls;
}

void secant_free(void** psc) {
  linesearch_secant_struct *sc;

  sc = *psc;  
  free(sc);
  *psc = NULL;
}


void addalphap(int ndim, double *PIM, double alpha, double* pk) {
  int i;
  if (alpha == 0) {
    return;
  }
  for (i=0;i<ndim;i++) {
    PIM[i] += alpha * pk[i];
  }
  return;
}

double deriv_along(distribution *dist, double *PIM, double *pk, double *gk,int update, error **err) {
  int i;
  double r;
  r = 0;
  if (update == 0 && gk ==NULL) {
    r = distribution_deriv_along(dist,PIM,pk,err);
    forwardError(*err,__LINE__,-1); 
    return r;
  }
  for (i=0;i<dist->ndim;i++) {
    double g;
    if (gk==NULL || update == 1) {
      if((pk[i]==0) && (update == 0 || gk==NULL)) {
        g = 0;
      } else {
	//_DEBUGHERE_("der %d %g",i,PIM[i]);
        g = distribution_first_derivative(dist,i,PIM,err);
        forwardError(*err,__LINE__,-1); 
       	//_DEBUGHERE_("der %d %g %g",i,PIM[i],g);

      }
      if (gk!=NULL) {
        gk[i] = g;
      }
    } else {
      g = gk[i];
    }
    //_DEBUGHERE_("up r %g %g %g %g",r,g,pk[i],r+g*pk[i]);
    r +=  g * pk[i];
  }
  _DEBUGHERE_("deriv : %g",r);
  return r;
}

double dicho_along(distribution *dist, double* PIM, double *pk,double *gk, int update, double *pxc, double *pxm, double *pxp, error **err) {
  double xplus,xmoins,xmid, xcur,delta,yn;
  
  xplus = *pxp;
  xmoins = *pxm;
  xmid = (xplus+xmoins)/2.;
  xcur = *pxc;
  
  //_DEBUGHERE_("%g %g %g",xmoins,xmid,xplus);
  delta = xmid - xcur;
  addalphap(dist->ndim,PIM,delta,pk);
  *pxc = xmid;

  yn = deriv_along(dist,PIM,pk,gk,update,err);
  forwardError(*err,__LINE__,-1);
  
  //_DEBUGHERE_(" --> %g",yn);

  
  if (yn>0) {
    *pxp = xmid;
    *pxm = xmoins;
    return yn;
  } else {
    *pxm = xmid;
    *pxp = xplus;
    return yn;
  }
}

#define GOLD_PHI (1.+sqrt(5.))/2.
#define GOLD_RPHI (2. - GOLD_PHI)

double brent_along(distribution *dist, double *PIM, double *pk, double xm,double xc,double xp, double fm, double fc, double fp, double xtol, double ftol, int istep, int maxstep, error **err) {
  double xn,xn2;
  double fn,xdum,r,xco,fbest,fnew;
  
  xco = xc;
  fbest = fm;
  if (fc>fbest) {
    fbest = fc;
  }
  if (fp>fbest) {
    fbest = fp;
  }
  
  while (istep < maxstep) {
    if (fabs(xm-xp)<xtol*(fabs(xm)+fabs(xp))/2.) {
      return xc;
    }
    xn = xp - .5* ((xp-xm)*(xp-xm) * (fp - fc) - (xp-xc)*(xp-xc) * (fp - fm)) / ((xp-xm) * (fp - fc) - (xp-xc) * (fp - fm));
    if (fabs(xn-xp)>fabs(xm-xp) || fabs(xn-xm)>fabs(xm-xp) || fabs(xn-xco)<1e-10) {
      // je tape en dehors, refuse ce deplacement !
      if (fabs(xc-xm) > fabs(xc-xp)) {
        // je choisit le petit cote
        xn = xp;
        xp = xm;
        xm = xn;
        xn = fp;
        fp = fm;
        fm = xn;
      }
      xn2 = 1./GOLD_PHI * (xc-xm) +xc; // saut gold like
    _DEBUGHERE_("refuse brent %g - > %g ",xn,xn2);
    xn = xn2;
    }
   
    addalphap(dist->ndim, PIM, xn-xco, pk);
    fn = distribution_lkl(dist,PIM,err);
    forwardError(*err,__LINE__,-1);    
    fnew = fn;
        
    xco = xn;
    if ((xn-xc)/(xp-xm)<0) {
      xdum = xc;
      xc = xn;
      xn = xdum;
      xdum = fc;
      fc = fn;
      fn = xdum;
    }
    _DEBUGHERE_("%g ** %g ** %g ** %g", xm,xc,xn,xp);
    
    _DEBUGHERE_("%g %g %g ",fc,fn,fn-fc);
   
    if (fabs(fnew-fbest)<ftol) {
      return xco;
    }

    if (fn>fbest) {
      fbest = fn;
    }

    if (fn<fc) {
      //addalphap(dist->ndim, PIM, -xn+xco, pk);
      xdum = xm;
      xm = xn;
      xp = xdum;
      //r = golden_along(dist, PIM, pk, xn, xc,xm,fc,xtol,ftol,err);
      //forwardError(*err,__LINE__,-1);
    } else {
      xm = xc;
      xc = xn;
      fc = fn;
      //r=golden_along(dist, PIM, pk, xc, xn,xp,fn,xtol,ftol,err);
      //forwardError(*err,__LINE__,-1);
    }
    istep++;
  }
  testErrorRetVA(1==1,-111000111,"too many steps (%d) and optimization is still poor (%g %g)",*err,__LINE__,-1,maxstep,fn,fc);
}

double golden_along(distribution *dist, double *PIM, double *pk, double xm,double xc, double xp, double fc, double xtol, double ftol,int istep, int maxstep, error **err) {
  double xn,fn,r,xdum;
  
  while (istep < maxstep) {
    if (fabs(xm-xp)<xtol*(fabs(xm)+fabs(xp))/2.) {
      return xc;
    }
    xn = xc + GOLD_RPHI * (xp-xc);
    _DEBUGHERE_("%g ** %g ** %g ** %g", xm,xc,xn,xp);
    
    addalphap(dist->ndim, PIM, xn-xc, pk);
    
    fn = distribution_lkl(dist,PIM,err);
    forwardError(*err,__LINE__,-1);
    
    _DEBUGHERE_("%g %g",fc,fn);
    
    if (fabs(fn-fc)<ftol) {
      return xn;
    }
    if (fn<fc) {
      addalphap(dist->ndim, PIM, -xn+xc, pk);
      xdum = xm;
      xm = xn;
      xp = xdum;
      //r = golden_along(dist, PIM, pk, xn, xc,xm,fc,xtol,ftol,err);
      //forwardError(*err,__LINE__,-1);
    } else {
      xm = xc;
      xc = xn;
      fc = fn;
      //r=golden_along(dist, PIM, pk, xc, xn,xp,fn,xtol,ftol,err);
      //forwardError(*err,__LINE__,-1);
    }
    istep++;
  }
  testErrorRetVA(1==1,-111000111,"too many steps (%d) and optimization is still poor (%g %g)",*err,__LINE__,-1,maxstep,fn,fc);
}
  
double secant_linesearch(linesearch_struct *lst, distribution *dist, double *PIM, double *pk, double *gk,error **err) {
  double yn,ynm;
  double xn,xnm;
  double delta;
  linesearch_secant_struct *sc;
  int istep;
  double yplus,xplus,ymoins,xmoins,fplus,fmoins;
  double fnm,fn;
  int candic;
  
  sc = lst->data;
  
  memcpy(lst->PIM,PIM,sizeof(double)*lst->ndim);
  {
    char spk[4000];
    char spim[4000];
    int i;
    spk[0]='\0';
    spim[0]='\0';
    for(i=0;i<dist->ndim;i++) {
      sprintf(spk,"%s %g",spk,pk[i]);
      sprintf(spim,"%s %g",spim,PIM[i]);
    }
    _DEBUGHERE_("linsearch from [%s ] toward [%s ]",spim,spk);
  }
  ynm = deriv_along(dist,lst->PIM,pk,gk,0,err);
  forwardError(*err,__LINE__,-1);
  
  fnm = distribution_lkl(dist,lst->PIM,err);
  forwardError(*err,__LINE__,-1);
  
  xn = 0;
  delta = sc->step;
  istep = 0;
  yplus = -1;
  xplus = 0;
  ymoins = 1;
  xmoins = 0;
  while(istep < sc->maxstep) {
    //_DEBUGHERE_("%d -> %g %g %g",istep,xn,ynm,xn+delta);
    //_DEBUGHERE_("+- : %g %g %g %g",xplus,yplus,xmoins,ymoins);

    if (ynm>0) {
      if (ynm<yplus || yplus<0) {
        yplus = ynm;
        xplus = xn;
        fplus = fnm;
      }
    }
    if (ynm<0) {
      if (ynm>ymoins || ymoins>0) {
        ymoins = ynm;
        xmoins = xn;
        fmoins = fnm;
      }
    }
    xnm = xn;
    
    candic = (yplus>=0 && ymoins<=0); 
  
    if (candic && sc->golden==1) {
      double xp,xm,fc,xc,fp,fm;
      if (fplus<fmoins) {
        xp = xplus;
        fp = fplus;
        xm = xmoins;
        fm = fmoins;
      }  else {
        xm = xplus;
        fm = fplus;
        xp = xmoins;
        fp = fmoins;
      } 
      xc = (xp + GOLD_PHI * xm) / (1. + GOLD_PHI);
      addalphap(lst->ndim,lst->PIM,xc-xn,pk);
      fc = distribution_lkl(dist,lst->PIM,err);
      forwardError(*err,__LINE__,-1);
      _DEBUGHERE_("%g ** %g ** %g",xm,xc,xp);      
      xn = brent_along(dist, lst->PIM, pk, xm,xc, xp,fm,fc,fp, 1e-8,sc->tol, istep,sc->maxstep,err);
      //xn = golden_along(dist, lst->PIM, pk, xm,xc, xp,fc, 1e-8,sc->tol, istep,sc->maxstep,err);
      forwardError(*err,__LINE__,-1);
      yn = deriv_along(dist,lst->PIM,pk,gk,1,err);
      forwardError(*err,__LINE__,-1);
      return xn;
    }
    
    if ( candic && ((xn+delta-xplus)/(xplus-xmoins)>= 0 || (xn+delta-xmoins)/(xplus-xmoins)<= 0) ) {
      double xp;
      _DEBUGHERE_("DICHOTOMIZE !","");
      xp = xn;
      yn = dicho_along(dist, lst->PIM, pk, gk, 0, &xn, &xmoins, &xplus, err);
      forwardError(*err,__LINE__,-1);
      //_DEBUGHERE_("xn,yn : %g %g",xn,yn);
      delta = xn - xp;
    } else {
      int imv;
      imv = 0;
      while(imv<10) {
        purgeError(err);
        xn += delta;
        addalphap(lst->ndim,lst->PIM,delta,pk);
        yn = deriv_along(dist,lst->PIM,pk,NULL,0,err);
        if (isError(*err)) {
          _DEBUGHERE_("delta %g too big",delta);
          addalphap(lst->ndim,lst->PIM,-delta,pk);
          xn -=delta;
          delta = delta/2.;
          imv++;
        } else {
          break;
        }
      }
      forwardError(*err,__LINE__,-1);      
    }

    fn = distribution_lkl(dist,lst->PIM,err);
    forwardError(*err,__LINE__,-1);
    
    delta = -delta/(yn-ynm) * yn;

    _DEBUGHERE_("%d -> lkl increase phi(%g) = %g -> phi(%g) = %g ",istep, xnm,fnm,xn,fn);
    
    
    if ((candic && fabs((xmoins-xplus)/xn) < 1e-8) || 
        (fabs(fnm-fn)<sc->tol) || 
        (fabs(delta)<1e-10)) {
      yn = deriv_along(dist,lst->PIM,pk,gk,1,err);
      forwardError(*err,__LINE__,-1);
      return xn;
    }

    ynm = yn;    
    
    fnm = fn;
        
    istep++;
  }
  testErrorRetVA(1==1,-111000111,"too many steps (%d) and braketing is still poor (%g)",*err,__LINE__,-1,sc->maxstep,fabs(delta/(xn+delta)));
}

optimize_struct* bydir_init(distribution *dist, double* guess,double *scale,linesearch_struct *ls, double tol, int maxstep, error **err) {
  bydir_optimize_struct *bd;
  optimize_struct *opt;
  
  bd = malloc_err(sizeof(bydir_optimize_struct),err);
  forwardError(*err,__LINE__,NULL);
  
  bd->pk = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  bd->scale = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  if (scale==NULL) {
    int i;
    for(i=0;i<dist->ndim;i++) {
      bd->scale[i] = 1;
    }
  } else{
    memcpy(bd->scale,scale,sizeof(double)*dist->ndim);
  }
  bd->PIM = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(bd->PIM,guess,sizeof(double)*dist->ndim);
  memset(bd->pk,0,sizeof(double)*dist->ndim);
  
  bd->tol = tol;  
  bd->ls = ls;
  bd->maxstep = maxstep;
  
  opt = optimize_init(dist,bd, bydir_func, bydir_get_p, NULL,bydir_free, err);
  forwardError(*err,__LINE__,NULL);
  
  return opt;
}

void bydir_func(optimize_struct *opt, error **err) {
  bydir_optimize_struct *bd;
  int id,istep;
  double alpha;
  double dd,d0,fnm,fn;
  distribution *target;
  int ndim;
  
  bd = opt->data;
  target = opt->target;
  ndim = opt->ndim;
  
  istep = 0;
  fnm = distribution_lkl(target, bd->PIM,err);
  forwardError(*err,__LINE__,);

  while(istep<bd->maxstep) {
    dd = 0;
    d0 = 0;
    for(id=0;id<ndim;id++) {
      d0 += bd->PIM[id];
      memset(bd->pk,0,sizeof(double)*ndim);

      bd->pk[id]=1.*bd->scale[id];
      alpha = linesearch(bd->ls, target, bd->PIM,bd->pk,NULL,err);
      forwardError(*err,__LINE__,);
      addalphap(ndim,bd->PIM,alpha,bd->pk);
      dd += alpha*alpha;
    }
    fn = distribution_lkl(target, bd->PIM,err);
    forwardError(*err,__LINE__,);

    {
      int i;
      char whr[10000];
      whr[0]='\0';
      for(i=0;i<ndim;i++) {
        sprintf(whr,"%s %g",whr,bd->PIM[i]);
      }
      _DEBUGHERE_("step %d residu %g %g | %g -> %g",istep,sqrt(dd/d0),sqrt(dd),fnm,fn);
      _DEBUGHERE_("pos [ %s ]\n\n\n",whr);      
    }    
    if (sqrt(dd)/sqrt(d0)<1e-8 || fabs(fnm-fn)<bd->tol) {
      return;
    }
    fnm=fn;
    istep++;
  }
  testErrorRetVA(1==1,-111000111,"too many steps (%d) and braketing is still poor (%g)",*err,__LINE__,,bd->maxstep,sqrt(dd)/sqrt(d0));
}

void bydir_free(void **pbd) {
  bydir_optimize_struct *bd;
  
  bd = *pbd;
  free(bd->PIM);
  free(bd->pk);
  free(bd->scale);
  free_linesearch(&(bd->ls));
  free(bd);
  *pbd = NULL;
}

double *bydir_get_p(optimize_struct *opt, error **err) {
  double* res;
  bydir_optimize_struct *bd;
  double alpha;
  int i,istep;
  
  bd = opt->data;
  
  res = malloc_err(sizeof(double)*opt->ndim,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(res,bd->PIM,sizeof(double)*opt->ndim);
  return res;
}

#ifdef _WITH_LAPACK_

optimize_struct* bfgs_init(distribution *dist, double* guess,double *scale,linesearch_struct *ls, double tol, int maxstep, error **err) {
  bfgs_optimize_struct *bd;
  optimize_struct *opt;
  
  bd = malloc_err(sizeof(bfgs_optimize_struct),err);
  forwardError(*err,__LINE__,NULL);
  
  bd->bk = malloc_err(sizeof(double)*dist->ndim*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  bd->bk_buf = malloc_err(sizeof(double)*dist->ndim*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  bd->gk = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  bd->sk = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  bd->yk = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);

  bd->pk = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  bd->scale = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  if (scale==NULL) {
    int i;
    for(i=0;i<dist->ndim;i++) {
      bd->scale[i] = 1;
    }
  } else {
    memcpy(bd->scale,scale,sizeof(double)*dist->ndim);
  }
  bd->PIM = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(bd->PIM,guess,sizeof(double)*dist->ndim);
  memset(bd->pk,0,sizeof(double)*dist->ndim);
  
  bd->tol = tol;  
  bd->ls = ls;
  bd->maxstep = maxstep;
  
  opt = optimize_init(dist,bd, bfgs_func, bfgs_get_p, NULL,bfgs_free, err);
  forwardError(*err,__LINE__,NULL);
  
  return opt;
}

void bfgs_free(void **pbd) {
  bfgs_optimize_struct *bd;
  
  bd = *pbd;
  free(bd->PIM);
  free(bd->pk);
  free(bd->scale);
  free(bd->bk);
  free(bd->bk_buf);
  free(bd->gk);
  free(bd->sk);
  free(bd->yk);
  free_linesearch(&(bd->ls));
  
  free(bd);
  *pbd = NULL;
}

double *bfgs_get_p(optimize_struct *opt, error **err) {
  double* res;
  bfgs_optimize_struct *bd;
  double alpha;
  int i,istep;
  

  bd = opt->data;
  res = malloc_err(sizeof(double)*opt->ndim,err);
  forwardError(*err,__LINE__,NULL);
  memcpy(res,bd->PIM,sizeof(double)*opt->ndim);
  return res;
}

void bfgs_func(optimize_struct *opt, error **err) {
  bfgs_optimize_struct *bfgs;
  int ndim,id;
  char uplo;
  int i,istep,info,one,imx;
  double alpha,d0,dd,yts,done,dzero,sbs,nrm,fnm,fn;
  distribution *target;
  
  one = 1;
  done = 1;
  dzero = 0;
  
  bfgs = opt->data;
  target = opt->target;
  ndim = opt->ndim;
    
  memset(bfgs->bk,0,sizeof(double)*ndim*ndim);
  for(i=0;i<ndim;i++) {
    bfgs->bk[i*ndim+i] = 1/(bfgs->scale[i]*bfgs->scale[i]);
  }
  // compute gradient
  for(i=0;i<ndim;i++) {
    bfgs->gk[i] = distribution_first_derivative(target,i,bfgs->PIM,err);
    forwardError(*err,__LINE__,);
  }
  istep = 0;
  uplo  = 'L';
  fnm = distribution_lkl(target, bfgs->PIM,err);
  forwardError(*err,__LINE__,);
  while(istep<bfgs->maxstep) {
    int j;
    // chol
    memcpy(bfgs->bk_buf,bfgs->bk,sizeof(double)*ndim*ndim);
    dpotrf(&uplo,&ndim,bfgs->bk_buf,&ndim,&info);
    testErrorRetVA(info!=0,-101001010,"Could not cholesky decompose using dpotrf (status %d)",*err,__LINE__,,info);
    
    // solve Bk sk = - gk 
    for(i=0;i<ndim;i++) {
      bfgs->sk[i] = - bfgs->gk[i];
      forwardError(*err,__LINE__,);
    }
    dpotrs( &uplo, &ndim, &one,bfgs->bk_buf,&ndim, bfgs->sk, &ndim, &info );
    testErrorRetVA(info!=0,-101001010,"Could not solve using dpotrs (status %d)",*err,__LINE__,,info);
    
    // pk = sk * scales
    
    nrm=0;
    imx=0;
    for(i=0;i<ndim;i++) {
      nrm += bfgs->sk[i]*bfgs->sk[i];
      if (fabs(bfgs->sk[i])>fabs(bfgs->sk[imx])) {
        imx=i;
      }
    }
    nrm = sqrt(nrm);
    
    for(i=0;i<ndim;i++) {
      bfgs->pk[i] = - bfgs->sk[i]*bfgs->scale[imx]/nrm;
      bfgs->yk[i] = bfgs->gk[i];
    }
    alpha = linesearch(bfgs->ls, target, bfgs->PIM,bfgs->pk, bfgs->yk,err);
    forwardError(*err,__LINE__,);
    _DEBUGHERE_("alpha %g",alpha);
    
    // check convergence
    d0 = 0;
    dd = 0;
    for(i=0;i<ndim;i++) {
      dd += alpha * alpha;
      d0 += bfgs->PIM[i]*bfgs->PIM[i]/bfgs->scale[i]/bfgs->scale[i];
    }

    // PIM -> PIM + alpha pk
    addalphap(ndim,bfgs->PIM,alpha,bfgs->pk);

    fn = distribution_lkl(target, bfgs->PIM,err);
    forwardError(*err,__LINE__,);
    
    {
      char whr[10000];
      whr[0]='\0';
      for(i=0;i<ndim;i++) {
        sprintf(whr,"%s %g",whr,bfgs->PIM[i]);
      }
      _DEBUGHERE_("step %d residu %g %g | %g -> %g ",istep,sqrt(dd/d0),sqrt(dd),fnm,fn);

      _DEBUGHERE_("pos [ %s ]\n\n\n",whr);      
    }
    
    
    if ((d0 !=0 && sqrt(dd/d0)< 1e-8) || (d0==0 && sqrt(dd) < 1e-8) || fabs(fnm-fn)<bfgs->tol) {
      return;
    }
    fnm = fn;

    // compute gradient and yk and ys
    yts = 0; 
    for(i=0;i<ndim;i++) {
      double grd;
      grd = bfgs->yk[i];
      bfgs->yk[i] +=  -bfgs->gk[i];
      bfgs->gk[i] = grd;
      yts += bfgs->yk[i] * bfgs->sk[i];
    }
    yts = 1./yts;
    
    // pk->bksk
    dsymv(&uplo, &ndim, &done, bfgs->bk, &ndim, bfgs->sk, &one, &dzero, bfgs->pk, &one);
    // compute skBksk = sk pk
    sbs = 0;
    for(i=0;i<ndim;i++) {
      sbs += bfgs->pk[i] * bfgs->sk[i];
    }
    sbs = 1./sbs;
    // bk -> bk + yyt/yts
    dsyr(&uplo, &ndim, &yts, bfgs->yk, &one, bfgs->bk, &ndim);
    
    
    // bk-> bk - bksk sktbk / sbs 
    sbs = -sbs;
    dsyr(&uplo, &ndim, &sbs, bfgs->pk, &one, bfgs->bk, &ndim);
    
    istep++;
  }

  fnm = fn;

  testErrorRetVA(1==1,-111000111,"too many steps (%d) and braketing is still poor (%g)",*err,__LINE__,,bfgs->maxstep,sqrt(dd)/sqrt(d0));
}
#endif
