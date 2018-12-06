import numpy as nm
import numpy.ma as ma
import numpy.random as ra
import numpy.linalg as la

## first building histograms

def ctr_level(his,lvl,infinite = False):

  mis=his.flat[:]*1.
  try:
    msk=(mis.mask==False)
    mis=nm.array(nm.compress(msk,mis))
  except Exception,e:
    pass
  mis.sort()
  cis=nm.cumsum(mis[::-1])
  cis/=cis[-1]

  alvl=nm.searchsorted(cis,lvl)[::-1]
  #print alvl
  clist=[0]+[mis[-ii] for ii in alvl]+[nm.max(mis)]
  if not infinite:
    return clist[1:]
  return clist

def sigma_level(his,half=False,infinite=False):
  if half:
    lvl=nm.array((68.26,90.,95.4,99,99.7))/100.
  else:
    lvl=nm.array((68.26,95.4,99.7))/100.
  clist = ctr_level(his,lvl,Infinite)

  return clist

def _covar(pars,w):
  
  if nm.sum(w[1:]-w[:-1])==0:
    return 
def smooth_scott(smooth=1):
  return lambda p,w:scott_binning(p,w,3.49*smooth)
  
def scott_binning(pars,w,cst=3.49):
  mean = nm.sum(pars*w[:,nm.newaxis])
  epars = pars - mean
  uh = nm.array([nm.outer(ep,ep) for ep in epars])
  var=nm.sum(uh*w[:,nm.newaxis,nm.newaxis],0)
  nu = cst*nm.dot(la.cholesky(var),[pars.shape[0]**(-1/(2.+pars.shape[1]))]*pars.shape[1])
  min = nm.min(pars,0)
  max = nm.max(pars,0)
  nb = nm.floor((max-min)/nu)
  return [nm.linspace(m,M,b+1) for m,M,b in zip(min,max,nb)]
    
class histogram(object):
  def __init__(self,pars,w=None,select=None,limits=None,bins=20):
    if w==None:
      w = nm.ones(pars.shape[0])*1.
    if select==None:
      select = nm.arange(pars.shape[1])
    opars = pars[:,select]
    if limits!=None:
      opars = opars[(opars<=limits[0]) * (opars>=limits[1])]
      w = w[(opars<=limits[0]) * (opars>=limits[1])]
    
    w = w/nm.sum(w)
    obins=bins
    if callable(bins):
      obins = bins(opars,w)
    self.histo,edges = nm.histogramdd(opars, bins=obins, normed=True, weights=w) 
    #self.var,edges = nm.histogramdd(opars**2, bins=obins, normed=True, weights=w) 
    #self.var-=self.histo
    self.bins = ([(ee[1:]+ee[:-1])/2. for ee in edges])
    self.binlims = [
      [ee[:-1] for ee in edges],
      [ee[1:] for ee in edges]
      ]
    self.extent = nm.array([(ee[0],ee[-1]) for ee in edges]).flat[:]
    self.shape = self.histo.shape
    
  def ctr_level(self,lvl,infinite = False):
    return ctr_level(self.histo,lvl,infinite)
    
  def sigma_level(self,half=False,infinite=False):
    return sigma_level(self.histo,half,infinite)

  def argpeak(self):
    uu=self.histo.argmax()
    rr = []
    for d in self.shape[::-1]:
      rr=[uu%d]+rr
      uu = (uu-uu%d)/d
    return nm.array(rr)
    
  def peak(self):
    ap=self.argpeak()
    rr = nm.array([self.bins[i][ap[i]] for i in range(len(self.histo.shape))])
    rr.shape = (len(self.histo.shape),)
    return rr
  def get_level_1d(self,lvl):
    if len(self.shape)!=1:
      raise TypeError("available only for 1D histograms")
    vv = self.ctr_level((lvl,),False)[0]
    bns = nm.argwhere(self.histo>=vv).flat[:]
    return (vv,self.bins[0][bns[0]],self.bins[0][bns[-1]])
    
import h5py
class grid(object):
  def __init__(self,fromfile):
    fi =  h5py.File(fromfile)
    self.log_lkl = ma.masked_array(fi["grid"]["log_lkl"][:],mask = fi["grid"]["flag"][:]==False)
    self.pos = [fi["grid"].attrs["pos_%d"%(i+1)][:] for i in range(fi["grid"].attrs["ndim"][0])]
    self.extent = nm.array([(pp[0],pp[-1]) for pp in self.pos]).flat[:]
    self.lkl = nm.exp(self.log_lkl-nm.max(self.log_lkl))
  def pretty_contour(self,select=(0,1),loc=[],levels=(68.26,95.4,99.7),cmap = plt.cm.Reds,linecolor = "k",alpha=1):
    try:
      import pylab as plt
      extent = [x[0],x[-1],y[0],y[-1]]
      if cmap!=None:
        plt.imshow(vals.T,origin="lower",extent=extent,aspect="auto",interpolation="gaussian",cmap=cmap,zorder=1,alpha=alpha)
      lvl = nm.array(levels)/100.
      
      sls = [loc[i] for i in range(min(select))]+[slice(0,len(self.pos[min(select)]))] + [loc[i] for i in range(min(select),max(select)-1)] + [slice(0,len(self.pos[max(select)]))] + [loc[i] for i in range(max(select)-1,len(loc))]
      
      cs = plt.contour(self.pos[select[0]],self.pos[select[1]],self.lkl[sls].T,ctr_level(vals,lvl),colors=linecolor,zorder=5)
      plt.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvl[::-1]])), fontsize=6)
    except Exception,e:
      print "cannot plot"
      print e
  
class chain(object):
  def __init__(self,file=None,pars=None,w=None,flag=None,log_lkl=None,deds=None):
    
    if file:
      fi = h5py.File(file)
      if "mcmc" in fi:
        self.__initmcmc(fi)
        return
      else:
        self.__initpmc(fi)
        return 
        
    self._pars=pars*1.

    if w!=None:
      self._w=(w)/nm.sum(w)
      self.has_w=True
    else:
      self.has_w=False

    if flag!=None:
      self.flag=flag.copy()
      self.has_flag=True
    else:
      self.has_flag=False
      self.flag=nm.ones(len(self._pars))

    if log_lkl!=None:
      self._log_lkl=log_lkl*1.
      self.has_log_lkl=True
    else:
      self.has_log_lkl=False

    if deds!=None:
      self._deds=deds*1.
      self.has_ded=True
    else:
      self.has_ded=False

  def _getAcceptedParams(self):
    if self.has_flag:
      return nm.compress(self.flag,self._pars,0)
    return self._pars
  acceptedParams=property(_getAcceptedParams)
  pars=acceptedParams  
  def _getRefusedParams(self):
    if self.has_flag:
      return nm.compress(self.flag==0,self._pars,0)
    return nm.zeros((0,self._pars.shape[1]))
  refusedParams=property(_getRefusedParams)

  def _getAcceptedW(self):
    if not self.has_w:
      return nm.ones((nm.sum(self.flag),))/nm.sum(self.flag)
      raise Exception("no W")
    if self.has_flag:
      ww=nm.compress(self.flag,self._w,0)
      return ww/nm.sum(ww)
    return self._w
  acceptedW=property(_getAcceptedW)
  w=acceptedW  
  def _getRefusedW(self):
    if not self.has_w:
      return nm.ones((len(self),))
    if self.has_flag:
      return nm.compress(self.flag==0,self._w,0)
    return nm.array([])
  refusedW=property(_getRefusedW)

  def _getAcceptedlog_lkl(self):
    if not self.has_log_lkl:
      raise Exception("no log_lkl")
    if self.has_flag:
      return nm.compress(self.flag,self._log_lkl,0)
    return self._log_lkl
  acceptedlog_lkl=property(_getAcceptedlog_lkl)
  log_lkl=acceptedlog_lkl  
  def _getRefusedlog_lkl(self):
    if not self.has_log_lkl:
      raise Exception("no log_lkl")
    if self.has_flag:
      return nm.compress(self.flag==0,self._log_lkl,0)
    return nm.array([])
  refusedlog_lkl=property(_getRefusedlog_lkl)

  def __len__(self):
    return self.pars.shape[0]
  def _getShape(self):
    return self.pars.shape
  shape=property(_getShape)


  def __inithdf(self,fi,name):
    pars = fi[name]["pars"][:]
    deds=None
    if fi[name].attrs["nded"][0]>0:
      deds = fi[name]["deds"][:]
    log_lkl = fi[name]["log_lkl"][:]
    return pars,deds,log_lkl
    
  def __initmcmc(self,fi):
    pars,deds,log_lkl = self.__inithdf(fi,"mcmc")
    self.__init__(pars = pars, log_lkl = log_lkl, deds = deds)
    
  def __initpmc(self,fi):
    pars,deds,log_lkl = self.__inithdf(fi,"pmc")
    flag = fi["pmc"]["flag"][:]
    w = fi["pmc"]["weight"][:]
    self.__init__(pars = pars, log_lkl = log_lkl, deds = deds,w=w,flag=flag)

  def repars(self,lc):
    pars=self.pars
    if callable(lc):
      epars=lc(pars)
    else:
      epars=pars
    return epars

  def integrate(self,func,lc=None):
    epars=self.repars(lc)
    return nm.sum(func(epars)*self.w[:,nm.newaxis],axis=0)

  def mean(self,lc=None):
    return self.integrate(lambda x:x,lc)

  def moment(self,order,lc=None,getmean=False,pmean=None):
    if pmean!=None:
      men=pmean
    else:
      men=self.mean(lc)
    moment=self.integrate(lambda x:(x-men)**order,lc)
    if getmean:
      return moment,men
    return moment

  def correlation(self,lc=None):
    epars = self.repars(lc)
    mean = self.mean(lc)
    epars = epars - mean
    uh = nm.array([nm.outer(ep,ep) for ep in epars])
    return nm.sum(uh*self.w[:,nm.newaxis,nm.newaxis],0)

  def resample(self,n,method="clever",lc=None):
    cw=nm.cumsum(self.w)
    cw=cw/cw[-1]
    if method.lower()=="clever":
      alea=ra.uniform()
      wids=(nm.arange(1,n+1)*alea)%1
    else:
      wids=ra.sample(n)
    nids=nm.searchsorted(cw,wids)
    npars=nm.take(self.pars,nids,axis=0)
    if callable(lc):
      npars = lc(npars)
    deds = None
    if self.has_deds:
      ndeds = nm.take(self.deds,nids)
    
    return chain(pars=npars,log_lkl = nm.take(self.log_lkl,nids),deds=deds)

  def K(self,prop=None,psi=None):
    w1=self.w
    w1=w1/nm.sum(w1)
    if psi==None:
      wpsi=1./len(w1)
    else:
      if prop == None:
        rprop = self.log_lkl-nm.log(self.w)
      else:
        rprop = prop(self.pars)
      lwpsi=psi(self.pars)-rprop
      lwpsi-=max(lwpsi)-100
      wpsi=nm.exp(lwpsi)
      wpsi=wpsi/nm.sum(wpsi)
    return nm.sum(w1*nm.log(w1/wpsi))

  def perp(self,prop=None,psi=None):
    return nm.exp(-self.K(prop,psi))

  def ess(self):
    return 1./nm.sum(self.w**2)

  def plot_triangle(self,select=None,scales=(),legend=(),levels=(68.26,95.4,99.7),show_prop=True,fill=68.26,show_mean=True,show_peak=True,show_extra=None,add_legend=r"$=%(peak).3g^{+%(up).3g}_{-%(down).3g}$",aspect=(8,8),fig=None,tick_at_peak=False,show_peak_2d=True,bins=smooth_scott(1.8),usetex=True):
    import pylab as plt
    if usetex:
      plt.rc('text', usetex=True)
    plt.rc('xtick',labelsize='8')
    plt.rc('ytick',labelsize='8')
    lvl = nm.array(levels)/100.

    if fig:
      plt.figure(fig,aspect)
    else:
      plt.figure(figsize=aspect)

    plt.clf()

    if select==None:
      select = range(self.shape[1])
    if not scales:
      scales = nm.ones(self.shape[1])
    scales=nm.array(scales)

    if show_mean:
      if not isinstance(show_mean,(list,tuple)): 
        show_mean=(True,)*self.shape[1]
    else:
      if not isinstance(show_mean,(list,tuple)): 
        show_mean=(False,)*self.shape[1]

    if show_peak:
      if not isinstance(show_peak,(list,tuple)): 
        show_peak=(True,)*self.shape[1]
    else:
      if not isinstance(show_peak,(list,tuple)): 
        show_peak=(False,)*self.shape[1]


    n = len(select)
    mns = self.mean()*scales
    var = self.moment(2)*scales**2
    pmax = nm.max(self.pars,axis=0)*scales
    pmin = nm.min(self.pars,axis=0)*scales
    dlt = (pmax-pmin)
    h1ds = [histogram(self.pars,select=(isel,),w=self.w,bins=bins) for isel in select]
    pks = nm.array([h1d.peak()[0] for h1d in h1ds])*scales
    if fill:
      lr = [h1d.get_level_1d(fill/100.) for h1d in h1ds]
    if tick_at_peak:
      ticks = nm.array((pmin+dlt*.1,pks,pmax-dlt*.1)).T
    else:
      ticks = nm.array((pmin+dlt*.1,(pmax+pmin)/2.,pmax-dlt*.1)).T

    dcts = [{"peak":pks[i],"up":lr[i][-1]*scales[i]-pks[i],"down":-lr[i][1]*scales[i]+pks[i],"mean":mns[i],"std":nm.sqrt(var[i]),"sigma":var[i]} for i in range(self.shape[1])]

    if show_extra:
      extra = []
      for ext in show_extra:
        if ext==None:
          extra+=[2*pmax[ii]]
        else:
          extra+=[ext]
    else:
      extra = 2*pmax

    for ii in range(len(select)):
      isel = select[ii]
      h1d = h1ds[ii]
      if show_prop  and self.has_w:
        h1d_prop = histogram(self.pars,select=(isel,),bins=bins)
      pos = ticks[ii]
      pos = [float("%.3g"%pp) for pp in pos]
      ax = plt.subplot(n,n,ii*n+ii+1)
      if fill:
        plt.fill_between(h1d.bins[0]*scales[ii],0,h1d.histo,where=((h1d.bins[0]>=lr[ii][1])*(h1d.bins[0]<=lr[ii][2])),color="red",alpha=.2)
      if show_prop and  self.has_w:
        plt.plot(h1d_prop.bins[0]*scales[ii],h1d_prop.histo,color="#669933",lw=0.5)
      plt.plot(h1d.bins[0]*scales[ii],h1d.histo,"black")
      plt.xticks(pos)
      plt.yticks(())
      if show_mean[ii]:
        plt.axvline(mns[ii],c="grey")
      if show_peak[ii]:
        plt.axvline(pks[ii],c="red")
      if extra[ii]>pmin[ii] and extra[ii]<pmax[ii]:
        plt.axvline(extra[ii],c="blue")

      ax.set_xlim((pmin[ii],pmax[ii]))
      if legend:
        ax.set_title((legend[ii]+add_legend)%dcts[ii],size="small")
        #plt.text(0.95, 0.9,legend[ii],
        #  horizontalalignment='right',
        #  verticalalignment='top',transform = ax.transAxes,size="small")
      for jj in range(ii):
        plt.subplot(n,n,ii*n+jj+1)
        if jj==0:
          plt.yticks(pos)
        else:
          plt.yticks(pos,("","",""))
      for jj in range(ii+1,len(select)):
        axs = plt.subplot(n,n,ii+jj*n+1)
        jsel = select[jj]
        h2d = histogram(self.pars,select=(isel,jsel),w=self.w,bins=bins)
        if show_prop:
          h2d_prop = histogram(self.pars,select=(isel,jsel),bins=bins)
        extent = nm.array(h2d.extent)
        extent[:2]*=scales[ii]
        extent[2:]*=scales[jj]
        plt.imshow(h2d.histo.T,origin="lower",extent=extent,aspect="auto",interpolation="gaussian",cmap=plt.cm.Reds,zorder=1)
        if show_peak[ii]:
          plt.axvline(pks[ii],c="red",zorder=1)
        if show_peak[jj]:
          plt.axhline(pks[jj],c="red",zorder=1)
        if show_peak_2d:
          pk2d = h2d.peak()*(scales[ii],scales[jj])
          plt.axvline(pk2d[0],c="red",ls="--",zorder=1)
          plt.axhline(pk2d[1],c="red",ls="--",zorder=1)
          
        if show_prop and self.has_w:
          plt.contour(h2d_prop.bins[0]*scales[ii],h2d_prop.bins[1]*scales[jj],h2d_prop.histo.T,h2d.ctr_level(lvl),colors="#669933",linewidths=0.5,zorder=3)
        cs = plt.contour(h2d.bins[0]*scales[ii],h2d.bins[1]*scales[jj],h2d.histo.T,h2d.ctr_level(lvl),colors="k",zorder=5)
        plt.clabel(cs, cs.levels[:-1], inline=True,inline_spacing=0, fmt=dict(zip(cs.levels[:-1],[r"%d \%%"%int(l*100) for l in lvl[::-1]])), fontsize=6)

        axs.set_xlim((pmin[ii],pmax[ii]))
        axs.set_ylim((pmin[jj],pmax[jj]))

        if jj==n-1:
          plt.xticks(pos)
        else:
          plt.xticks(pos,("","",""))
