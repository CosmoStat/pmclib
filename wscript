VERSION = '0.0'
APPNAME = "pmclib"

srcdir = '.'
blddir = 'build'

import os.path as osp
from Utils import pprint,exec_command
import sys

def add_lib(conf,prefix,include,libpath,libname, funcname="",headername="",libs = [], uselib=[],defines=[]):
  if include:
    include_prefix = include
  else:
    include_prefix = osp.join(prefix,"include")
  if libpath:
    libpath_prefix = libpath.split(":")
  else:
    libpath_prefix = osp.join(prefix,"lib")
  conf.env.append_value("CPPPATH",include_prefix)
  if libs == []:
    libs = [libname]
  if type(uselib)==type(""):
    uselib = [uselib]
  if type(funcname)==type(""):
    funcname=[funcname]
  if type(defines)==type(""):
    defines=[defines]

    
  conf.check_cc(lib=libs, libpath = libpath_prefix, rpath=libpath_prefix ,uselib_store=libname,mandatory=1,uselib=uselib,defines=defines)
  for fnc in funcname:
    conf.check_cc(
      errmsg="failed (check whether lib is compiled in 32 or 64bits)",
      function_name=fnc,header_name=headername,uselib=" ".join([libname]+uselib),mandatory=1)
      
def add_lib_option(libname,opt,default="/usr/local/lib"):
  opt.add_option("--%s_islocal"%libname,action="store_true",default=False,help="%s has been installed with install%s"%(libname,libname))
  opt.add_option("--%s_prefix"%libname,action="store",default=default,help="%s include/lib path prefix"%libname)
  #opt.add_option("--%s_prefix"%libname,action="store",default="./",help="%s include/lib path prefix"%libname)
  opt.add_option("--%s_include"%libname,action="store",default="",help="%s include path"%libname)
  opt.add_option("--%s_lib"%libname,action="store",default="",help="%s lib path"%libname)
  opt.add_option("--%s_link"%libname,action="store",default="",help="%s link line"%libname)
  
def libsfromlinkline(ll):
  return [lb.strip() for lb in ll.split("-l") if lb.strip()]
  
def set_options(opt):
  #add_lib_option("gsl",opt,default="/media/mkilbing/DATA/miniconda3/envs/cosmopmc")
  add_lib_option("gsl",opt,default="MK_TO_REPLACE")
  #add_lib_option("gsl",opt)
  add_lib_option("lua",opt)
  add_lib_option("hdf5",opt)
  add_lib_option("fftw3",opt)
  add_lib_option("lapack",opt)

  opt.add_option("--m32",action="store_true",default=False,help="compile & link in 32bits")
  opt.add_option("--m64",action="store_true",default=False,help="compile & link in 64bits")
  opt.add_option("--local",action="store_true",default=False,help="install in current directory")
  opt.add_option("--gcc",action="store_true",default=False,help="Do not test for icc and only use gcc")
  
def configure(conf):
  import Options
  if not Options.options.gcc:
    try:
      conf.check_tool('icc')
    except:
      pprint("PINK", "icc not found, defaulting to gcc")
      conf.check_tool('gcc')
  else:
    conf.check_tool('gcc')

  #32 bits
  if sys.platform.lower()=="darwin":
    mopt = ""
    if Options.options.m64:
      mopt += "-arch x86_64 "
    if Options.options.m32:    
      mopt += "-arch i386 "    
  else:
    mopt = "-m32"
    if Options.options.m64:
      mopt = "-m64"
    
  
  conf.env.mopt=mopt
  conf.env.append_value('CCFLAGS',mopt.split())
  conf.env.append_value('LINKFLAGS',mopt.split())
 
  #install where ?
  conf.env.mprefix=conf.env.PREFIX

  import os
  localpref = os.getcwd()
  if Options.options.local:
    conf.env.PREFIX=localpref
    conf.env.mprefix=localpref
  conf.env.install_examples = osp.join(conf.env.PREFIX,"share","pmclib","examples")

  #do something about -rpath
  #import Utils
  #arch = Utils.cmd_output("uname").lower().strip()
  #if arch=="linux":
  #  conf.env.RPATH_ST="-Wl,-R,%s"

  # rpath
  conf.env.append_value("RPATH",conf.env.PREFIX+"/lib")
  conf.env.append_value("RPATH",conf.env.PREFIX+"/share/pmclib/examples")

  conf.env.shsuffix = "so"
  if sys.platform.lower()=="darwin":
    conf.env.shsuffix = "dylib"
    
  conf.store()
  # configure gsl location
  if Options.options.gsl_islocal:
    gsl_prefix=localpref
    gsl_include=""
    gsl_lib=""
    gsl_link=[]
  else :
    gsl_prefix = Options.options.gsl_prefix
    gsl_include = Options.options.gsl_include
    gsl_lib = Options.options.gsl_lib
    gsl_link = libsfromlinkline(Options.options.gsl_link)
  try:
    libs=["gsl","gslcblas"]
    if gsl_link:
      libs = gsl_link
    
    add_lib(conf,gsl_prefix,gsl_include,gsl_lib,"gsl",
          ["gsl_ran_gaussian","gsl_blas_dgemm","gsl_linalg_cholesky_invert"],["gsl/gsl_randist.h","gsl/gsl_blas.h","gsl/gsl_linalg.h"],
          libs=libs)
  except Exception as e:
    pprint("RED","gsl not found")
    pprint("PINK", "check that gsl_prefix or gsl_lib and gsl_include command line options point toward your gsl install")
    pprint("PINK", "or check that gsl is compiled in %d bit (as you have specified for pmclib)"%{True:64}.get(Options.options.m64,32))
    pprint("PINK","alternatively, I can also install gsl for you, type './waf installgsl'")
    #pprint("ORANGE", "continuing for now...")
    raise e
  
  # dl
  conf.check_cc(lib="dl",mandatory=0,defines=["HAS_RTLD_DEFAULT"],fragment="#include <dlfcn.h> \nint main() {void* tt = RTLD_DEFAULT;}",msg="checking for RTLD_DEFAULT in dl",uselib_store="dl")
    
  #configure lua
  try: 
    if Options.options.lua_islocal:
      lua_prefix=localpref
      lua_include=""
      lua_lib=""
      lua_link=[]
    else:
      lua_prefix = Options.options.lua_prefix
      lua_include = Options.options.lua_include
      lua_lib = Options.options.lua_lib
      lua_link = libsfromlinkline(Options.options.lua_link)
    libs = ["lua"]
    if lua_link:
      libs = lua_link
    libs += ["m","dl"]
    add_lib(conf,lua_prefix,lua_include,lua_lib,"lua","lua_newstate","lua.h", libs=libs,defines="_WITH_RC_")
    conf.env.dolua = True
  except:
    conf.env.dolua = False
    pprint("PINK","lua not found, will skip building things requiring lua")
    pprint("PINK", "check that lua_prefix or lua_lib and lua_include command line options point toward your lua install")
    pprint("PINK","alternatively, I can also install lua for you, type './waf installlua'")

  #configure hdf5
  try: 
    if Options.options.hdf5_islocal:
      hdf5_prefix=localpref
      hdf5_include=""
      hdf5_lib=""
      hdf5_link=[]
    else:
      hdf5_prefix = Options.options.hdf5_prefix
      hdf5_include = Options.options.hdf5_include
      hdf5_lib = Options.options.hdf5_lib
      hdf5_link = libsfromlinkline(Options.options.hdf5_link)

    libs = ["hdf5","hdf5_hl"]
    if hdf5_link:
      libs = hdf5_link
      
    #add_lib(conf,hdf5_prefix,hdf5_include,hdf5_lib,"hdf5","H5Fcreate","hdf5.h", libs=libs,defines="_WITH_HDF5_")
    #conf.env.hdf5 = True
  except:
    conf.env.hdf5 = False
    pprint("PINK", "hdf5 not found, will skip building things requiring hdf5")
    pprint("PINK", "check that hdf5_prefix or hdf5_lib and hdf5_include command line options point toward your hdf5 install")
    pprint("PINK","alternatively, I can also install lua for you, type './waf installhdf5'")
  
  
  #configure fftw3
  try: 
    if Options.options.fftw3_islocal:
      fftw3_prefix=localpref
      fftw3_include=""
      fftw3_lib=""
      fftw3_link=[]
    else:
      fftw3_prefix = Options.options.fftw3_prefix
      fftw3_include = Options.options.fftw3_include
      fftw3_lib = Options.options.fftw3_lib
      fftw3_link = libsfromlinkline(Options.options.fftw3_link)

    libs = ["fftw3"]
    if fftw3_link:
      libs = fftw3_link

    add_lib(conf,fftw3_prefix,fftw3_include,fftw3_lib,"fftw3","fftw_execute","fftw3.h", libs=libs,defines="_WITH_FFTW3_")
    conf.env.fftw3 = True
  except Exception as e:
    print(e)
    conf.env.fftw3 = False
    pprint("PINK", "fftw3 not found, will skip building things requiring fftw3")
    pprint("PINK", "check that fftw3_prefix or fftw3_lib and fftw3_include command line options point toward your fftw3 install")
    pprint("PINK","alternatively, I can also install lua for you, type './waf installfftw3'")


  #configure lapack
  try: 
    lapack_prefix = Options.options.lapack_prefix
    lapack_include = Options.options.lapack_include
    lapack_lib = Options.options.lapack_lib
    lapack_link = libsfromlinkline(Options.options.lapack_link)

    libs = ["BLAS","LAPACK"]
    includes = ["lapack.h","blas.h"]
    conf.env.mkl = False
    extradefs = []

    if lapack_link:
      libs = lapack_link

    if "mkl" in Options.options.lapack_link:
      conf.env.mkl = True
      extradefs += ["_HAS_MKL_"]
      includes = ["mkl_lapack.h","mkl_blas.h"]

    if "vecLib" in Options.options.lapack_lib:
      #must be darwin !
      raise Exception("Accelerate framework is not supported at this time")
      includes = ["cblas.h","clapack.h"]
      extradefs += ["_HAS_APPLELAPACK_"]

    add_lib(conf,lapack_prefix,lapack_include,lapack_lib,"lapack","dpotrf",includes, libs=libs,defines=["_WITH_LAPACK_"]+extradefs)
    conf.env.lapack = True
    
  except Exception as e:
    conf.env.lapack = False
    pprint("PINK", "lapack not found, will skip building things requiring lapack")
    pprint("PINK", "check that lapack_prefix or lapack_lib and lapack_include command line options point toward your lapack install")
    pprint("PINK", "check that lapack_link contains the correct link options")
    pprint("PINK","error was : %s"%e)
    
  # add tools
  conf.env.append_value("CPPPATH","pmctools/include")
  
  
   
  # configure for mpi
  conf.env.dompi = False
  if conf.find_program("mpicc"):
    conf.env.dompi = True
    envmpi = conf.env.copy() 
    conf.set_env_name('mpi', envmpi) 
    conf.setenv('mpi')
    conf.env.CC = "mpicc"
    conf.env.LINK_CC = "mpicc"
    envmpibld = envmpi = conf.env.copy() 
    conf.set_env_name('mpibld', envmpibld) 
    conf.setenv('default')
  
  print("configure ok\n\nnow run './waf build install'")

def build(bld):
  print("build")
  bld.add_subdirs("pmctools")
  bld.add_subdirs("pmclib")  
  bld.add_subdirs("pmcexec")  
  bld.add_subdirs("examples")

def _installsmthg_pre(ctx,where,what):

  import Options, Environment,Utils
  import urllib2
  import re
  import os.path as osp
  import tarfile
  import shutil
  
  ctx.env = Environment.Environment(filename="build/c4che/default.cache.py")
  
  if osp.exists("build/"+what):
    pprint("PINK","%s already downloaded"%what)
  else:
    pprint("PINK","download from "+where)
    luaf = urllib2.urlopen(where)
    #if luaf.code!=200 and luaf.code!=None:
    #  raise Utils.WscriptError("Cannot install : %d reported error %d"%(luaf.code,where))
    f=open("build/"+what,"w")
    print >>f,luaf.read(),
    luaf.close()
    f.close()
  tf = tarfile.open("build/"+what)
  for ff in [ff.path for ff in tf.getmembers()[0:2] if not osp.split(ff.path)[0]]:
    if osp.exists("build/"+ff):
      if osp.isdir("build/"+ff):
        shutil.rmtree("build/"+ff)
      else:
        os.remove("what/"+ff)
  tf.close()
  pprint("PINK","untar "+what)
  if exec_command("cd build/;tar -zxf "+what)!=0:
    raise Utils.WscriptError("Cannot untar "+what)

def _installsmthg_post(ctx,where,what,extra_config=""):
  import Options, Environment,Utils
  CCMACRO = "\"gcc %s\""%ctx.env.mopt
  CCMACRO = "CC=%s CXX=%s "%(CCMACRO,CCMACRO)
  CPPMACRO = "CPP=\"gcc -E\" CXXCPP=\"g++ -E\" "
  cmdline = "cd build/%s; ./configure --prefix=%s %s  %s %s; make clean;make;make install"%(where,ctx.env.mprefix,extra_config,CCMACRO, CPPMACRO)
  pprint("PINK",cmdline)
  if exec_command(cmdline)!=0:
    raise Utils.WscriptError("Cannot build %s"%what)
  pprint("GREEN","You can now run ./waf configure, adding the following option '--%s_islocal'"%what)
  
def installhdf5(ctx):
  _installsmthg_pre(ctx,"http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.5-patch1.tar.gz","hdf5-1.8.5-patch1.tar.gz")
  _installsmthg_post(ctx,"hdf5-1.8.5-patch1","hdf5")

  
def installgsl(ctx):
  _installsmthg_pre(ctx,"ftp://ftp.gnu.org/gnu/gsl/gsl-1.14.tar.gz","gsl-1.14.tar.gz")
  _installsmthg_post(ctx,"gsl-1.14","gsl")


def installfftw(ctx):
  _installsmthg_pre(ctx,"http://www.fftw.org/fftw-3.2.2.tar.gz","fftw-3.2.2.tar.gz")
  _installsmthg_post(ctx,"fftw-3.2.2","fftw","--enable-shared")

def installlua(ctx):

  _installsmthg_pre(ctx,"http://www.lua.org/ftp/lua-5.1.4.tar.gz","lua-5.1.4.tar.gz")
  import Options, Environment
  import re
  import sys
  plat = "linux"
  if sys.platform.lower()=="darwin":
    plat = "macosx"
    
  pprint("PINK","patching lua-5.1.4/src/Makefile")
  f=open("build/lua-5.1.4/src/Makefile")
  mk=f.read()
  f=open("build/lua-5.1.4/src/Makefile","w")
  mk = re.sub("PLAT= none","PLAT= "+plat,mk)
  mk = re.sub("CFLAGS= -O2 -Wall ","CFLAGS= -O2 -Wall %s -fPIC "%ctx.env.mopt,mk)
  mk = re.sub("CC= gcc","CC= gcc %s "%ctx.env.mopt,mk)
  print >>f,mk
  f.close()
  pprint("PINK","patching lua-5.1.4/Makefile")
  f=open("build/lua-5.1.4/Makefile")
  mk=f.read()
  f=open("build/lua-5.1.4/Makefile","w")
  mk = re.sub("PLAT= none","PLAT= "+plat,mk)
  mk = re.sub("/usr/local","%s"%ctx.env.mprefix,mk)
  print >>f,mk
  f.close()
  pprint("PINK","patching lua-5.1.4/src/luaconf.h")
  f=open("build/lua-5.1.4/src/luaconf.h")
  mk=f.read()
  f=open("build/lua-5.1.4/src/luaconf.h","w")
  mk = re.sub("/usr/local","%s"%ctx.env.mprefix,mk)
  print >>f,mk
  f.close()
  
  _installsmthg_post(ctx,"lua-5.1.4","lua")
  


