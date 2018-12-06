-- init bentGauss target
bentGauss = {}
bentGauss.ndim = 4
bentGauss.bent = .1
bentGauss.sig1_square = 10
bentGauss.conditional = {}
bentGauss.conditional.id = {2,3}
bentGauss.conditional.value = {1,2}
bentGauss.name = "bentGauss"
bentGauss.libpath = "%(path)s"

-- parabox
pb = {}
pb.min = {-20,-20}
pb.max = {20,20}

-- adaptMH proposal
adaptMH = {}
adaptMH.ndim = 2
adaptMH.sig0 = {10,100}
adaptMH.nadapt = 1000
adaptMH.adaptifaccepted = 1
adaptMH.k_damp = .5
--adaptMH.c_enlarge = 1.
adaptMH.name = "adaptMH"

-- mcmc run
mcmc = {}
mcmc.nsample = 30000*4
mcmc.target = bentGauss
mcmc.kernel = adaptMH
mcmc.pb = pb
mcmc.nbatch = 1000
mcmc.file = "exploration.mcmc"
mcmc.pars0 = {0,0}


rc={}
rc.verbose = 1