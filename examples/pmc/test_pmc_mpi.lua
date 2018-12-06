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
pb.min = {-10,-10}
pb.max = {10,10}

-- mix_mvdens proposal
mix_mvdens = {}
mix_mvdens.ndim = 2
mix_mvdens.ncomp = 9
mix_mvdens.draw_from = { mean= {0,0}, var={10,10}}
mix_mvdens.name = "mix_mvdens_mpi"
mix_mvdens.save = "proposal.mix_mvdens"

-- pmc run
pmc = {}
pmc.nsample = 30000
pmc.niter = 4
pmc.target = bentGauss
pmc.proposal = mix_mvdens
pmc.pb = pb
pmc.print_pc = 20
pmc.perpfile = "pmc_exploration.perp"
pmc.resfile = "pmc_exploration_%%d.pmc"
pmc.propfile = "pmc_exploration_%%d.mix_mvdens"


rc={}
rc.verbose = 1