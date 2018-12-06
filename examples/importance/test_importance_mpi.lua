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

-- importance run
importance = {}
importance.nsample = 30000
importance.target = bentGauss
importance.proposal = mix_mvdens
importance.pb = pb
importance.print_pc = 20

importance.file = "exploration_mpi.pmc"

rc={}
rc.verbose = 1