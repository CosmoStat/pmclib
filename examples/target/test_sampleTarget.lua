-- init bentGauss lkl
bentGauss = {}
bentGauss.ndim = 4
bentGauss.bent = .001
bentGauss.sig1_square = 100
bentGauss.conditional = {}
bentGauss.conditional.id = {"extra_1","extra_2"}
bentGauss.conditional.value = {1,2}
bentGauss.name = "bentGauss"
bentGauss.libpath = "%(path)s"

bentGauss2 = {}
bentGauss2.ndim = 4
bentGauss2.bent = .1
bentGauss2.sig1_square = 1000
bentGauss2.conditional = {}
bentGauss2.conditional.id = {2,3}
bentGauss2.conditional.value = {1,2}
bentGauss2.name = "bentGauss"
bentGauss2.libpath = "%(path)s"

with_prior = {}
with_prior.distribution = bentGauss
with_prior.mean = {1,1}
with_prior.var = {1,10}
with_prior.prior_idx = {1,0}
with_prior.name = "add_gaussian_prior"

combine = {}
combine.name = "combine_distributions"
combine.ndim = 2
combine.distributions = {}
combine.distributions[1] = {}
combine.distributions[1].distribution = bentGauss
combine.distributions[2] = {}
combine.distributions[2].distribution = bentGauss
combine.distributions[2].dim_idx = {1,0}

-- test parameters
testbed = {}
testbed.pars = {4,1}
testbed.result = -2.99953
testbed.target = combine

-- verbosity 
--rc={}
rc.verbose = -1