-- init bentGauss target
bentGauss = {}
bentGauss.ndim = 4
bentGauss.bent = .1
bentGauss.sig1_square = 100
bentGauss.conditional = {}
bentGauss.conditional.id = {2,3}
bentGauss.conditional.value = {1,2}
bentGauss.name = "bentGauss"
bentGauss.libpath = "%(path)s"



linesearch = {}
linesearch.name="secant_linesearch"
linesearch.ndim = 6
linesearch.tol = 1e-5
-- optimize
optimize = {}
optimize.guess = {1.,-1.}
optimize.name="bfgs"
optimize.tol = 1e-5
optimize.scales = {1,100}
-- optimize.scales = optimize.guess
rc={}
rc.verbose = 1