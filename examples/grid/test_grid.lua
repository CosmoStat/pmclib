-- init bentGauss lkl
bentGauss = {}
bentGauss.ndim = 4
bentGauss.bent = .1
bentGauss.sig1_square = 10
bentGauss.conditional = {}
bentGauss.conditional.id = {2,3}
bentGauss.conditional.value = {1,2}
bentGauss.name = "bentGauss"
bentGauss.libpath = "%(path)s"

grid.file = "exploration.grid"

-- init grid
grid = {}

grid.target = bentGauss
grid.islog = 0

grid.nbin = 200 --same number of bins in each direction
-- use
-- grid.nbin = {50,30} 
-- to set different number of bin in each direction

grid.usepb = 1 -- use min/max from pb
-- use
-- grid.usepb = 0
-- grid.min = {0,0}
-- grid.max = {10,10}
-- to set min/max with different values than pb

-- or use
-- grid.pos = somehdfgridfile
-- to compute a grid at the same location that another one

--or use
-- grid.pos = { {1,2,3,4,5}, {100,200,300}}
-- to give a vector of indices for each location

-- init pb
pb = {}
pb.min = {-10,-10}
pb.max = {10,10}

grid.pb = pb

rc.verbose = 1