#!/bin/bash
# Regression case eq3: travelling wave in x at Ri = 0.0025 (buoyant settling scalar), 16x35x32,
# sztxz-symmetric. Reconverges eq3_16x35x32.nc from its converged shift (eq3_sigma.asc); a code
# change that leaves the physics unchanged should start at L2Norm(G) ~ 4e-15 and stay there.
# Run from this directory.
tier=quick
findsoln=${FINDSOLN:-../../build/programs/findsoln}
Re=400
T=20
vs=2.5e-3
kappa=2.5e-3
Ri=0.0025
newt_tol=1e-12
infile=eq3_16x35x32.nc
velfile=eq3_16x35x32.nc
symms=sztxz.asc
sigma=eq3_sigma.asc
outdir=${OUTDIR:-eq3_findsoln/}

mkdir -p $outdir
$findsoln -eqb -xrel -R $Re -vsett $vs -kap $kappa -T $T \
    -ri $Ri -es $newt_tol -opt 'hookstep' \
    -symms $symms -sigma $sigma -od $outdir \
    -is 'sbdf1' -Ng 2000 $velfile $infile
