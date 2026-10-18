#!/bin/bash
# Regression case eq1: equilibrium at Ri = 0 (passive settling scalar), 16x35x32, s1-symmetric.
# Reconverges eq1_16x35x32.nc; a code change that leaves the physics unchanged should start at
# L2Norm(G) ~ 1e-13 and stay there. Run from this directory.
findsoln=${FINDSOLN:-../../build/programs/findsoln}
Re=400
T=20
vs=2.5e-3
kappa=2.5e-3
Ri=0.0
newt_tol=1e-12
infile=eq1_16x35x32.nc
velfile=eq1_16x35x32.nc
symms=s1.asc
outdir=${OUTDIR:-eq1_findsoln/}

mkdir -p $outdir
$findsoln -eqb -R $Re -vs $vs -kap $kappa -T $T \
    -ri $Ri -es $newt_tol -opt 'hookstep' \
    -symms $symms -od $outdir \
    -is 'sbdf1' -Ng 2000 $velfile $infile
