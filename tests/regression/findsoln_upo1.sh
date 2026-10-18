#!/bin/bash
# Regression case upo1: periodic orbit at Ri = 0, vs/kappa = 0.56, period T ~ 87.9, 32x35x48,
# sztxz-symmetric, box Lx = 1.75 pi, Lz = 1.2 pi. Reconverges upo1_32x35x48.nc; a code change
# that leaves the physics unchanged should start at L2Norm(G) ~ 2e-11. About 1-2 minutes per
# f^T evaluation. Run from this directory.
tier=full
findsoln=${FINDSOLN:-../../build/programs/findsoln}
Re=400
T=87.919012419852208
vs=1.4e-3
kappa=2.5e-3
Ri=0.0
newt_tol=1e-10
g0_ref=2.0e-11  # residual of the stored reference under the code that made it
infile=upo1_32x35x48.nc
velfile=upo1_32x35x48.nc
symms=sztxz.asc
outdir=${OUTDIR:-upo1_findsoln/}

mkdir -p $outdir
$findsoln -orb -R $Re -vsett $vs -kap $kappa -T $T \
    -ri $Ri -es $newt_tol -opt 'hookstep' \
    -symms $symms -od $outdir \
    -is 'sbdf1' -Ng 2000 $velfile $infile
