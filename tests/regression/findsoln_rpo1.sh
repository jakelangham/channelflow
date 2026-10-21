#!/bin/bash
# Regression case rpo1: relative periodic orbit at Ri = 5.24e-4 (buoyant settling scalar), vs/kappa =
# 0.56, T ~ 88.09, 24x35x40, sztxz-symmetric, drifting in x (Ri > 0 breaks the y-reflecting
# symmetries). upo1 continued in Ri with continuesoln -orb -xrel (2026-10-21, see README). Reconverges
# rpo1_24x35x40.nc from its converged shift (rpo1_sigma.asc); a code change that leaves the physics
# unchanged should start at L2Norm(G) ~ 2e-11. About 1 min. Run from this directory.
tier=full
findsoln=${FINDSOLN:-../../build/programs/findsoln}
Re=400
T=88.094540922300993
vs=1.4e-3
kappa=2.5e-3
Ri=0.00052412185045274274
newt_tol=1e-10
g0_ref=1.9e-11  # residual of the stored reference under the code that made it
infile=rpo1_24x35x40.nc
velfile=rpo1_24x35x40.nc
symms=sztxz.asc
sigma=rpo1_sigma.asc
outdir=${OUTDIR:-rpo1_findsoln/}

mkdir -p $outdir
$findsoln -orb -xrel -R $Re -vsett $vs -kap $kappa -T $T \
    -ri $Ri -es $newt_tol -opt 'hookstep' \
    -symms $symms -sigma $sigma -od $outdir \
    -is 'sbdf1' -Ng 2000 $velfile $infile
