#!/bin/bash
# Regression case straty1: equilibrium of plane Couette flow with a linear stable stratification along y
# (wall-normal gravity), Neumann walls for the density perturbation (its mean is pinned at 0 in the
# Newton vector), Ri = 0.001, Sc = 2, 16x35x32, s1-symmetric (s1 leaves y, and so the density
# perturbation, unchanged). Continued from eq1's velocity at Ri = 0. Reconverges straty1_16x35x32.nc.
# Run from this directory.
tier=quick
findsoln=${FINDSOLN:-../../build/programs/findsoln}
Re=400
T=20
Ri=0.001
Sc=2
newt_tol=1e-12
g0_ref=2.2e-16  # residual of the stored reference under the code that made it
infile=straty1_16x35x32.nc
velfile=straty1_16x35x32.nc
symms=s1.asc
outdir=${OUTDIR:-straty1_findsoln/}

mkdir -p $outdir
$findsoln -eqb -R $Re -T $T -ri $Ri -strat y -sc $Sc -sbc neum -es $newt_tol -edx 1e-9 -opt 'hookstep' \
    -symms $symms -od $outdir \
    -is 'sbdf1' -Ng 2000 $velfile $infile
