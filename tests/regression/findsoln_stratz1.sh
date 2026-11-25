#!/bin/bash
# Regression case stratz1: equilibrium of plane Couette flow with a linear stable stratification along z
# (spanwise gravity, Facchini et al. 2018), Dirichlet walls for the density perturbation, Ri = 0.001,
# Sc = 1, 16x35x32, s1-symmetric (s1 reflects z, so it flips the density perturbation). Continued from
# eq1's velocity at Ri = 0. Reconverges stratz1_16x35x32.nc. Run from this directory.
tier=quick
findsoln=${FINDSOLN:-../../build/programs/findsoln}
Re=400
T=20
Ri=0.001
Sc=1
newt_tol=1e-12
g0_ref=2.0e-14  # residual of the stored reference under the code that made it
infile=stratz1_16x35x32.nc
velfile=stratz1_16x35x32.nc
symms=s1.asc
outdir=${OUTDIR:-stratz1_findsoln/}

mkdir -p $outdir
$findsoln -eqb -R $Re -T $T -ri $Ri -strat z -sc $Sc -sbc diri -es $newt_tol -edx 1e-9 -opt 'hookstep' \
    -symms $symms -od $outdir \
    -is 'sbdf1' -Ng 2000 $velfile $infile
