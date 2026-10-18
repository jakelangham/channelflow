#!/bin/bash
# Regression case vsett1: continuation in the settling velocity. Continues eq1 (Ri = 0, vs/kappa = 1)
# in vsett for two steps, then evaluates the last solution with a fresh findsoln at its vs, which
# builds the no-flux BC from scratch. If the BC followed vs during continuation, that residual is
# at the Newton-tolerance level (~1e-14); with the BC frozen at the starting vs it would not be
# (the same check at the previous step's vs gives ~1e-3). About 9 minutes. Run from this directory.
tier=full
here=$(pwd)
continuesoln=${CONTINUESOLN:-$here/../../build/programs/continuesoln}
findsoln=${FINDSOLN:-$here/../../build/programs/findsoln}
Re=400
T=20
vs=2.5e-3
kappa=2.5e-3
Ri=0.0
nsteps=2
infile=$here/eq1_16x35x32.nc
symms=$here/s1.asc
outdir=${OUTDIR:-vsett1_continuesoln/}

mkdir -p $outdir && cd $outdir || exit 1
# continuesoln writes its initial-*/ and search-*/ directories to the working directory
$continuesoln -eqb -cont vsett -R $Re -vsett $vs -kap $kappa -ri $Ri -T $T \
    -symms $symms -is sbdf1 -es 1e-12 -Ng 2000 -ns $nsteps $infile > continuesoln.log 2>&1

# Check the last solution with a fresh findsoln at its vs (no Newton steps, residual only)
last=$(awk '!/^ *%/ { d = $6; v = $1 } END { print v, d }' MuD.asc)
vs_last=${last% *}
dir_last=${last#* }
mkdir -p check
$findsoln -eqb -R $Re -vsett $vs_last -kap $kappa -ri $Ri -T $T -symms $symms -is sbdf1 \
    -es 1 -Nn 0 -od check/ $dir_last/ubest.nc $dir_last/ubest.nc > check/findsoln.log 2>&1
echo "$vs_last" > check/vs_last

# Read by run_regression.sh: the fresh residual must be below this
tol_fresh=1e-11
