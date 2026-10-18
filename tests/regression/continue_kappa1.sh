#!/bin/bash
# Regression case kappa1: continuation in the diffusivity kappa. Continues eq1 (Ri = 0, vs/kappa = 1)
# in kappa for two steps, then evaluates the last solution with a fresh findsoln at its kappa, which
# builds the no-flux BC (alpha = vs/kappa) and the diffusion from scratch. About 9 minutes. Run from this directory.
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
param=kappa        # continuation parameter
start=$kappa       # its starting value (also used by run_regression.sh to check it advanced)
extra=""           # extra search flags for continuesoln and the fresh findsoln
infile=$here/eq1_16x35x32.nc
symms=$here/s1.asc
outdir=${OUTDIR:-kappa1_continuesoln/}

mkdir -p $outdir && cd $outdir || exit 1
# continuesoln writes its initial-*/ and search-*/ directories to the working directory
$continuesoln -eqb $extra -cont $param -R $Re -vsett $vs -kap $kappa -ri $Ri -T $T \
    -symms $symms -is sbdf1 -es 1e-12 -Ng 2000 -ns $nsteps $infile > continuesoln.log 2>&1

# Check the last solution with a fresh findsoln at its parameter value (no Newton steps, residual
# only), so that everything parameter-dependent (e.g. the no-flux BC) is rebuilt from scratch
# last search directory (last column of MuD.asc, which has an extra ax column with -xrel), and the
# parameter at full precision from its mu.asc (MuD.asc rounds to 8 digits)
dir_last=$(awk '!/^ *%/ { d = $NF } END { print d }' MuD.asc)
mu_last=$(cat $dir_last/mu.asc)
case $param in
    vsett) vs=$mu_last ;;
    kappa) kappa=$mu_last ;;
    ri) Ri=$mu_last ;;
esac
sigma=""
[ -f $dir_last/sigmabest.asc ] && sigma="-sigma $dir_last/sigmabest.asc"
mkdir -p check
$findsoln -eqb $extra $sigma -R $Re -vsett $vs -kap $kappa -ri $Ri -T $T -symms $symms -is sbdf1 \
    -es 1 -Nn 0 -od check/ $dir_last/ubest.nc $dir_last/ubest.nc > check/findsoln.log 2>&1
echo "$mu_last" > check/mu_last

# Read by run_regression.sh: the fresh residual must be below this
tol_fresh=1e-11
