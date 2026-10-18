#!/bin/bash
# Regression case error_nonlinearity: simulateflow with the settling scalar and "-nl conv" is not supported and
# must stop with an error rather than produce a wrong result. run_regression.sh requires a non-zero
# exit and the expected message. Run from this directory.
tier=quick
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
settlingprofile=${SETTLINGPROFILE:-../../build/tools/settlingprofile}
outdir=${OUTDIR:-error_nonlinearity/}

mkdir -p $outdir
$settlingprofile -Nx 16 -Ny 35 -Nz 32 -vsok 1 $outdir/initial.nc
$simulateflow -R 400 -vsett 2.5e-3 -kap 2.5e-3 -ri 0.0025 -nl conv \
    -T 0.1 -dT 0.1 -dt 0.02 -vdt false -is sbdf1 -o $outdir $outdir/initial.nc
status=$?

# Read by run_regression.sh
expect="support only the rotational nonlinearity"
exit $status
