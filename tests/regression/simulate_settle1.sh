#!/bin/bash
# Regression case settle1: exact steady settling profile, vs/kappa = 1, Ri = 0.0025.
# Writes C(y) = A exp(-(vs/kappa) y) with mean(C) = 1 and zero velocity perturbation
# (tools/settlingprofile), integrates it for T = 10 and requires it to stay put: the profile is an
# exact steady state (advection by the base flow vanishes, buoyancy is balanced hydrostatically).
# Checks the scalar operator, the no-flux Robin BC, the mass constraint and the hydrostatic
# pressure against the analytic solution rather than stored output. Run from this directory.
tier=quick
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
settlingprofile=${SETTLINGPROFILE:-../../build/tools/settlingprofile}
Re=400
T=10
dt=0.02
vs=2.5e-3
kappa=2.5e-3
vsok=1
Ri=0.0025
Lx=5.51156605892946
Lz=2.51327412287183
outdir=${OUTDIR:-settle1_simulateflow/}

mkdir -p $outdir
$settlingprofile -Nx 16 -Ny 35 -Nz 32 -Lx $Lx -Lz $Lz -vsok $vsok $outdir/initial.nc
$simulateflow -R $Re -vsett $vs -kap $kappa -ri $Ri \
    -T $T -dT $T -dt $dt -vdt false -ts sbdf3 -is sbdf1 \
    -o $outdir $outdir/initial.nc

# Compared by run_regression.sh with the initial (exact) field; it must not move.
reference=initial
tol_strict=1e-12
tol_loose_velocity=1e-12
tol_loose_scalar=1e-12
