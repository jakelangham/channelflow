#!/bin/bash
# Regression case decay1: exact decaying scalar perturbation, mode m = 1, n = 1, vs/kappa = 1, Ri = 0.
# tools/settlingprofile writes the base settling profile plus eps phi_m(y) cos(2 pi n z/Lz) e^{sigma t}
# (see its header for the solution), with zero velocity. Integrated to T = 10 and compared with the
# exact field at T. Checks the non-mean scalar modes, the no-flux Robin BC and the diffusion and
# settling terms against an analytic solution. Run from this directory.
tier=quick
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
settlingprofile=${SETTLINGPROFILE:-../../build/tools/settlingprofile}
Re=400
T=10
dt=0.02
vs=2.5e-3
kappa=2.5e-3
vsok=1
m=1
n=1
eps=0.1
Ri=0.0
Lx=5.51156605892946
Lz=2.51327412287183
outdir=${OUTDIR:-decay1_simulateflow/}

mkdir -p $outdir
for t in 0 $T; do
    name=$( [ $t = 0 ] && echo initial || echo exact )
    $settlingprofile -Nx 16 -Ny 35 -Nz 32 -Lx $Lx -Lz $Lz -vsok $vsok -kap $kappa \
        -eps $eps -m $m -n $n -t $t $outdir/$name.nc
done
$simulateflow -R $Re -vsett $vs -kap $kappa -ri $Ri \
    -T $T -dT $T -dt $dt -vdt false -ts sbdf3 \
    -o $outdir $outdir/initial.nc

# Compared by run_regression.sh with the exact field at T. With the default SMRK2 start-up, SBDF3 is
# third order: the scalar error is 3.3e-13 relative at dt = 0.02 (2.6e-12 at 0.04); with -is sbdf1 it
# was 2.8e-8 (second order). Tolerance 10x that (3e-12): a lower-order start-up or linear() fails it.
# The velocity must stay zero.
reference=exact
tol_strict=1e-12
tol_loose_velocity=1e-12
tol_loose_scalar=3e-12
