#!/bin/bash
# Regression case decay2: exact decaying scalar perturbation, mode m = 2, n = 2, vs/kappa = 4, Ri = 0.
# tools/settlingprofile writes the base settling profile plus eps phi_m(y) cos(2 pi n z/Lz) e^{sigma t}
# (see its header for the solution), with zero velocity. Integrated to T = 10 and compared with the
# exact field at T. Checks the non-mean scalar modes, the no-flux Robin BC and the diffusion and
# settling terms against an analytic solution. Run from this directory.
tier=full
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
settlingprofile=${SETTLINGPROFILE:-../../build/tools/settlingprofile}
Re=400
T=10
dt=0.02
vs=0.01
kappa=2.5e-3
vsok=4
m=2
n=2
eps=0.1
Ri=0.0
Lx=5.51156605892946
Lz=2.51327412287183
outdir=${OUTDIR:-decay2_simulateflow/}

mkdir -p $outdir
for t in 0 $T; do
    name=$( [ $t = 0 ] && echo initial || echo exact )
    $settlingprofile -Nx 16 -Ny 35 -Nz 32 -Lx $Lx -Lz $Lz -vsok $vsok -kap $kappa \
        -eps $eps -m $m -n $n -t $t $outdir/$name.nc
done
$simulateflow -R $Re -vsett $vs -kap $kappa -ri $Ri \
    -T $T -dT $T -dt $dt -vdt false -ts sbdf3 \
    -o $outdir $outdir/initial.nc

# Compared by run_regression.sh with the exact field at T. With the default SMRK2 start-up the scalar
# error is 1.6e-10 relative at dt = 0.02 (third order; 2.2e-7 with -is sbdf1); tolerance 10x that.
# The velocity must stay zero.
reference=exact
tol_strict=1e-12
tol_loose_velocity=1e-12
tol_loose_scalar=2e-9
