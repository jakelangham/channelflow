#!/bin/bash
# Regression case traj1: time-stepping check. Starts from the eq1 reference (a Ri = 0 equilibrium)
# but integrates at Ri = 0.0025, where it is not steady, for T = 20 with fixed dt. The field at
# T = 20 is compared with traj1_T20.nc. The tolerance comes from the time-stepping error, measured
# as the distance between runs at dt and dt/2 (see README). Run from this directory.
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
Re=400
T=20
dt=${DT:-0.02}
vs=2.5e-3
kappa=2.5e-3
Ri=0.0025
infile=eq1_16x35x32.nc
outdir=${OUTDIR:-traj1_simulateflow/}

mkdir -p $outdir
$simulateflow -R $Re -vs $vs -kap $kappa -ri $Ri \
    -T $T -dT $T -dt $dt -vdt false -ts sbdf3 -is sbdf1 \
    -o $outdir $infile

# Reference and tolerances used by run_regression.sh (relative L2 distances at T = 20):
#   strict: unchanged equations and discretisation
#   loose:  10x the time-stepping error at dt = 0.02, estimated from runs at dt = 0.02, 0.01, 0.005
#           (errors 3.9e-7 velocity, 6.2e-8 scalar; observed order 2 from the sbdf1 start-up)
reference=traj1_T20.nc
tol_strict=1e-10
tol_loose_velocity=4e-6
tol_loose_scalar=6e-7
