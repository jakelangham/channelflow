#!/bin/bash
# Regression case steady: the exact steady settling profile (vs/kappa = 1, Ri = 0.0025, zero
# velocity perturbation on plane Couette flow) under every time-stepping scheme and both mean
# constraints. Every combination must keep it steady (to rounding) for T = 2. Extends upstream's
# dnsParabolaTest/dnsZeroTest (exact steady flow under all schemes) to the settling scalar.
# run_regression.sh runs this once per combination with TS (scheme), IS (start-up scheme for
# multistep schemes), MC (gradp or bulkv), NL (nonlinearity form) and ROT (rotation rate) set.
# Run from this directory.
tier=full
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
settlingprofile=${SETTLINGPROFILE:-../../build/tools/settlingprofile}
Re=400
T=2
dt=0.02
vs=2.5e-3
kappa=2.5e-3
vsok=1
Ri=0.0025
Lx=5.51156605892946
Lz=2.51327412287183
outdir=${OUTDIR:-steady_schemes/}
ts=${TS:-sbdf3}
is=${IS:-sbdf1}
mc=${MC:-gradp}
nl=${NL:-rot}
rot=${ROT:-0}

mkdir -p $outdir
$settlingprofile -Nx 16 -Ny 35 -Nz 32 -Lx $Lx -Lz $Lz -vsok $vsok $outdir/initial.nc
$simulateflow -R $Re -vsett $vs -kap $kappa -ri $Ri -mc $mc -Ubulk 0 -nl $nl -rot $rot \
    -T $T -dT $T -dt $dt -vdt false -ts $ts -is $is \
    -o $outdir $outdir/initial.nc

# Combinations run by run_regression.sh, as scheme:startup:constraint[:nonlinearity[:rotation]]
# (defaults rot and 0). Laminar Couette flow with the settling profile stays exactly steady under
# every form of the nonlinear term, and under rotation about z (the Coriolis force on the base flow is
# balanced by the pressure).
combinations="cnfe1:sbdf1:gradp sbdf1:sbdf1:gradp sbdf2:sbdf1:gradp sbdf3:sbdf1:gradp sbdf4:sbdf1:gradp
              sbdf3:cnfe1:gradp sbdf3:smrk2:gradp sbdf3:cnrk2:gradp
              cnab2:sbdf1:gradp cnrk2:sbdf1:gradp smrk2:sbdf1:gradp
              sbdf3:sbdf1:bulkv cnrk2:sbdf1:bulkv
              sbdf3:sbdf1:gradp:conv sbdf3:sbdf1:gradp:skew sbdf3:sbdf1:gradp:div sbdf3:sbdf1:gradp:alt
              sbdf3:sbdf1:gradp:rot:0.5"
# Combinations known not to work yet (see TODO.md); reported as XFAIL, or XPASS once they work.
# (Before the settling module: CNAB2/CNRK2/SMRK2 and the RK start-ups ignored the scalar in
# NSE::linear, and -mc bulkv, non-rotational nonlinear forms and rotation stopped with an error.)
xfail=""
reference=initial
tol_strict=1e-12
