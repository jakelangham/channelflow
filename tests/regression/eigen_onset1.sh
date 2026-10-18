#!/bin/bash
# Regression case onset1: linear-stability onset of the laminar settling state against Langham & Hogg
# (2026), Appendix B: at Re = 400, Sc = 1, vs/kappa = 1, in the box Lx = 2 pi/1.14, Lz = 2 pi/2.5, the
# streamwise-invariant (m, n) = (0, 1) mode becomes unstable at Ri_b = -9.08e-4 (computed there
# independently with Chebfun). findeigenvals is run on the laminar state (exact settling profile, zero
# velocity perturbation) at Ri_b on either side; the onset is where the leading growth rate crosses
# zero, by linear interpolation. Checks the coupled linear operator, buoyancy, the scalar BCs and the
# linearisation/Arnoldi path. About 10 minutes. Run from this directory.
tier=full
findeigenvals=${FINDEIGENVALS:-../../build/programs/findeigenvals}
settlingprofile=${SETTLINGPROFILE:-../../build/tools/settlingprofile}
Re=400
vs=2.5e-3
kappa=2.5e-3
T=10
Ri_lo=-8.5e-4   # stable side
Ri_hi=-9.5e-4   # unstable side
outdir=${OUTDIR:-onset1_findeigenvals/}

mkdir -p $outdir
$settlingprofile -Nx 16 -Ny 35 -Nz 32 -Lx 5.51156605892946 -Lz 2.51327412287183 -vsok 1 $outdir/laminar.nc
for Ri in $Ri_lo $Ri_hi; do
    mkdir -p $outdir/ri$Ri
    $findeigenvals -R $Re -vsett $vs -kap $kappa -ri $Ri -T $T -dt 0.02 -vdt false -is sbdf1 -N 40 \
        -o $outdir/ri$Ri/ $outdir/laminar.nc > $outdir/ri$Ri/findeigenvals.log 2>&1
done
# leading growth rate (largest Re(1/T log Lambda)) at each Ri, then the interpolated onset
s_lo=$(awk 'NR == 2 { print $1 }' $outdir/ri$Ri_lo/lambda.asc)
s_hi=$(awk 'NR == 2 { print $1 }' $outdir/ri$Ri_hi/lambda.asc)
awk -v a=$Ri_lo -v b=$Ri_hi -v sa=$s_lo -v sb=$s_hi \
    'BEGIN { printf "%s %s %s %s %.6e\n", a, sa, b, sb, a - sa * (b - a) / (sb - sa) }' > $outdir/onset.txt

# Read by run_regression.sh
onset_ref=-9.08e-4
tol_rel=5e-3
