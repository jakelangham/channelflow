#!/bin/bash
# Regression case stratgrowth1: linear instability of plane Couette flow with a stable stratification
# along z (Facchini et al. 2018, JFM 853, their reference case Re = 966, Fr = 0.82, i.e. Ri = 1/Fr^2,
# Sc = 7, kx = 0.96, kz = 5.16; our U = y is their U = -y with x -> -x). DNS in a one-wavelength box
# from small random noise, to T = 600 (~30 s each); the growth rate of L2Norm over t = 300..600 (still
# linear: L2Norm < 3e-4) must match the leading eigenvalue of an independent Chebyshev-collocation
# solver of the linearised equations (scratch/strat/linstab.py, N = 60 and 100 agree to 10 digits): 0.0086343687 with
# Dirichlet walls for the density perturbation, 0.0077769007 with Neumann walls (the DNS is converged in
# Ny and dt to 1e-8). Full tier. Run from this directory.
tier=full
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
randomfield=${RANDOMFIELD:-$(dirname $simulateflow)/../tools/randomfield}
outdir=${OUTDIR:-stratgrowth1_simulateflow/}

mkdir -p $outdir
$randomfield -Nx 16 -Ny 49 -Nz 16 -Lx 6.544984694978736 -Lz 1.2176715389552284 -m 1e-4 -sd 3 \
    $outdir/u0 > $outdir/randomfield.log 2>&1
for bc in diri neum; do
    $simulateflow -R 966 -ri 1.48720999405116 -strat z -sc 7 -sbc $bc -T 600 -dT 10 -dt 0.02 -vdt false \
        -s 0 -sprof 0 -o $outdir/$bc/ $outdir/u0.nc > $outdir/$bc.log 2>&1
done

# label value reference tolerance (read by run_regression.sh: |value - reference| <= tolerance)
rate() {
    awk '!/^ *#/ { if ($1 == 300) a = $2; if ($1 == 600) b = $2 } END { printf "%.8f", log(b / a) / 300 }' "$1"
}
{
    echo "rate_diri $(rate $outdir/diri/energy.asc) 0.0086343687 1e-7"
    echo "rate_neum $(rate $outdir/neum/energy.asc) 0.0077769007 1e-7"
} > $outdir/checks.txt
