#!/bin/bash
# Regression case growth1: long-time boundedness of the settling scalar at high vs/kappa. Coupled DNS
# from upo1 (24x35x40) at vs/kappa = 10, Ri = 0, dt = 0.02, to T = 400 (~4 min). Until the fix of
# 2026-11-07 the module's scalar drifted from the conjugate symmetry of a real field (kz = 0 modes kx
# and -kx), and that invisible part grew from t ~ 290 (chebyNorm 16 at t = 300, 1e3 by t = 320); the
# original code (e9b7541) stays bounded. Checks chebyNorm at t = 400 against the fixed code (which
# matches the original code to all printed digits up to t = 360), and its maximum over the run.
# Long-running: full tier only. Run from this directory.
tier=full
simulateflow=${SIMULATEFLOW:-../../build/programs/simulateflow}
outdir=${OUTDIR:-growth1_simulateflow/}

mkdir -p $outdir
$simulateflow -R 400 -vsett 0.025 -kap 2.5e-3 -ri 0 -T 400 -dT 10 -dt 0.02 -vdt false -is sbdf1 \
    -s 0 -sprof 0 -o $outdir upo1_24x35x40.nc > $outdir/simulateflow.log 2>&1

# label value reference tolerance (read by run_regression.sh: |value - reference| <= tolerance)
grep 'chebyNorm(u) ==' $outdir/simulateflow.log | awk '
    { n = $3; if (n > max) max = n; last = n }
    END {
        printf "final %.6g 11.6351 0.01\n", last
        printf "max %.6g 11.6351 0.01\n", max
    }' > $outdir/checks.txt
