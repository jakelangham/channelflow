#!/bin/bash
# Regression case floquet1: Floquet multipliers of upo1 (Ri = 0, vs/kappa = 0.56, T ~ 87.8, in the
# sztxz subspace) with findeigenvals (Arnoldi, 40 Krylov vectors; one f^T each, ~15 min). Checks the
# linearisation about a periodic orbit in the settling module: the leading multiplier (real, negative,
# |Lambda| = 14.32), the number of unstable multipliers (5: the real one and two complex pairs, |Lambda|
# = 2.33 and 1.65) and the two marginal ones from time and x-translation (|Lambda - 1| < 2e-3; the
# finite-difference Jacobian puts them 5e-4 and 7e-5 from 1). Run from this directory.
tier=full
findeigenvals=${FINDEIGENVALS:-../../build/programs/findeigenvals}
outdir=${OUTDIR:-upo1_floquet/}

mkdir -p $outdir
$findeigenvals -R 400 -vsett 1.4e-3 -kap 2.5e-3 -ri 0 -T 87.84095659533736 -is sbdf1 -symms sztxz.asc \
    -N 40 -o $outdir/ upo1_24x35x40.nc > $outdir/findeigenvals.log 2>&1

# The last table of abs(Lambda), arg(Lambda) in the log
awk '/abs\(Lambda\)/ { n = 0; delete m; delete a; on = 1; next }
     on && NF == 2 { m[++n] = $1; a[n] = $2; next }
     on && NF != 2 { on = 0 }
     END {
         lead = m[1] * (a[1] > 3.14 ? -1 : 1)          # real negative multiplier
         nu = 0; d1 = 1e9; d2 = 1e9
         for (i = 1; i <= n; ++i) {
             if (m[i] > 1.002) ++nu
             d = m[i] - 1; if (d < 0) d = -d
             if (a[i] == 0 && d < d1) { d2 = d1; d1 = d } else if (a[i] == 0 && d < d2) d2 = d
         }
         # label value reference tolerance (read by run_regression.sh: |value - reference| <= tolerance)
         printf "leading %.10g -14.31902963 1.4e-5\n", lead
         printf "unstable %d 5 0\n", nu
         printf "marginal1 %.3e 0 2e-3\n", d1
         printf "marginal2 %.3e 0 2e-3\n", d2
     }' $outdir/findeigenvals.log > $outdir/checks.txt
