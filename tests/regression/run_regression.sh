#!/bin/bash
# Run the regression cases against their reference solutions and report PASS/FAIL.
#
# Usage: ./run_regression.sh [-b BUILD_DIR] [-t TOL] [case ...]
#   -b BUILD_DIR  channelflow build to test (default: ../../build)
#   -t TOL        tolerance on relative L2 distances to the reference and on |ax - ax_ref|
#                 (default: 1e-8)
#   case          any of: eq1 eq3 (default: all)
#
# For each findsoln case the reference is re-run through Newton and checked for:
#   G0     initial residual L2Norm(G) of the reference under the code being tested.
#          ~1e-13 or below means the discretised equations are unchanged; a larger value is
#          reported (INFO) but is not a failure by itself, e.g. after a change of time
#          discretisation, which moves f^T but not the steady solution.
#   G      final residual, must be below the case's Newton tolerance.
#   dist   relative L2 distance of the converged velocity and scalar to the reference.
#   ax     converged shift, for travelling waves.
# Output goes to run/<case>/ (git-ignored). Run serially; each case uses 1 thread.

set -u
cd "$(dirname "$0")"

build=../../build
tol=1e-8
while getopts "b:t:" opt; do
    case $opt in
        b) build=$OPTARG ;;
        t) tol=$OPTARG ;;
        *) sed -n '2,9p' "$0"; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
cases=${*:-"eq1 eq3"}

build=$(cd "$build" && pwd) || exit 2
export FINDSOLN=$build/programs/findsoln
fieldcompare=$build/tools/fieldcompare
for exe in "$FINDSOLN" "$fieldcompare"; do
    [ -x "$exe" ] || { echo "missing executable: $exe"; exit 2; }
done
export OMP_NUM_THREADS=1

# le A B: true if A <= B (floating point)
le() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a <= b) }'; }

# column NAME FILE ROW: value of the named column (header row 1) in row ROW ("last" for last row)
column() {
    awk -v name="$1" -v row="$3" '
        NR == 1 { for (i = 1; i <= NF; ++i) { h = $i; sub(/^%-/, "", h); if (h == name) c = i } next }
        { v = $c } NR == row + 1 { print v; found = 1; exit }
        END { if (!found && row == "last") print v }' "$2"
}

nfail=0
report() {  # case check status detail
    printf "  %-6s %-10s %-5s %s\n" "$1" "$2" "$3" "$4"
    [ "$3" = FAIL ] && nfail=$((nfail + 1))
}

echo "build: $build"
echo "tolerance: $tol"
for c in $cases; do
    script=findsoln_$c.sh
    [ -f "$script" ] || { echo "unknown case: $c"; nfail=$((nfail + 1)); continue; }
    ref=$(sed -n 's/^infile=//p' "$script")
    newt_tol=$(sed -n 's/^newt_tol=//p' "$script")
    sigma=$(sed -n 's/^sigma=//p' "$script")
    out=run/$c/
    rm -rf "$out" && mkdir -p "$out"

    echo "== $c (reference $ref)"
    OUTDIR=$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/findsoln.log" 2>&1
    read -r wall rss < "$out/time"
    conv=$out/convergence.asc
    if [ ! -s "$conv" ] || [ ! -f "$out/ubest.nc" ]; then
        report "$c" findsoln FAIL "no convergence.asc/ubest.nc; see $out/findsoln.log"
        continue
    fi

    g0=$(column L2Norm\(Gx\) "$conv" 1)
    g=$(column L2Norm\(Gx\) "$conv" last)
    nsteps=$(($(wc -l < "$conv") - 2))
    if le "$g0" 1e-12; then report "$c" G0 ok "$g0 (equations unchanged)"
    else report "$c" G0 INFO "$g0 (f^T differs from the reference code)"; fi
    if le "$g" "$newt_tol"; then report "$c" G PASS "$g after $nsteps Newton steps"
    else report "$c" G FAIL "$g > $newt_tol after $nsteps Newton steps"; fi

    while read -r label d rel n nref; do
        [ "$label" = all ] && continue
        if le "$rel" "$tol"; then report "$c" "$label" PASS "rel dist $rel"
        else report "$c" "$label" FAIL "rel dist $rel > $tol (abs $d)"; fi
    done < <("$fieldcompare" "$out/ubest.nc" "$ref")

    if [ -n "$sigma" ]; then
        ax_ref=$(awk '{print $5}' "$sigma")
        ax=$(awk '{print $5}' "$out/sigmabest.asc")  # full precision, unlike convergence.asc
        dax=$(awk -v a="$ax" -v b="$ax_ref" 'BEGIN { d = a - b; print (d < 0 ? -d : d) }')
        if le "$dax" "$tol"; then report "$c" ax PASS "$ax (|dax| = $dax)"
        else report "$c" ax FAIL "$ax vs $ax_ref (|dax| = $dax)"; fi
    fi
    echo "  ($wall s, peak RSS $((rss / 1024)) MB)"
done

if [ $nfail -eq 0 ]; then echo "ALL PASSED"; else echo "$nfail FAILED"; fi
exit $((nfail > 0))
