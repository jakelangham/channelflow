#!/bin/bash
# Run the regression cases against their reference solutions and report PASS/FAIL.
#
# Usage: ./run_regression.sh [-a] [-b BUILD_DIR] [-m MPI_BUILD_DIR] [-t TOL] [case ...]
#   -a            run every tier (quick and full); default is the quick tier only
#   -b BUILD_DIR  channelflow build to test (default: ../../build)
#   -m MPI_BUILD  an MPI build of channelflow, for the "mpi" case (default: $MPI_BUILD if set)
#   -t TOL        tolerance on relative L2 distances to the reference and on |ax - ax_ref|
#                 (default: 1e-8)
#   case          run only these cases, whatever their tier: unit, or the name of any
#                 findsoln_<case>.sh / simulate_<case>.sh / continue_<case>.sh script here
#
# Tiers: each case script sets tier=quick (run after every change, ~1 min in total) or
# tier=full (minutes to an hour; before merges, restructure milestones and production runs).
# The pseudo-case "unit" (quick tier) runs the CTest unit/integration tests that currently pass.
# The pseudo-case "mpi" (full tier) re-runs eq1, eq3 and traj1 on 2 MPI processes (needs -m);
# results must match the serial references to the same tolerances.
#
# For each findsoln case the reference is re-run through Newton and checked for:
#   G0     initial residual L2Norm(G) of the reference under the code being tested.
#          Within 10x the case's g0_ref (the stored reference's own residual) means the
#          discretised equations are unchanged; a larger value is
#          reported (INFO) but is not a failure by itself, e.g. after a change of time
#          discretisation, which moves f^T but not the steady solution.
#   G      final residual, must be below the case's Newton tolerance.
#   dist   relative L2 distance of the converged velocity and scalar to the reference.
#   ax     converged shift, for travelling waves.
#   T      converged period, for periodic orbits (relative change).
# For each simulateflow (traj) case the final field is compared with the reference:
#   PASS (unchanged) if every relative distance is <= the case's tol_strict;
#   PASS (within time-stepping error) if <= its tol_loose_* (10x the dt error); FAIL otherwise.
#   The -t option does not apply to traj cases.
# For each continuesoln (continue_*) case the last continuation solution is re-evaluated with a
#   fresh findsoln at its parameter value; that residual must be below the case's tol_fresh.
# Every case also checks that its output field keeps the reference's file layout: NetCDF
# dimensions, variables (4 components, scalar as Component_3) and geometry attributes. Scripts
# outside channelflow read these files, so the layout is part of the interface.
# Output goes to run/<case>/ (git-ignored). Run serially; each case uses 1 thread.

set -u
cd "$(dirname "$0")"

build=../../build
tol=1e-8
all_tiers=0
mpibuild=${MPI_BUILD:-}
mpirun=${MPIRUN:-mpirun}
while getopts "ab:m:t:" opt; do
    case $opt in
        a) all_tiers=1 ;;
        m) mpibuild=$OPTARG ;;
        b) build=$OPTARG ;;
        t) tol=$OPTARG ;;
        *) sed -n '2,20p' "$0"; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
# CTest tests run by the "unit" case; extend as more of the legacy tests are fixed
unit_tests='^(unitTests|tridiag|cheby|laminar|helmholtz|tausolver|poisson|io|vector2field)$'

# Default case list: unit plus every case script whose tier is selected, in name order
if [ $# -gt 0 ]; then
    cases="$*"
else
    cases=unit
    [ $all_tiers = 1 ] && cases="$cases mpi"
    for f in findsoln_*.sh simulate_*.sh continue_*.sh; do
        c=${f#*_}; c=${c%.sh}
        t=$(sed -n 's/^tier=//p' "$f")
        { [ "$t" = quick ] || [ $all_tiers = 1 ]; } && cases="$cases $c"
    done
fi

build=$(cd "$build" && pwd) || exit 2
[ -n "$mpibuild" ] && { mpibuild=$(cd "$mpibuild" && pwd) || exit 2; }
export FINDSOLN=$build/programs/findsoln
export SIMULATEFLOW=$build/programs/simulateflow
export CONTINUESOLN=$build/programs/continuesoln
export SETTLINGPROFILE=$build/tools/settlingprofile
# Older builds may predate settlingprofile; it only writes the exact field, so use this tree's
[ -x "$SETTLINGPROFILE" ] || SETTLINGPROFILE=$(cd ../../build && pwd)/tools/settlingprofile
fieldcompare=$build/tools/fieldcompare
# Older builds may predate fieldcompare; it only measures, so use this tree's copy instead
[ -x "$fieldcompare" ] || fieldcompare=$(cd ../../build && pwd)/tools/fieldcompare
for exe in "$FINDSOLN" "$SIMULATEFLOW" "$fieldcompare"; do
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

# layout FILE: the NetCDF structure that external scripts rely on (no run-specific attributes)
layout() {
    ncdump -h "$1" | awk '
        /^dimensions:/ { on = 1 }
        /global attributes/ { on = 0; glob = 1; next }
        on { print; next }
        glob && /:(Nx|Ny|Nz|Lx|Lz|a|b|format_version|title|Conventions) =/ { print }'
}

nfail=0
report() {  # case check status detail
    printf "  %-6s %-10s %-5s %s\n" "$1" "$2" "$3" "$4"
    [ "$3" = FAIL ] && nfail=$((nfail + 1))
}

# check_format CASE FILE REF: the output field must keep the reference's file layout
check_format() {
    if ! command -v ncdump > /dev/null; then
        report "$1" format INFO "ncdump not found, layout not checked"
    elif diff <(layout "$2") <(layout "$3") > /dev/null; then
        report "$1" format PASS "same NetCDF layout as the reference"
    else
        report "$1" format FAIL "NetCDF layout differs from the reference: diff <(ncdump -h $2) <(ncdump -h $3)"
    fi
}

# run_cont CASE: continuation, then a fresh residual check of its last solution
run_cont() {
    local c=$1 script=continue_$1.sh out=run/$1/
    local tol vs0 vs_last g wall rss
    tol=$(sed -n 's/^tol_fresh=//p' "$script")
    vs0=$(sed -n 's/^vs=//p' "$script")
    rm -rf "$out" && mkdir -p "$out"
    echo "== $c (continuation)"
    OUTDIR=$PWD/$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/run.log" 2>&1
    read -r wall rss < "$out/time"
    if [ ! -s "$out/check/convergence.asc" ]; then
        report "$c" continue FAIL "no fresh residual; see $out/continuesoln.log and $out/check/"
        return
    fi
    vs_last=$(cat "$out/check/vs_last")
    if le "$vs_last" "$vs0"; then
        report "$c" continue FAIL "continuation did not advance (last vs = $vs_last); see $out/continuesoln.log"
    else
        report "$c" continue ok "reached vs = $vs_last from $vs0"
    fi
    g=$(column L2Norm\(Gx\) "$out/check/convergence.asc" 1)
    if le "$g" "$tol"; then report "$c" fresh PASS "residual $g at vs = $vs_last with a freshly built BC"
    else report "$c" fresh FAIL "residual $g > $tol at vs = $vs_last: BC did not follow the continuation?"; fi
    check_format "$c" "$out/check/../$(awk '!/^ *%/ { d = $6 } END { print d }' "$out/MuD.asc")/ubest.nc" eq1_16x35x32.nc
    echo "  ($wall s, peak RSS $((rss / 1024)) MB)"
}

# run_traj CASE: integrate a simulateflow case and compare its final field with the reference
run_traj() {
    local c=${2:-$1} script=simulate_$1.sh out=run/${2:-$1}/
    local ref T tol_strict wall rss label d rel n nref loose
    ref=$(sed -n 's/^reference=//p' "$script")
    T=$(sed -n 's/^T=//p' "$script")
    case $ref in initial | exact) ref=$out/$ref.nc ;; esac  # fields the case writes itself
    tol_strict=$(sed -n 's/^tol_strict=//p' "$script")
    rm -rf "$out" && mkdir -p "$out"
    echo "== $c (reference $ref)"
    OUTDIR=$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/simulateflow.log" 2>&1
    read -r wall rss < "$out/time"
    if [ ! -f "$out/u$T.nc" ]; then
        report "$c" simulate FAIL "no u$T.nc; see $out/simulateflow.log"
        return
    fi
    while read -r label d rel n nref; do
        [ "$label" = all ] && continue
        loose=$(sed -n "s/^tol_loose_$label=//p" "$script")
        if le "$rel" "$tol_strict"; then report "$c" "$label" PASS "rel dist $rel (unchanged)"
        elif le "$rel" "$loose"; then report "$c" "$label" PASS "rel dist $rel (within time-stepping error, <= $loose)"
        else report "$c" "$label" FAIL "rel dist $rel > $loose (abs $d)"; fi
    done < <("$fieldcompare" "$out/u$T.nc" "$ref")
    check_format "$c" "$out/u$T.nc" "$ref"
    echo "  ($wall s, peak RSS $((rss / 1024)) MB)"
}

# run_unit: the CTest subset, run in the build directory
run_unit() {
    local log=$PWD/run/unit.log
    mkdir -p run
    echo "== unit (ctest -R '$unit_tests')"
    if (cd "$build" && ctest -j1 -R "$unit_tests" > "$log" 2>&1); then
        report unit ctest PASS "$(grep -E 'tests passed' "$log")"
    else
        report unit ctest FAIL "$(grep -E 'tests passed' "$log"); see $log"
    fi
}

echo "build: $build"
echo "tolerance: $tol"
echo "cases: $cases"
# run_findsoln CASE [LABEL]: re-run a findsoln case from its reference; output in run/LABEL/
run_findsoln() {
    local c=$1 label=${2:-$1} script=findsoln_$1.sh out=run/${2:-$1}/
    local ref newt_tol sigma wall rss conv g0 g nsteps lab d rel n nref ax_ref ax dax
    ref=$(sed -n 's/^infile=//p' "$script")
    newt_tol=$(sed -n 's/^newt_tol=//p' "$script")
    sigma=$(sed -n 's/^sigma=//p' "$script")
    rm -rf "$out" && mkdir -p "$out"

    echo "== $label (reference $ref)"
    OUTDIR=$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/findsoln.log" 2>&1
    read -r wall rss < "$out/time"
    conv=$out/convergence.asc
    if [ ! -s "$conv" ] || [ ! -f "$out/ubest.nc" ]; then
        report "$label" findsoln FAIL "no convergence.asc/ubest.nc; see $out/findsoln.log"
        return
    fi

    g0=$(column L2Norm\(Gx\) "$conv" 1)
    g=$(column L2Norm\(Gx\) "$conv" last)
    nsteps=$(($(wc -l < "$conv") - 2))
    local g0_ref g0_lim
    g0_ref=$(sed -n 's/^g0_ref=\([^ ]*\).*/\1/p' "$script")
    g0_lim=$(awk -v r="${g0_ref:-1e-13}" 'BEGIN { print 10 * r }')
    if le "$g0" "$g0_lim"; then report "$label" G0 ok "$g0 (equations unchanged; reference ${g0_ref:-?})"
    else report "$label" G0 INFO "$g0 > 10x reference ${g0_ref:-?}: f^T differs from the reference code"; fi
    if le "$g" "$newt_tol"; then report "$label" G PASS "$g after $nsteps Newton steps"
    else report "$label" G FAIL "$g > $newt_tol after $nsteps Newton steps"; fi

    while read -r lab d rel n nref; do
        [ "$lab" = all ] && continue
        if le "$rel" "$tol"; then report "$label" "$lab" PASS "rel dist $rel"
        else report "$label" "$lab" FAIL "rel dist $rel > $tol (abs $d)"; fi
    done < <("$fieldcompare" "$out/ubest.nc" "$ref")
    check_format "$label" "$out/ubest.nc" "$ref"

    if grep -q -- '-orb' "$script" && [ -f "$out/Tbest.asc" ]; then
        local T_ref T_new dT
        T_ref=$(sed -n 's/^T=//p' "$script")
        T_new=$(cat "$out/Tbest.asc")
        dT=$(awk -v a="$T_new" -v b="$T_ref" 'BEGIN { d = (a - b) / b; print (d < 0 ? -d : d) }')
        if le "$dT" "$tol"; then report "$label" T PASS "$T_new (relative change $dT)"
        else report "$label" T FAIL "$T_new vs $T_ref (relative change $dT)"; fi
    fi

    if [ -n "$sigma" ]; then
        ax_ref=$(awk '{print $5}' "$sigma")
        ax=$(awk '{print $5}' "$out/sigmabest.asc")  # full precision, unlike convergence.asc
        dax=$(awk -v a="$ax" -v b="$ax_ref" 'BEGIN { d = a - b; print (d < 0 ? -d : d) }')
        if le "$dax" "$tol"; then report "$label" ax PASS "$ax (|dax| = $dax)"
        else report "$label" ax FAIL "$ax vs $ax_ref (|dax| = $dax)"; fi
    fi
    echo "  ($wall s, peak RSS $((rss / 1024)) MB)"
}

# run_mpi: re-run eq1, eq3 and traj1 on 2 MPI processes against the same references
run_mpi() {
    if [ -z "$mpibuild" ]; then
        report mpi build INFO "no MPI build given (-m DIR), skipped"
        return
    fi
    local fs=$FINDSOLN sf=$SIMULATEFLOW
    export FINDSOLN="$mpirun -np 2 $mpibuild/programs/findsoln"
    export SIMULATEFLOW="$mpirun -np 2 $mpibuild/programs/simulateflow"
    run_findsoln eq1 mpi_eq1
    run_findsoln eq3 mpi_eq3
    run_traj traj1 mpi_traj1
    export FINDSOLN=$fs SIMULATEFLOW=$sf
}

for c in $cases; do
    if [ "$c" = unit ]; then run_unit; continue; fi
    if [ "$c" = mpi ]; then run_mpi; continue; fi
    if [ -f simulate_$c.sh ]; then run_traj "$c"; continue; fi
    if [ -f continue_$c.sh ]; then run_cont "$c"; continue; fi
    if [ -f findsoln_$c.sh ]; then run_findsoln "$c"; continue; fi
    echo "unknown case: $c"; nfail=$((nfail + 1))
done

if [ $nfail -eq 0 ]; then echo "ALL PASSED"; else echo "$nfail FAILED"; fi
exit $((nfail > 0))
