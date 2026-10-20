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
#                 findsoln_<case>.sh / simulate_<case>.sh / continue_<case>.sh /
#                 schemes_<case>.sh / error_<case>.sh / eigen_<case>.sh script here
#
# Tiers: each case script sets tier=quick (run after every change, ~1 min in total) or
# tier=full (minutes to an hour; before merges, restructure milestones and production runs).
# The pseudo-case "unit" (quick tier) runs the CTest unit/integration tests that currently pass.
# The pseudo-case "mpi" (full tier) re-runs eq1, eq3 and traj1 on 2 MPI processes, with process
# grids 2x1 and 1x2 (needs -m);
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
# A schemes_* case runs one simulateflow script for each scheme:startup:constraint[:nl[:rot]] combination it
#   lists, requiring the final field to equal its reference to tol_strict. Combinations listed in
#   its xfail are known not to work yet: they report XFAIL (not a failure) while broken, and XPASS
#   (also not a failure, but a prompt to update xfail) once they start working.
# An eigen_* case writes onset.txt ("Ri_lo sigma_lo Ri_hi sigma_hi onset"); the leading growth rate
#   must change sign across the bracket and the onset must match onset_ref to tol_rel.
# An error_* case runs an unsupported configuration and requires a non-zero exit with the
#   message given by its expect=.
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
# every CTest test (upstream's, including the DNS tests, and the settling module's): ~20 s
unit_tests='.'

# Default case list: unit plus every case script whose tier is selected, in name order
if [ $# -gt 0 ]; then
    cases="$*"
else
    cases=unit
    [ $all_tiers = 1 ] && cases="$cases mpi"
    for f in findsoln_*.sh simulate_*.sh continue_*.sh schemes_*.sh error_*.sh eigen_*.sh floquet_*.sh; do
        [ -f "$f" ] || continue  # a pattern with no matching scripts stays unexpanded
        c=${f#*_}; c=${c%.sh}
        t=$(sed -n 's/^tier=//p' "$f")
        { [ "$t" = quick ] || [ $all_tiers = 1 ]; } && cases="$cases $c"
    done
fi

build=$(cd "$build" && pwd) || exit 2
[ -n "$mpibuild" ] && { mpibuild=$(cd "$mpibuild" && pwd) || exit 2; }
# RR_FINDSOLN etc. override individual programs (e.g. to test a module's versions)
export FINDSOLN=${RR_FINDSOLN:-$build/programs/findsoln}
export SIMULATEFLOW=${RR_SIMULATEFLOW:-$build/programs/simulateflow}
export CONTINUESOLN=${RR_CONTINUESOLN:-$build/programs/continuesoln}
export FINDEIGENVALS=${RR_FINDEIGENVALS:-$build/programs/findeigenvals}
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

# run_eigen CASE: stability onset from the leading growth rates at two Ri, against a reference
run_eigen() {
    local c=$1 script=eigen_$1.sh out=run/$1/ ref tol wall rss a sa b sb onset rel
    ref=$(sed -n 's/^onset_ref=//p' "$script")
    tol=$(sed -n 's/^tol_rel=//p' "$script")
    rm -rf "$out" && mkdir -p "$out"
    echo "== $c (stability onset)"
    OUTDIR=$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/run.log" 2>&1
    read -r wall rss < "$out/time"
    if ! read -r a sa b sb onset < "$out/onset.txt" 2> /dev/null || [ -z "$onset" ]; then
        report "$c" onset FAIL "no onset.txt; see $out"
        return
    fi
    if le "$sa" 0 && le 0 "$sb"; then report "$c" bracket PASS "sigma($a) = $sa < 0 < sigma($b) = $sb"
    else report "$c" bracket FAIL "no sign change: sigma($a) = $sa, sigma($b) = $sb"; fi
    rel=$(awk -v x="$onset" -v r="$ref" 'BEGIN { d = (x - r) / r; print (d < 0 ? -d : d) }')
    if le "$rel" "$tol"; then report "$c" onset PASS "Ri_b = $onset vs reference $ref (relative $rel)"
    else report "$c" onset FAIL "Ri_b = $onset vs reference $ref (relative $rel > $tol)"; fi
    echo "  ($wall s, peak RSS $((rss / 1024)) MB)"
}

# run_floquet CASE: the script writes checks.txt, lines "label value reference tolerance"; each must
# satisfy |value - reference| <= tolerance
run_floquet() {
    local c=$1 script=floquet_$1.sh out=run/floquet_$1/ wall rss label value ref tol n=0
    rm -rf "$out" && mkdir -p "$out"
    echo "== $c (Floquet multipliers)"
    OUTDIR=$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/run.log" 2>&1
    read -r wall rss < "$out/time"
    while read -r label value ref tol; do
        n=$((n + 1))
        if awk -v v="$value" -v r="$ref" -v t="$tol" 'BEGIN { d = v - r; if (d < 0) d = -d; exit !(d <= t) }'; then
            report "$c" "$label" PASS "$value (reference $ref, tolerance $tol)"
        else
            report "$c" "$label" FAIL "$value (reference $ref, tolerance $tol)"
        fi
    done < "$out/checks.txt"
    [ $n -gt 0 ] || report "$c" checks FAIL "no checks.txt; see $out"
    echo "  ($wall s, peak RSS $((rss / 1024)) MB)"
}

# run_error CASE: an unsupported configuration must stop with the expected message
run_error() {
    local c=$1 script=error_$1.sh out=run/error_$1/ expect status
    expect=$(sed -n 's/^expect="\(.*\)"/\1/p' "$script")
    rm -rf "$out" && mkdir -p "$out"
    echo "== $c (must stop with an error)"
    OUTDIR=$out ./"$script" > "$out/run.log" 2>&1
    status=$?
    if [ $status -ne 0 ] && grep -q -- "$expect" "$out/run.log"; then
        report "$c" error PASS "stopped (exit $status): $expect"
    elif [ $status -eq 0 ]; then
        report "$c" error FAIL "ran to completion instead of stopping; see $out/run.log"
    else
        report "$c" error FAIL "stopped (exit $status) without the expected message; see $out/run.log"
    fi
}

# run_schemes CASE: one run per scheme:startup:constraint combination, with expected failures
run_schemes() {
    local c=$1 script=schemes_$1.sh
    local combos xfail T tol combo ts is mc label out ref ok worst lab d rel n nref wall rss
    combos=$(sed -n '/^combinations="/,/"/p' "$script" | tr -d '"' | sed 's/^combinations=//')
    # the xfail="..." assignment, which may span lines (or be empty, xfail="")
    xfail=$(awk '/^xfail="/ { on = 1; sub(/^xfail="/, "") } on { if (sub(/".*/, "")) { print; exit } print }' \
        "$script" | tr -s ' \n' '  ')
    T=$(sed -n 's/^T=//p' "$script")
    tol=$(sed -n 's/^tol_strict=//p' "$script")
    echo "== $c (scheme:startup:constraint)"
    for combo in $combos; do
        IFS=: read -r ts is mc nl rot <<< "$combo"
        label=$c/$combo
        out=run/$c/${combo//:/_}/
        rm -rf "$out" && mkdir -p "$out"
        TS=$ts IS=$is MC=$mc NL=${nl:-rot} ROT=${rot:-0} OUTDIR=$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/simulateflow.log" 2>&1
        ok=1
        worst=""
        if [ ! -f "$out/u$T.nc" ]; then
            ok=0
            worst="no u$T.nc: $(grep -m1 -iE 'error|not supported' "$out/simulateflow.log" | cut -c1-90)"
        else
            ref=$out/initial.nc
            while read -r lab d rel n nref; do
                [ "$lab" = all ] && continue
                le "$rel" "$tol" || { ok=0; worst="$worst $lab $rel"; }
                [ -z "$worst" ] && worst="max rel dist $rel"
            done < <("$fieldcompare" "$out/u$T.nc" "$ref")
        fi
        case " $xfail " in
            *" $combo "*) if [ $ok = 1 ]; then report "$c" "$combo" XPASS "now works: remove it from xfail ($worst)"
                          else report "$c" "$combo" XFAIL "$worst"; fi ;;
            *)            if [ $ok = 1 ]; then report "$c" "$combo" PASS "$worst"
                          else report "$c" "$combo" FAIL "$worst; see $out"; fi ;;
        esac
    done
}

# run_cont CASE: continuation, then a fresh residual check of its last solution
run_cont() {
    local c=$1 script=continue_$1.sh out=run/$1/
    local tol param mu0 mu_last g wall rss
    tol=$(sed -n 's/^tol_fresh=//p' "$script")
    param=$(sed -n 's/^param=\([^ ]*\).*/\1/p' "$script")
    case $param in
        vsett) mu0=$(sed -n 's/^vs=//p' "$script") ;;
        kappa) mu0=$(sed -n 's/^kappa=//p' "$script") ;;
        ri) mu0=$(sed -n 's/^Ri=//p' "$script") ;;
    esac
    rm -rf "$out" && mkdir -p "$out"
    echo "== $c (continuation in $param)"
    OUTDIR=$PWD/$out /usr/bin/time -f "%e %M" -o "$out/time" ./"$script" > "$out/run.log" 2>&1
    read -r wall rss < "$out/time"
    if [ ! -s "$out/check/convergence.asc" ]; then
        report "$c" continue FAIL "no fresh residual; see $out/continuesoln.log and $out/check/"
        return
    fi
    mu_last=$(cat "$out/check/mu_last")
    if le "$mu_last" "$mu0"; then
        report "$c" continue FAIL "continuation did not advance (last $param = $mu_last); see $out/continuesoln.log"
    else
        report "$c" continue ok "reached $param = $mu_last from $mu0"
    fi
    g=$(column L2Norm\(Gx\) "$out/check/convergence.asc" 1)
    if le "$g" "$tol"; then report "$c" fresh PASS "residual $g at $param = $mu_last, rebuilt from scratch"
    else report "$c" fresh FAIL "residual $g > $tol at $param = $mu_last: parameter not followed in continuation?"; fi
    check_format "$c" "$out/$(awk '!/^ *%/ { d = $NF } END { print d }' "$out/MuD.asc")/ubest.nc" eq1_16x35x32.nc
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

# run_unit: the CTest suite, run in the build directory
run_unit() {
    local log=$PWD/run/unit.log
    mkdir -p run
    echo "== unit (ctest)"
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

# run_mpi: the MPI build's CTest suite, then eq1, eq3 and traj1 on 2 MPI processes against the
# same references
run_mpi() {
    if [ -z "$mpibuild" ]; then
        report mpi build INFO "no MPI build given (-m DIR), skipped"
        return
    fi
    local fs=$FINDSOLN sf=$SIMULATEFLOW layout np0 np1 tag
    # the MPI build's CTest suite (the serial tests on 1 process, the mpi_* tests on 2): ~2 min
    local log=$PWD/run/mpi_ctest.log
    mkdir -p run
    if (cd "$mpibuild" && ctest -j1 > "$log" 2>&1); then
        report mpi ctest PASS "$(grep -E 'tests passed' "$log")"
    else
        report mpi ctest FAIL "$(grep -E 'tests passed' "$log"); see $log"
    fi
    # RR_MPI_FINDSOLN / RR_MPI_SIMULATEFLOW override the MPI programs (e.g. a module's versions)
    local mfs=${RR_MPI_FINDSOLN:-$mpibuild/programs/findsoln}
    local msf=${RR_MPI_SIMULATEFLOW:-$mpibuild/programs/simulateflow}
    # two process grids: 2x1 (the default) and 1x2 (a different transpose decomposition)
    for layout in 2:1 1:2; do
        IFS=: read -r np0 np1 <<< "$layout"
        tag=mpi${np0}x${np1}
        export FINDSOLN="$mpirun -np 2 $mfs -np0 $np0 -np1 $np1"
        export SIMULATEFLOW="$mpirun -np 2 $msf -np0 $np0 -np1 $np1"
        run_findsoln eq1 ${tag}_eq1
        run_findsoln eq3 ${tag}_eq3
        run_traj traj1 ${tag}_traj1
    done
    export FINDSOLN=$fs SIMULATEFLOW=$sf
}

for c in $cases; do
    if [ "$c" = unit ]; then run_unit; continue; fi
    if [ "$c" = mpi ]; then run_mpi; continue; fi
    if [ -f simulate_$c.sh ]; then run_traj "$c"; continue; fi
    if [ -f continue_$c.sh ]; then run_cont "$c"; continue; fi
    if [ -f schemes_$c.sh ]; then run_schemes "$c"; continue; fi
    if [ -f error_$c.sh ]; then run_error "$c"; continue; fi
    if [ -f floquet_$c.sh ]; then run_floquet "$c"; continue; fi
    if [ -f eigen_$c.sh ]; then run_eigen "$c"; continue; fi
    if [ -f findsoln_$c.sh ]; then run_findsoln "$c"; continue; fi
    echo "unknown case: $c"; nfail=$((nfail + 1))
done

if [ $nfail -eq 0 ]; then echo "ALL PASSED"; else echo "$nfail FAILED"; fi
exit $((nfail > 0))
