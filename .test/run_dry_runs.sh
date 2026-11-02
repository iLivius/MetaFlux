#!/usr/bin/env bash
# Dry-run every config of the CI matrix (.test/config/) and check the outcome.
#
# What it is for: a dry run (snakemake -n) reads the workflow, runs every start-up
# check and lists the jobs it would run, without running any. So it proves, in a few
# seconds per config and without conda environments or reference databases, that
# every marker, taxonomy method, SILVA release and shotgun profiler still builds its
# plan after a change, and that the start-up checks still stop what they must stop.
# Two groups of configs:
#   .test/config/*.yaml                must pass: exit code 0 and a job table;
#   .test/config/expected_stop/*.yaml  must stop at start-up with the message written
#                                      in the config's "# EXPECT:" line; passing, or
#                                      stopping for another reason, is a failure.
# Then one setting typed inline with --config, switched off and on, must change the
# plan as it would written in a config file (see that section below).
# How: each config is dry-run from the repository root (the Snakefile needs that),
# with the template config/config.yaml underneath it, as in a real run. The snakemake
# command can be given with the SNAKEMAKE environment variable (default: the one on
# PATH). No conda is needed: the --sdm conda flag is left out, which only matters
# when jobs run, and a dry run runs none.
# Output: one line per config on the terminal (PASS/FAIL, the number of jobs or the
# stop message) and one log per config under .test/out/dry_run_logs/ (git-ignored).
# Exit code: the number of failures, so 0 means everything passed.
# Used by .github/workflows/ci.yml; run it locally the same way:
#     bash .test/run_dry_runs.sh
set -u

# Work from the repository root, wherever the script was called from.
cd "$(dirname "$0")/.." || exit 1
SNAKEMAKE="${SNAKEMAKE:-snakemake}"
LOG_DIR=.test/out/dry_run_logs
mkdir -p "$LOG_DIR"

failures=0
echo "snakemake: $("$SNAKEMAKE" --version 2>/dev/null || echo "not found: $SNAKEMAKE")"
echo

# ── Configs that must pass ──
for config in .test/config/*.yaml; do
    name=$(basename "$config" .yaml)
    log="$LOG_DIR/$name.log"
    start=$(date +%s)
    if "$SNAKEMAKE" -n --cores 1 --configfile "$config" > "$log" 2>&1; then
        # The job table ends with a "total N" line (Snakemake prints it twice; the
        # first is enough).
        jobs=$(grep -m1 -E '^total[[:space:]]+[0-9]+' "$log" | awk '{print $2}')
        if [ -z "$jobs" ]; then
            echo "FAIL  $name: the dry run passed but printed no job table (see $log)"
            failures=$((failures + 1))
        else
            echo "PASS  $name: $jobs jobs planned ($(( $(date +%s) - start )) s)"
        fi
    else
        echo "FAIL  $name: the dry run stopped (see $log); last lines:"
        tail -n 5 "$log" | sed 's/^/        /'
        failures=$((failures + 1))
    fi
done
echo

# ── Configs that must stop at start-up, with the expected message ──
for config in .test/config/expected_stop/*.yaml; do
    name=$(basename "$config" .yaml)
    log="$LOG_DIR/expected_stop_$name.log"
    expected=$(grep -m1 '^# EXPECT:' "$config" | sed 's/^# EXPECT:[[:space:]]*//')
    if [ -z "$expected" ]; then
        echo "FAIL  $name: no '# EXPECT:' line in the config"
        failures=$((failures + 1))
        continue
    fi
    if "$SNAKEMAKE" -n --cores 1 --configfile "$config" > "$log" 2>&1; then
        echo "FAIL  $name: the dry run passed, but it should have stopped with: $expected"
        failures=$((failures + 1))
    elif grep -qF -- "$expected" "$log"; then
        echo "PASS  $name: stopped at start-up as expected ($expected)"
    else
        echo "FAIL  $name: stopped, but not with the expected message ($expected); last lines:"
        tail -n 5 "$log" | sed 's/^/        /'
        failures=$((failures + 1))
    fi
done
echo

# ── A setting typed inline with --config ──
# Snakemake hands MetaFlux every value of a block typed after --config as text, and
# 00_common.smk reads it again as YAML reads a config file ("Settings typed inline
# with --config" there). Checked here end to end with the phylogeny switch, which is
# read with bool(): typed off inline over a config that turns it on, the plan must
# hold no phylogeny job (the text "false" used to count as true); typed on inline over
# a config that leaves it off, the plan must hold phylogeny jobs, which also shows
# that the job names searched for exist. The phylogeny rules are named phylo_*, and
# each one has a line with its job count in the job table.
check_inline_switch() {
    local name=$1 config=$2 inline=$3 wanted=$4    # wanted: "none" or "some"
    local log="$LOG_DIR/inline_$name.log"
    if ! "$SNAKEMAKE" -n --cores 1 --configfile "$config" --config "$inline" > "$log" 2>&1; then
        echo "FAIL  inline $name: the dry run stopped (see $log); last lines:"
        tail -n 5 "$log" | sed 's/^/        /'
        failures=$((failures + 1))
        return
    fi
    local phylo_lines
    phylo_lines=$(grep -cE '^phylo_[a-z_]+[[:space:]]+[0-9]+' "$log")
    local passed=no
    local described="no phylogeny job"
    if [ "$wanted" = none ] && [ "$phylo_lines" -eq 0 ]; then passed=yes; fi
    if [ "$wanted" = some ]; then described="phylogeny jobs"; fi
    if [ "$wanted" = some ] && [ "$phylo_lines" -gt 0 ]; then passed=yes; fi
    if [ "$passed" = yes ]; then
        echo "PASS  inline $name: --config '$inline' plans $described"
    else
        echo "FAIL  inline $name: --config '$inline' should plan $described" \
             "($phylo_lines phylo_* lines in the job table, see $log)"
        failures=$((failures + 1))
    fi
}
check_inline_switch phylogeny_off .test/config/16S_sintax_silva144.yaml \
    'amplicon={phylogeny: {enabled: false}}' none
check_inline_switch phylogeny_on .test/config/16S_rdp_silva144.yaml \
    'amplicon={phylogeny: {enabled: true}}' some
echo

if [ "$failures" -eq 0 ]; then
    echo "all dry runs behaved as expected"
else
    echo "$failures dry run(s) did not behave as expected"
fi
exit "$failures"
