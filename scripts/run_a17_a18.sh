#!/bin/sh
# Amendments A17 and A18, from the third review.
#
# A18 runs first because it is cheap and because it carries a prediction fixed in the
# pre-registration before the run: across the 13 modules where the linear arm beats its
# own random null, the full ridge map clears the reduced-rank map on at most 4, and the
# median ratio of reduced-rank to full-rank held-out DE-overlap is at least 0.90. That
# converts the low-rank account from a diagnostic into a prediction, and a refutation is
# a problem for claim 2 rather than something to be reinterpreted.
#
# A17 replicates the generator sweep of section 2.34 on six further modules, because the
# load-bearing claim currently rests on Cytokine_production alone. The six were named in
# the pre-registration before this ran and are not revisited afterwards.
#
# Idempotent on content, not on the presence of a file, since a killed run leaves a
# partial one behind. Safe to restart after a VPN drop or a host reboot.
set -u
WORK="$HOME/mmc_work"
STATUS="$WORK/a17_a18_status.txt"
LOGS="$WORK/logs/a17_a18"
CPUS=12

A18_MODULES="coresponse_ACTR2 coresponse_CFAP298 coresponse_HCCS coresponse_KIF20A \
coresponse_MOV10 coresponse_PIM1 regulon_AHR regulon_HIF1A regulon_NFE2L2 regulon_STAT3 \
regulon_YY1 Cytokine_production TCR_signalosome"

A17_MODULES="coresponse_PIM1 coresponse_HCCS coresponse_ACTR2 TCR_signalosome \
CD4_lineage_TFs regulon_STAT3"

mkdir -p "$LOGS" "$WORK/results/step12" "$WORK/results/a17"

note() { printf '%s %s\n' "$(date -u +%H:%M:%S)" "$1" >> "$STATUS"; }

# built-in modules carry no definition file; generated ones live in step7_defs
defarg() {
    case "$1" in
        Cytokine_production) printf -- '--module-def /work/cytokine_module_def.json' ;;
        TCR_signalosome|CD4_lineage_TFs|Th2_GATA3) printf '' ;;
        *) printf -- '--module-def /work/step7_defs/%s.json' "$1" ;;
    esac
}

complete_key() {
    python3 - "$1" "$2" <<'CHECK' 2>/dev/null
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(1)
sys.exit(0 if d.get(sys.argv[2]) else 1)
CHECK
}

note "A17/A18 runner started"

# ---------------------------------------------------------------- A18
for name in $A18_MODULES; do
    out="$WORK/results/step12/${name}.json"
    if complete_key "$out" tests; then note "skip A18 $name (complete)"; continue; fi
    note "start A18 $name"
    docker rm -f "mmc_a18_$name" >/dev/null 2>&1
    docker run --rm --name "mmc_a18_$name" --cpus="$CPUS" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        mmc:v4 python -u scripts/step12_reduced_rank.py \
            --module "$name" --condition Stim8hr $(defarg "$name") \
            --out "/work/results/step12/${name}.json" \
        > "$LOGS/a18_${name}.log" 2>&1
    if complete_key "$out" tests; then note "done A18 $name"; else note "FAILED A18 $name"; fi
done
note "A18 finished"

# ---------------------------------------------------------------- A17
for name in $A17_MODULES; do
    out="$WORK/results/a17/${name}.json"
    if complete_key "$out" variants; then note "skip A17 $name (complete)"; continue; fi
    note "start A17 $name"
    docker rm -f "mmc_a17_$name" >/dev/null 2>&1
    docker run --rm --name "mmc_a17_$name" --cpus="$CPUS" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        mmc:v4 python -u scripts/step3_calibrate.py \
            --module "$name" --condition Stim8hr $(defarg "$name") \
            --edges 48 --seeds 3 \
            --out "/work/results/a17/${name}.json" \
        > "$LOGS/a17_${name}.log" 2>&1
    if complete_key "$out" variants; then note "done A17 $name"; else note "FAILED A17 $name"; fi
done
note "A17 finished"
note "A17/A18 runner finished"
