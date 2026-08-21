#!/bin/sh
# Amendments A19 and A20, from the third review.
#
# A19 is diagnostics only and runs first because it is cheap: effective rank, leading-PC
# fraction and the perturbation-specific ratio on the two atlas states this project has not
# used, across all 27 modules. It tests no hypothesis. It establishes whether the signature
# behind claim 2 is a property of this module set or of perturbation response matrices more
# generally, and a state that does not show it is reported unadjusted.
#
# A20 re-runs the comparator on the 13 modules where the linear arm beats its own random
# null, this time retaining per-seed per-fold scores, which section 2.31 recorded as never
# stored. Only those 13, because the other 14 cannot support the comparison and re-running
# them would not change any reported count.
#
# Seeds are 5 rather than the 3 the Step 7 scale-up used. A median over 3 seeds is a weak
# statistic and the median is now the primary summary, so the extra depth is the point of
# the re-run. This changes the ceiling estimate for reasons beyond the retention fix, and
# that is stated rather than glossed.
#
# Idempotent on content. Safe to restart after a VPN drop or a host reboot.
set -u
WORK="$HOME/mmc_work"
STATUS="$WORK/a19_a20_status.txt"
LOGS="$WORK/logs/a19_a20"
CPUS=12
SEEDS=5
NRAND=300

ALL_MODULES="coresponse_ACTR2 coresponse_CFAP298 coresponse_ELAVL1 coresponse_HCCS \
coresponse_KIF20A coresponse_MBD5 coresponse_MOV10 coresponse_MTA2 coresponse_NATD1 \
coresponse_PIM1 coresponse_RPRD2 coresponse_SHOC2 coresponse_UBE2I regulon_AHR \
regulon_CREB1 regulon_DNMT1 regulon_EGR1 regulon_ETS1 regulon_HIF1A regulon_NFE2L2 \
regulon_STAT3 regulon_USF1 regulon_YY1 Cytokine_production TCR_signalosome Th2_GATA3 \
CD4_lineage_TFs"

SOUND_MODULES="coresponse_ACTR2 coresponse_CFAP298 coresponse_HCCS coresponse_KIF20A \
coresponse_MOV10 coresponse_PIM1 regulon_AHR regulon_HIF1A regulon_NFE2L2 regulon_STAT3 \
regulon_YY1 Cytokine_production TCR_signalosome"

mkdir -p "$LOGS" "$WORK/results/a19" "$WORK/results/a20"

note() { printf '%s %s\n' "$(date -u +%H:%M:%S)" "$1" >> "$STATUS"; }

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

note "A19/A20 runner started"

# ---------------------------------------------------------------- A19
for name in $ALL_MODULES; do
    out="$WORK/results/a19/zhu_${name}.json"
    if complete_key "$out" states; then note "skip A19 $name (complete)"; continue; fi
    docker rm -f "mmc_a19_$name" >/dev/null 2>&1
    docker run --rm --name "mmc_a19_$name" --cpus="$CPUS" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        mmc:v4 python -u scripts/step13_external_diagnostics.py \
            --source zhu --module "$name" $(defarg "$name") \
            --conditions Rest,Stim8hr,Stim48hr \
            --out "/work/results/a19/zhu_${name}.json" \
        > "$LOGS/a19_${name}.log" 2>&1
    if complete_key "$out" states; then note "done A19 $name"; else note "FAILED A19 $name"; fi
done
note "A19 finished"

# ---------------------------------------------------------------- A20
for name in $SOUND_MODULES; do
    out="$WORK/results/a20/step1_${name}.json"
    if complete_key "$out" table; then note "skip A20 $name (complete)"; continue; fi
    note "start A20 $name"
    docker rm -f "mmc_a20_$name" >/dev/null 2>&1
    docker run --rm --name "mmc_a20_$name" --cpus="$CPUS" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        mmc:v4 python -u scripts/step1_comparator.py \
            --module "$name" --condition Stim8hr $(defarg "$name") \
            --sources linear,mean,zero,oracle,random \
            --workers "$CPUS" --seeds "$SEEDS" --n-random "$NRAND" \
            --out "/work/results/a20/step1_${name}.json" \
        > "$LOGS/a20_${name}.log" 2>&1
    if complete_key "$out" table; then note "done A20 $name"; else note "FAILED A20 $name"; fi
done
note "A20 finished"
note "A19/A20 runner finished"
