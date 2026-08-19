#!/bin/sh
# The A5 extension of PREREG_v4 amendment A16, stratified and restartable.
#
# Two strata that are never pooled. The curated modules carry hand-written biological
# context, so A1 there proposes from gene symbols plus a paragraph of regulatory biology.
# The generated modules fall through to a single generic sentence, so A1 there proposes
# from symbols and almost nothing. Pooling them would measure A5 partly against a
# handicapped A1 and inflate the gain, which is why they run and report separately.
#
# Only A1 and A5 are drawn: the anonymisation arms are established and add nothing to the
# question of whether showing the proposer data changes what it proposes.
#
# Each module writes its own file and the underlying script skips any arm-seed already
# recorded, so this resumes after a reboot rather than repeating.
set -u
WORK="$HOME/mmc_work"
RES="$WORK/results/step5_ext"
LOGS="$WORK/logs/step5_ext"
STATUS="$WORK/step5_ext_status.txt"
SEEDS=5
CPUS=8

mkdir -p "$RES" "$LOGS"
note() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ)  $*" >> "$STATUS"; }

# a module is finished when its file carries an agreement matrix
complete() {
    [ -f "$1" ] || return 1
    python3 - "$1" <<'CHECK' 2>/dev/null
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(1)
sys.exit(0 if d.get("agreement") else 1)
CHECK
}

run_module() {
    name="$1"; defarg="$2"
    out="$RES/step5_$name.json"
    if complete "$out"; then
        note "skip $name (complete)"
        return 0
    fi
    if docker ps -a --format '{{.Names}}' | grep -qx "mmc_s5e_$name"; then
        note "skip $name (container exists)"
        return 0
    fi
    note "start $name"
    docker run --rm --name "mmc_s5e_$name" --cpus="$CPUS" \
        --env-file "$HOME/.verdict.env" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        mmc:v4 python -u scripts/step5_anonymization.py \
            --modules "$name" --seeds "$SEEDS" --arms A1,A5 $defarg \
            --out "/work/results/step5_ext/step5_$name.json" \
        >> "$LOGS/$name.log" 2>&1
    if complete "$out"; then note "done $name"; else note "INCOMPLETE $name"; fi
    return 0
}

note "A5 extension started, stratum ${STRATUM:-generated}, $SEEDS seeds, arms A1 and A5"

# Curated stratum, the primary one: real biological context in the prompt.
#
# Skipped while the four-module A5 job is still producing seeds 0 to 2 on these same
# modules, because running both would draw the same arm-seeds twice against the model API
# for no gain. Pass "curated" once that job is finished to extend them to five seeds.
if [ "${STRATUM:-generated}" = "curated" ]; then
    run_module Cytokine_production "--module-def /work/cytokine_module_def.json"
    for m in TCR_signalosome CD4_lineage_TFs Th2_GATA3; do
        run_module "$m" ""
    done
    note "curated stratum finished"
    exit 0
fi

# generated stratum, reported separately: generic context, chosen because their linear arm
# beats its own null so the comparison against linear stays meaningful on them
for m in coresponse_PIM1 coresponse_MOV10 coresponse_HCCS \
         regulon_AHR regulon_STAT3 regulon_YY1; do
    run_module "$m" "--module-def /work/step7_defs/$m.json"
done

note "A5 extension finished"
