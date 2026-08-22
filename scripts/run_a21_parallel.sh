#!/bin/sh
# A20, parallel form. One worker per shard of the machine, self-balancing by claim.
#
# The sequential runner used 12 of 32 cores and left the rest idle. The oracle search
# parallelises over candidate structures within a module, and that scales sublinearly, so
# several modules at moderate width finishes sooner than one module at full width.
#
# Workers are identical and each walks the whole module list. A module is taken by creating
# its claim directory, which is atomic under POSIX, so two workers never take the same one.
# A claim is released when the run finishes without producing a complete result, so a
# killed worker does not strand a module permanently.
#
# BLAS threading is pinned to one thread per process. The search parallelises by spawning a
# pool of processes, and each one otherwise starts its own BLAS thread pool sized to the
# whole machine, so the thread count multiplies out to the hundreds and the cores spend
# their time context switching. Measured on this machine, four containers at eight workers
# drove the load average past 125 before this was set.
#
#   sh run_a20_parallel.sh <worker-id>
#
set -u
WORKER="${1:-0}"
WORK="$HOME/mmc_work"
STATUS="$WORK/a21_status.txt"
LOGS="$WORK/logs/a21"
CLAIMS="$WORK/a21_claims"
CPUS=8
SEEDS=20
NRAND=300

# heaviest first, so the long 40-perturbation modules start before the short ones and the
# tail of the run is short rather than a single large module finishing alone
MODULES="regulon_AHR regulon_HIF1A regulon_NFE2L2 regulon_STAT3 regulon_YY1 \
Cytokine_production coresponse_ACTR2 coresponse_CFAP298 coresponse_HCCS coresponse_KIF20A \
coresponse_MOV10 coresponse_PIM1 TCR_signalosome"

mkdir -p "$LOGS" "$CLAIMS" "$WORK/results/a21"

note() { printf '%s [w%s] %s\n' "$(date -u +%H:%M:%S)" "$WORKER" "$1" >> "$STATUS"; }

defarg() {
    case "$1" in
        Cytokine_production) printf -- '--module-def /work/cytokine_module_def.json' ;;
        TCR_signalosome|CD4_lineage_TFs|Th2_GATA3) printf '' ;;
        *) printf -- '--module-def /work/step7_defs/%s.json' "$1" ;;
    esac
}

complete() {
    python3 - "$1" <<'CHECK' 2>/dev/null
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(1)
sys.exit(0 if d.get("table") and d.get("oracle_seeds") else 1)
CHECK
}

note "worker started"
for name in $MODULES; do
    # host-side path, used only for the completeness check. The container is given the
    # path under its own /work mount: passing this one through would make the container
    # create the tree inside itself and lose the result on --rm, which is what happened
    # on 2026-08-22 and cost four module-runs that had actually succeeded.
    out="$WORK/results/a21/step1_${name}.json"
    complete "$out" && continue
    mkdir "$CLAIMS/$name" 2>/dev/null || continue      # already taken by another worker
    note "start $name"
    docker rm -f "mmc_a21_$name" >/dev/null 2>&1
    docker run --rm --name "mmc_a21_$name" --cpus="$CPUS" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
        -e NUMEXPR_NUM_THREADS=1 \
        mmc:v4 python -u scripts/step1_comparator.py \
            --module "$name" --condition Stim8hr $(defarg "$name") \
            --sources linear,mean,zero,oracle,random \
            --workers "$CPUS" --seeds "$SEEDS" --n-random "$NRAND" \
            --out "/work/results/a21/step1_${name}.json" \
        > "$LOGS/${name}.log" 2>&1
    if complete "$out"; then
        note "done $name"
    else
        note "FAILED $name, releasing claim"
        rmdir "$CLAIMS/$name" 2>/dev/null
    fi
done
note "worker finished"
