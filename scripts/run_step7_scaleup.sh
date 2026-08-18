#!/bin/sh
# Step 7 scale-up: run the comparator across every module that passed the screen.
#
# Reduced relative to Step 1 on purpose. Step 1 needed ten oracle seeds because it was
# estimating a ceiling precisely enough to read a branch off it. Step 7 needs a ceiling
# per module accurate enough to regress against the diagnostics, across many modules, so
# the budget moves from depth to breadth: three seeds and a 300-structure null. At Step 1
# settings, 28 modules would be several days of compute.
#
# Sources are limited to what the regime map needs: the ceiling, the two baselines that
# define the bar, the floor, and the null. Textbook and proposal arms are not defined for
# generated modules anyway.
#
# Idempotent: a module with a complete result is skipped, so this is safe to restart
# after a VPN drop or a host reboot. Completeness is checked on content, not on the
# presence of a file, because a killed run leaves a partial one behind.
set -u
WORK="$HOME/mmc_work"
RES="$WORK/results/step7"
LOGS="$WORK/logs/step7"
STATUS="$WORK/step7_status.txt"
MODULES="$WORK/results/step7_modules.json"
CPUS=12
SEEDS=3
NRAND=300

mkdir -p "$RES" "$LOGS"
note() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ)  $*" >> "$STATUS"; }

# A result file appears as soon as the run flushes its first section, so file existence
# does not mean the module finished. A host reboot part-way through leaves a parseable
# but empty file, and testing existence alone would skip that module for good on the
# restart. Completeness is therefore read off the content: a finished run carries a
# populated comparison table and the random null it is judged against.
complete() {
    python3 - "$1" <<'CHECK' 2>/dev/null
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(1)
sys.exit(0 if d.get("table") and d.get("random_null") else 1)
CHECK
}


# names of the candidates that passed the screen, and a gene-list file for each
python3 - "$MODULES" "$WORK/step7_defs" <<'PY'
import json, os, sys
src, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
d = json.load(open(src))
passed = [c for c in d["candidates"] if c["passed"]]
for c in passed:
    with open(os.path.join(outdir, c["name"] + ".json"), "w") as f:
        json.dump({"module": c["name"], "regulators": c["genes"],
                   "targets": c["genes"]}, f)
print("\n".join(c["name"] for c in passed))
PY

note "step7 scale-up started"
for name in $(python3 -c "
import json,sys
d=json.load(open('$MODULES'))
print('\n'.join(c['name'] for c in d['candidates'] if c['passed']))
"); do
    out="$RES/step1_${name}.json"
    if complete "$out"; then
        note "skip $name (complete)"
        continue
    fi
    note "start $name"
    docker rm -f "mmc_s7_$name" >/dev/null 2>&1
    docker run --rm --name "mmc_s7_$name" --cpus="$CPUS" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        mmc:v4 python -u scripts/step1_comparator.py \
            --module "$name" --condition Stim8hr \
            --module-def "/work/step7_defs/${name}.json" \
            --sources linear,mean,zero,oracle,random \
            --workers "$CPUS" --seeds "$SEEDS" --n-random "$NRAND" \
            --out "/work/results/step7/step1_${name}.json" \
        > "$LOGS/${name}.log" 2>&1
    if complete "$out"; then note "done $name"; else note "FAILED $name (no usable result)"; fi
done
note "step7 scale-up finished"
