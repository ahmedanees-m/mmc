#!/bin/sh
# A second Step 7 worker, taking the queue from the far end.
#
# The host has 32 cores and the original queue uses 12 of them, so most of the machine
# sits idle while a module runs. This worker walks the module list in reverse so it meets
# the original queue in the middle rather than racing it from the start.
#
# The two never work on the same module because both refuse to start one whose output
# path already exists, and this worker creates that file as a claim before it launches
# anything. The original queue tests the same path, so a claimed module reads to it as
# already handled. The claim is replaced by the real result on success and removed on
# failure, so a module is never left permanently claimed and unrun.
set -u
WORK="$HOME/mmc_work"
RES="$WORK/results/step7"
LOGS="$WORK/logs/step7"
STATUS="$WORK/step7_status.txt"
CPUS=10
SEEDS=3
NRAND=300

note() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ)  [w2] $*" >> "$STATUS"; }

note "worker 2 started, taking the queue in reverse"

# A claim only means something while the worker that wrote it is alive. If this instance
# is starting, no claim of its own can be in flight, so any left over came from a run the
# host killed and would otherwise hide its module from both workers permanently. They are
# cleared here. Held-out placeholders carry a different marker and are left alone.
python3 - <<'CLEAN' >> "$STATUS" 2>/dev/null
import json, glob, os, datetime
stamp = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
for f in glob.glob(os.path.expanduser("~/mmc_work/results/step7/step1_*.json")):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    if d.get("claimed_by") and not d.get("table"):
        os.remove(f)
        print("%s  [w2] released stale claim on %s" % (stamp, d.get("module")))
for f in glob.glob(os.path.expanduser("~/mmc_work/results/step7/step1_*.tmp.json")):
    os.remove(f)
CLEAN

for name in $(python3 -c "
import json
d = json.load(open('$WORK/results/step7_modules.json'))
names = [c['name'] for c in d['candidates'] if c['passed']]
print('\n'.join(reversed(names)))
"); do
    out="$RES/step1_${name}.json"
    # anything with a file is finished, claimed, held out, or in flight elsewhere
    if [ -f "$out" ]; then
        continue
    fi
    if docker ps -a --format '{{.Names}}' | grep -qx "mmc_s7_$name"; then
        continue
    fi
    printf '{"claimed_by": "worker2", "module": "%s"}' "$name" > "$out"
    note "start $name"
    docker run --rm --name "mmc_s7_$name" --cpus="$CPUS" \
        -v "$HOME/mmc:/app" -v "$WORK:/work" \
        -v "$HOME/verdict/data/store:/store:ro" \
        -e MMC_ZHU_STORE=/store/zhu -e PYTHONPATH=/app -w /app \
        mmc:v4 python -u scripts/step1_comparator.py \
            --module "$name" --condition Stim8hr \
            --module-def "/work/step7_defs/${name}.json" \
            --sources linear,mean,zero,oracle,random \
            --workers "$CPUS" --seeds "$SEEDS" --n-random "$NRAND" \
            --out "/work/results/step7/step1_${name}.tmp.json" \
        > "$LOGS/${name}.log" 2>&1
    if python3 -c "
import json, sys
try:
    d = json.load(open('$RES/step1_${name}.tmp.json'))
except Exception:
    sys.exit(1)
sys.exit(0 if d.get('table') and d.get('random_null') else 1)
" 2>/dev/null; then
        mv "$RES/step1_${name}.tmp.json" "$out"
        note "done $name"
    else
        rm -f "$out" "$RES/step1_${name}.tmp.json"
        note "FAILED $name, claim released"
    fi
done
note "worker 2 finished"
