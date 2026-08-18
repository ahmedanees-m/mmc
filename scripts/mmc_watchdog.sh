#!/bin/sh
# Restart the MMC job queues if they are not running.
#
# The host reboots irregularly (twice on 2026-08-18, roughly four hours apart), which
# kills every container and both queues with them. The queues are idempotent and skip
# completed work, so relaunching them after a reboot resumes from where they stopped.
# Follows the same watchdog-plus-crontab pattern already used on this host.
#
# Each queue is checked independently and neither is touched while it is alive, so this
# is safe to run every few minutes. Worker 2 clears its own stale claims on startup,
# which is what makes a reboot mid-module recoverable rather than silently skipped.
WORK="$HOME/mmc_work"

alive() { [ "$(ps -eo args | grep "$1" | grep -v grep | wc -l)" -gt 0 ]; }

if [ -f "$WORK/orchestrator.sh" ] && ! alive orchestr; then
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ)  orchestrator not running, relaunching"
    (setsid nohup sh "$WORK/orchestrator.sh" > "$WORK/logs/orchestrator.log" 2>&1 < /dev/null &)
fi

# Worker 2 only has work while modules remain unrun; it exits on its own when the queue
# is exhausted, so it is relaunched only while the scale-up as a whole is unfinished.
if [ -f "$WORK/run_step7_worker2.sh" ] && ! alive worker2 \
   && ! grep -q "step7 scale-up finished" "$WORK/step7_status.txt" 2>/dev/null; then
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ)  worker 2 not running, relaunching"
    (setsid nohup sh "$WORK/run_step7_worker2.sh" > "$WORK/logs/worker2.log" 2>&1 < /dev/null &)
fi
