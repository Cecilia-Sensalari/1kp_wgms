#!/usr/bin/env bash
# Compact progress of a Nextflow run, one line per process, from its trace file (trace/trace_<date>_<time>.txt in
# the folder the run was launched from, written by the "trace" block of batch_ksrates/nextflow.config).
# The trace file has a line for every task that has ENDED, not for the running ones.
#
# Usage:   nextflow_progress.sh [TRACE_FILE]         (default: the newest trace/trace_*.txt of the current folder)
# Live:    watch -n 30 nextflow_progress.sh
# Output:  process, tasks completed, tasks failed (with errorStrategy 'ignore' they do not stop the run),
#          cached, other (e.g. aborted)

set -euo pipefail

trace=${1:-$(ls -t trace/trace_*.txt 2>/dev/null | head -n 1 || true)}
[ -n "$trace" ] && [ -r "$trace" ] || { echo "ERROR: no trace file found (give it as argument)" >&2; exit 1; }

echo "$trace"
awk -F'\t' '
    NR == 1 { for (i = 1; i <= NF; i++) column[$i] = i; next }
    {
        process = $column["name"]; sub(/ \(.*\)$/, "", process)
        status = $column["status"]
        if (!(process in seen)) { seen[process] = 1; order[++n] = process }
        if (status == "COMPLETED") completed[process]++
        else if (status == "FAILED") failed[process]++
        else if (status == "CACHED") cached[process]++
        else other[process]++
    }
    END {
        printf "%-24s %10s %8s %8s %8s\n", "process", "completed", "failed", "cached", "other"
        for (i = 1; i <= n; i++) { p = order[i]; printf "%-24s %10d %8d %8d %8d\n", p, completed[p], failed[p], cached[p], other[p] }
    }
' "$trace"
