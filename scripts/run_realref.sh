#!/usr/bin/env bash
# Whole chain on a scripts/make_realref_data.py dataset (check on a real reference genome, see README.md):
#   1. merged run: one strain on PacBio + ONT RNA002 + ONT RNA004 (+ an NGS-only library) -> `iflas all`, all reports;
#   2. per library, for diff_as: preproc + mapping of the four libraries, the long-read chain on the PacBio library,
#      rMATS condA vs condB with GO, the diff report.
# Usage (in the iflas env): scripts/run_realref.sh <dataset dir> [threads]; logs and /usr/bin/time reports in <dataset>/logs.
set -euo pipefail
cd "$1"
threads=${2:-16}
mkdir -p logs
step() {
    local name=$1
    shift
    echo "$(date '+%F %T') $name"
    /usr/bin/time -v -o "logs/$name.time" "$@" > "logs/$name.log" 2>&1
}
step all iflas all -c config.yaml -s samples.tsv -j 1
step report iflas report -all -c config.yaml -s samples.tsv
awk -F'\t' 'NR == 1 || $7 != "None"' samples_conditions.tsv > samples_conditions_tgs.tsv
step cond_mapping iflas all -c config_conditions.yaml -s samples_conditions.tsv --stages preproc,mapping -j 1
step cond_lr iflas all -c config_conditions.yaml -s samples_conditions_tgs.tsv --stages collapse,refine,pu_filter,find_as -j 1
step diff_as iflas diff_as -c config_conditions.yaml -s samples_conditions.tsv -d condA,condB --pu-filter --go -t "$threads"
step report_diff iflas report -diff -c config_conditions.yaml -s samples_conditions.tsv
echo "$(date '+%F %T') done"
