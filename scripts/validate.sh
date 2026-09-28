#!/usr/bin/env bash
# Validate nanorust against the R pipeline outputs on a (per-chromosome) BAM subset.
# Usage: validate.sh sample.bam ref_nanoT_alldata.rds ref_coverage.bw chrom1[,chrom2...]
set -euo pipefail
bam=$1 ref_rds=$2 ref_bw=$3 chroms=$4
here=$(cd "$(dirname "$0")" && pwd)
bin="$here/../target/release/nanorust"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$bin" run -b "$bam" -o "$tmp/test"
Rscript "$here/compare_alldata.R" "$tmp/test_nanoT_alldata.rds" "$ref_rds" "$chroms"
echo "--- coverage bigWig (only listed chromosomes are meaningful)"
"$bin" bw-compare "$tmp/test.bw" "$ref_bw" | grep -E "^($(echo "$chroms" | tr ',' '|'))\s"
