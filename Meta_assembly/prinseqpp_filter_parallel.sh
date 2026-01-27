#!/usr/bin/env bash
set -euo pipefail

THREADS=2
MAX_JOBS=2        # number of libraries to run in parallel
OUTDIR="prinseqpp_filtered"
mkdir -p "$OUTDIR"

process_library() {
    BASENAME="$1"
    R1="${BASENAME}_R1.unmapped.fastq.gz"
    R2="${BASENAME}_R2.unmapped.fastq.gz"
    OUT1="$OUTDIR/${BASENAME}_R1.filtered"
    OUT2="$OUTDIR/${BASENAME}_R2.filtered"

    echo "[INFO] Running PRINSEQ++ on $BASENAME"

    prinseq++ \
        -fastq "$R1" \
        -fastq2 "$R2" \
        -out_good "$OUT1" \
        -out_good2 "$OUT2" \
        -threads $THREADS

    echo "[INFO] Finished $BASENAME"
}

export -f process_library
export THREADS OUTDIR

# Get unique basenames from R1 files
ls {YG01,YG15,79RR,80R}*_R1.unmapped.fastq.gz | sed 's/_R1.unmapped.fastq.gz//' | sort -u | \
    xargs -n 1 -P $MAX_JOBS -I {} bash -c 'process_library "$@"' _ {}
