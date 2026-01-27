#!/usr/bin/env bash
set -euo pipefail

########################
# CONFIG
########################
INDEX="/raid_md0/Reference_Genomes/Ruminant_chimera/rumiv1"      # Bowtie2 index prefix
THREADS=2              # threads per job
MAX_JOBS=6              # parallel libraries
OUTDIR="host_depleted"

mkdir -p "$OUTDIR"

########################
# FUNCTION
########################
process_library() {
    BASENAME="$1"

    echo "[INFO] Processing $BASENAME"

    bowtie2 \
        -x "$INDEX" \
        -1 "${BASENAME}.pair1.truncated.gz" \
        -2 "${BASENAME}.pair2.truncated.gz" \
        --threads "$THREADS" \
        --very-sensitive \
        2> "$OUTDIR/${BASENAME}.bowtie2.log" | \
    samtools view -@ "$THREADS" -b | \
    samtools sort -@ "$THREADS" -o "$OUTDIR/${BASENAME}.sorted.bam"

    samtools view -@ "$THREADS" -b \
        -f 12 -F 256 \
        "$OUTDIR/${BASENAME}.sorted.bam" > \
        "$OUTDIR/${BASENAME}.unmapped_pairs.bam"

    samtools fastq -@ "$THREADS" \
        -1 "$OUTDIR/${BASENAME}_R1.unmapped.fastq.gz" \
        -2 "$OUTDIR/${BASENAME}_R2.unmapped.fastq.gz" \
        -0 /dev/null \
        -s /dev/null \
        -n "$OUTDIR/${BASENAME}.unmapped_pairs.bam"

    echo "[INFO] Finished $BASENAME"
}

export -f process_library
export INDEX THREADS OUTDIR

########################
# PARALLEL EXECUTION
########################
ls {YG01,YG15,80R,79RR}*.pair1.truncated.gz | \
sed 's/\.pair1\.truncated\.gz//' | \
sort -u | \
xargs -n 1 -P "$MAX_JOBS" -I {} bash -c 'process_library "$@"' _ {}
