import os
import sys
import csv
import pysam
from glob import glob
import argparse

def count_bam_reads(path):
    if not os.path.exists(path):
        return 0
    with pysam.AlignmentFile(path, "rb") as bam:
        return sum(1 for _ in bam.fetch(until_eof=True))

def count_fastq_reads(path):
    if not os.path.exists(path):
        return 0
    count = 0
    with pysam.FastxFile(path) as fq:
        for _ in fq:
            count += 1
    return count

def main(samples, output_file):
    results = []
    for sample in samples:
        merged_bam = f"results/{sample}/merged/{sample}_merged.bam"
        if not os.path.exists(merged_bam):
            # Skip samples without merged bam file
            continue

        # Gather all PCRs belonging to this sample
        pcr_bams = glob(f"results/*/final/*_q30_rmdup.bam")
        pcrs = [os.path.basename(p).split("_")[0] for p in pcr_bams if sample in p]

        if not pcrs:
            print(f"[WARNING] No PCR runs found for sample {sample}", file=sys.stderr)
            continue

        raw_reads_total = 0
        host_filtered_total = 0
        host_unaligned = 0
        sheeppox_f4 = 0
        sheeppox_q30 = 0
        sheeppox_rmdup = 0

        for pcr in pcrs:
            # 0. Raw R1 reads
            r1 = f"data/{pcr}_R1.fastq.gz"
            raw_reads_total += count_fastq_reads(r1)

            # 1. Host-filtered input (collapsed reads)
            collapsed = f"results/{pcr}/adapter_removal/{pcr}.collapsed.gz"
            host_filtered_total += count_fastq_reads(collapsed)

            # 2. Host-unaligned reads (f4.bam)
            f4_bam = f"results/{pcr}/host_cleaned/{pcr}_f4.bam"
            host_unaligned += count_bam_reads(f4_bam)

            # 3. Sheeppox reads (aligned F4)
            f4 = f"results/{pcr}/bam/{pcr}_Sheeppoxvirus_F4.bam"
            sheeppox_f4 += count_bam_reads(f4)

            # 4. Q30
            q30 = f"results/{pcr}/q30/{pcr}_Sheeppoxvirus_F4_q30.bam"
            sheeppox_q30 += count_bam_reads(q30)

            # 5. Deduplicated
            rmdup = f"results/{pcr}/final/{pcr}_Sheeppoxvirus_F4_q30_rmdup.bam"
            sheeppox_rmdup += count_bam_reads(rmdup)

        host_mapped = host_filtered_total - host_unaligned
        dup_rate = (sheeppox_q30 - sheeppox_rmdup) / sheeppox_q30 if sheeppox_q30 > 0 else 0

        results.append({
            "Sample": sample,
            "Raw R1 Reads": raw_reads_total,
            "Collapsed Reads (host input)": host_filtered_total,
            "Host-Unaligned Reads": host_unaligned,
            "Host-Mapped Reads": host_mapped,
            "Mapped to Sheeppox (F4)": sheeppox_f4,
            "Q30 Reads": sheeppox_q30,
            "Deduplicated Reads": sheeppox_rmdup,
            "Duplicate Rate": round(dup_rate, 4)
        })

    if not results:
        print("[ERROR] No metrics collected for the given samples.", file=sys.stderr)
        sys.exit(1)

    # Write to CSV
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect sample metrics")
    parser.add_argument("samples", nargs="+", help="Sample names")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    main(args.samples, args.output)
