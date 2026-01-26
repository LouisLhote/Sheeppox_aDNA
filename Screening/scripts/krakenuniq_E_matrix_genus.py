#!/usr/bin/env python3

import os
import glob
import csv
from collections import defaultdict

# =========================
# USER SETTINGS
# =========================

ROOT_DIR = "results"
OUT_CSV = "krakenuniq_E_matrix_genus.csv"

# Target genus to screen for
TARGET_GENUS = "Capripoxvirus"

FILL_VALUE = 0.0   # use None or "NA" if preferred


# =========================
# SCRIPT
# =========================

def safe_float(x):
    try:
        return float(x)
    except:
        return None


# sample -> genus -> E
# sample -> genus -> reads count
data = defaultdict(dict)
reads_data = defaultdict(dict)

reports = glob.glob(
    os.path.join(ROOT_DIR, "**/krakenuniq/kraken-report.txt"),
    recursive=True
)

for report in reports:
    sample = report.split(os.sep)[1]

    with open(report) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split()
            if len(parts) < 9:
                continue

            reads = safe_float(parts[1])
            kmers = safe_float(parts[3])
            cov = safe_float(parts[5])
            rank = parts[7]
            taxon_name = " ".join(parts[8:])

            if rank != "genus":
                continue

            if taxon_name != TARGET_GENUS:
                continue

            if reads is None or kmers is None or cov is None or reads == 0:
                continue

            E = (kmers / reads) * cov
            data[sample][TARGET_GENUS] = E
            reads_data[sample][TARGET_GENUS] = int(reads)


# =========================
# WRITE MATRIX CSV
# =========================

samples = sorted(data.keys())

with open(OUT_CSV, "w", newline="") as out_f:
    writer = csv.writer(out_f)

    # Header: sample, E value, reads count
    writer.writerow(["sample", f"{TARGET_GENUS}_E", f"{TARGET_GENUS}_reads"])

    for sample in samples:
        row = [sample]
        # E value
        row.append(data[sample].get(TARGET_GENUS, FILL_VALUE))
        # Reads count
        row.append(reads_data[sample].get(TARGET_GENUS, 0))
        writer.writerow(row)

print(f"✅ Matrix written to {OUT_CSV}")
