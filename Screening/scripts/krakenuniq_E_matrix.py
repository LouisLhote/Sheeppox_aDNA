#!/usr/bin/env python3

import os
import glob
import csv
from collections import defaultdict

# =========================
# USER SETTINGS
# =========================

ROOT_DIR = "results"
OUT_CSV = "krakenuniq_E_matrix.csv"

# List ONLY the species you want as columns
TARGET_SPECIES = [
    "Sheeppox virus",
    "Lumpy skin disease virus",
    "Goatpox virus",
]

FILL_VALUE = 0.0   # use None or "NA" if preferred


# =========================
# SCRIPT
# =========================

def safe_float(x):
    try:
        return float(x)
    except:
        return None


# sample -> species -> E
data = defaultdict(dict)

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
            species = " ".join(parts[8:])

            if rank != "species":
                continue

            if species not in TARGET_SPECIES:
                continue

            if reads is None or kmers is None or cov is None or reads == 0:
                continue

            E = (kmers / reads) * cov
            data[sample][species] = E


# =========================
# WRITE MATRIX CSV
# =========================

samples = sorted(data.keys())

with open(OUT_CSV, "w", newline="") as out_f:
    writer = csv.writer(out_f)

    # Header
    writer.writerow(["sample"] + TARGET_SPECIES)

    for sample in samples:
        row = [sample]
        for sp in TARGET_SPECIES:
            row.append(data[sample].get(sp, FILL_VALUE))
        writer.writerow(row)

print(f"✅ Matrix written to {OUT_CSV}")
