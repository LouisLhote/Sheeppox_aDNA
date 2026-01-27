import os
import glob
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Set

import numpy as np
import pandas as pd


def read_fasta(path: str) -> List[Tuple[str, str]]:
    """Minimal FASTA parser: returns list of (id, seq)."""
    records: List[Tuple[str, str]] = []
    header = None
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks)))

    return records


def load_metadata(metadata_file: str) -> Dict[str, str]:
    """
    Load metadata and create mapping from accession to virus/group.

    Uses the same structure as in analyze_shared_mutations.py:
    - column 'accession' for the accession
    - column 'Virus' for the group (e.g. GTPV, SPPV, LSDV)
    """
    df = pd.read_csv(metadata_file, sep="\t")
    meta_map: Dict[str, str] = {}
    for _, row in df.iterrows():
        acc = str(row["accession"]).strip()
        virus = str(row["Virus"]).strip()
        meta_map[acc] = virus
    return meta_map


def extract_accession_from_id(seq_id: str) -> str:
    """
    Extract accession from sequence ID.

    Example:
      'MN072625.1|LSDVgp002|inactivated|frame0' -> 'MN072625.1'
    """
    return seq_id.split("|")[0]


def is_terminal_stop(seq: str, pos0: int) -> bool:
    """
    Return True if a '*' at 0-based position pos0 should be considered
    the normal terminal stop of the protein and ignored.

    Logic:
      - Find the last non-gap character in the sequence.
      - If that character is '*' and its index == pos0, treat it as
        the normal terminal stop.
    """
    last_idx = None
    for i in range(len(seq) - 1, -1, -1):
        if seq[i] != "-":
            last_idx = i
            break
    if last_idx is None:
        return False
    return seq[last_idx] == "*" and last_idx == pos0


def build_disruptive_matrix_for_alignment(
    aln_path: str,
) -> Tuple[str, pd.DataFrame, Dict[str, int]]:
    """
    For one *_aa.fa alignment, build a presence/absence matrix for disruptive
    mutations ('*' = stop from SNP, '!' = frameshift) excluding the normal
    final stop codon.

    Returns:
      gene_name,
      matrix_df (columns: Sample, then mutation columns 0/1),
      summary_counts: dict with counts for this gene.
    """
    gene_name = os.path.basename(aln_path).replace("_aa.fa", "")
    records = read_fasta(aln_path)
    if not records:
        raise ValueError(f"No sequences found in {aln_path}")

    # Map: sample_id -> set of mutation IDs (initialize with all samples,
    # so even samples with no disruptive mutations are retained)
    sample_mutations: Dict[str, Set[str]] = defaultdict(set)
    for seq_id, _ in records:
        sample_mutations[seq_id] = set()

    # Summary counters
    n_stop_mut = 0
    n_fs_mut = 0
    positions_with_stop: Set[int] = set()
    positions_with_fs: Set[int] = set()

    for seq_id, seq in records:
        aln_len = len(seq)
        if aln_len == 0:
            continue

        for pos0 in range(aln_len):
            aa = seq[pos0]
            if aa not in ("*", "!"):
                continue

            # Exclude the normal terminal stop
            if aa == "*" and is_terminal_stop(seq, pos0):
                continue

            pos1 = pos0 + 1  # 1-based position in alignment
            if aa == "*":
                mut_type = "STOP"
                n_stop_mut += 1
                positions_with_stop.add(pos1)
            else:
                mut_type = "FSHIFT"
                n_fs_mut += 1
                positions_with_fs.add(pos1)

            mut_id = f"{mut_type}_{pos1}"
            sample_mutations[seq_id].add(mut_id)

    # Collect all mutation IDs for this gene
    all_mut_ids: List[str] = sorted(
        {m for muts in sample_mutations.values() for m in muts}
    )

    # Build presence/absence matrix
    # Order samples by group (if metadata is available later, we will sort
    # again in main; here we preserve the original order from the alignment).
    sample_order = [seq_id for seq_id, _ in records]

    rows: List[Dict[str, int]] = []
    for seq_id in sample_order:
        row: Dict[str, int] = {"Sample": seq_id}
        present = sample_mutations[seq_id]
        for mid in all_mut_ids:
            row[mid] = 1 if mid in present else 0
        rows.append(row)

    if all_mut_ids:
        cols = ["Sample"] + all_mut_ids
    else:
        cols = ["Sample"]

    matrix_df = pd.DataFrame(rows, columns=cols)

    summary_counts = {
        "gene": gene_name,
        "n_samples": len(records),
        "n_mutations_columns": len(all_mut_ids),
        "n_stop_events": n_stop_mut,
        "n_frameshift_events": n_fs_mut,
        "n_positions_with_stop": len(positions_with_stop),
        "n_positions_with_frameshift": len(positions_with_fs),
    }

    return gene_name, matrix_df, summary_counts


def summarize_shared_by_group(
    gene: str,
    matrix_df: pd.DataFrame,
    metadata_map: Dict[str, str],
) -> pd.DataFrame:
    """
    For one gene, identify disruptive mutations that are:
      - present in all SPPV samples
      - present in all GTPV samples
      - present in all SPPV and all GTPV samples
    based on the presence/absence matrix and metadata.
    """
    if matrix_df.shape[0] == 0:
        return pd.DataFrame()

    feature_cols = [c for c in matrix_df.columns if c != "Sample"]
    if not feature_cols:
        return pd.DataFrame()

    # Map samples to groups based on accession
    sppv_samples = []
    gtpv_samples = []
    for sid in matrix_df["Sample"]:
        acc = extract_accession_from_id(str(sid))
        g = metadata_map.get(acc, "Unknown")
        if g == "SPPV":
            sppv_samples.append(sid)
        elif g == "GTPV":
            gtpv_samples.append(sid)

    results = []
    for mut in feature_cols:
        col = matrix_df.set_index("Sample")[mut]

        n_sppv_total = len(sppv_samples)
        n_gtpv_total = len(gtpv_samples)

        if n_sppv_total > 0:
            sppv_vals = col.loc[sppv_samples]
            sppv_all = bool((sppv_vals == 1).all())
            sppv_with = int((sppv_vals == 1).sum())
        else:
            sppv_all = False
            sppv_with = 0

        if n_gtpv_total > 0:
            gtpv_vals = col.loc[gtpv_samples]
            gtpv_all = bool((gtpv_vals == 1).all())
            gtpv_with = int((gtpv_vals == 1).sum())
        else:
            gtpv_all = False
            gtpv_with = 0

        results.append(
            {
                "Gene": gene,
                "Mutation": mut,
                "SPPV_all": sppv_all,
                "SPPV_with_mutation": sppv_with,
                "SPPV_total": n_sppv_total,
                "GTPV_all": gtpv_all,
                "GTPV_with_mutation": gtpv_with,
                "GTPV_total": n_gtpv_total,
                "GTPV_and_SPPV_all": bool(sppv_all and gtpv_all),
            }
        )

    return pd.DataFrame(results)


def run_pca_on_matrix(
    gene: str,
    matrix_df: pd.DataFrame,
    metadata_map: Dict[str, str],
    out_dir: str,
) -> None:
    """
    Run PCA on the mutation presence/absence matrix and save coordinates.

    Output:
      - <gene>_disruptive_PCA.csv with columns:
          Sample, Group, PC1, PC2, and optionally the original mutation columns.
    """
    if matrix_df.shape[0] == 0:
        return

    feature_cols = [c for c in matrix_df.columns if c != "Sample"]
    if not feature_cols:
        return

    X = matrix_df[feature_cols].values.astype(float)

    # Drop constant columns (no variation)
    col_var = X.var(axis=0)
    keep_idx = np.where(col_var > 0)[0]
    if keep_idx.size == 0:
        return
    X = X[:, keep_idx]

    # Center data
    X_mean = X.mean(axis=0, keepdims=True)
    Xc = X - X_mean

    # PCA via SVD
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    scores = U * S  # rows = samples

    pc1 = scores[:, 0]
    pc2 = scores[:, 1] if scores.shape[1] > 1 else np.zeros_like(pc1)

    pca_df = matrix_df.copy()
    pca_df["PC1"] = pc1
    pca_df["PC2"] = pc2

    # Attach group from metadata
    groups: List[str] = []
    for sid in pca_df["Sample"]:
        acc = extract_accession_from_id(str(sid))
        groups.append(metadata_map.get(acc, "Unknown"))
    pca_df["Group"] = groups

    out_path = os.path.join(out_dir, f"{gene}_disruptive_PCA.csv")
    pca_df.to_csv(out_path, index=False)


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    align_dir = os.path.join(base_dir, "align")

    metadata_file = os.path.join(base_dir, "Dating_parhcments_beast - Sheet12.tsv")
    if os.path.exists(metadata_file):
        metadata_map = load_metadata(metadata_file)
    else:
        metadata_map = {}

    aa_paths = sorted(glob.glob(os.path.join(align_dir, "*_aa.fa")))
    if not aa_paths:
        raise SystemExit(f"No *_aa.fa alignments found in {align_dir}")

    out_dir = os.path.join(base_dir, "disruptive_matrices")
    os.makedirs(out_dir, exist_ok=True)

    all_summaries: List[Dict[str, int]] = []
    shared_dfs = []

    for path in aa_paths:
        gene, matrix_df, summary = build_disruptive_matrix_for_alignment(path)

        # If metadata is available, order samples in the matrix by group for readability
        if metadata_map and not matrix_df.empty:
            def group_key(sid: str) -> str:
                acc = extract_accession_from_id(str(sid))
                return metadata_map.get(acc, "Unknown")

            ordered_samples = sorted(
                matrix_df["Sample"].tolist(),
                key=lambda s: (group_key(str(s)), str(s)),
            )
            matrix_df = matrix_df.set_index("Sample").loc[ordered_samples].reset_index()

        all_summaries.append(summary)

        # Save presence/absence matrix for this gene (including samples with no mutations)
        matrix_out = os.path.join(out_dir, f"{gene}_disruptive_matrix.csv")
        matrix_df.to_csv(matrix_out, index=False)

        # Run PCA if we have metadata
        if metadata_map:
            run_pca_on_matrix(gene, matrix_df, metadata_map, out_dir)

            # Also compute which disruptive mutations are shared by all samples
            # in GTPV and/or SPPV for this gene
            shared_df = summarize_shared_by_group(gene, matrix_df, metadata_map)
            if not shared_df.empty:
                shared_gene_out = os.path.join(
                    out_dir, f"{gene}_disruptive_shared_by_group.csv"
                )
                shared_df.to_csv(shared_gene_out, index=False)
                shared_dfs.append(shared_df)

    # Save per-gene summary statistics for disruptive mutations
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        summary_out = os.path.join(out_dir, "disruptive_summary_per_gene.csv")
        summary_df.to_csv(summary_out, index=False)

    # Save combined summary of shared mutations across all genes
    if shared_dfs:
        combined_shared = pd.concat(shared_dfs, ignore_index=True)
        shared_out = os.path.join(out_dir, "disruptive_shared_by_group_all_genes.csv")
        combined_shared.to_csv(shared_out, index=False)


if __name__ == "__main__":
    main()


