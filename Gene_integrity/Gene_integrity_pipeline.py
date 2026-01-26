#!/usr/bin/env python3
"""
Automate Liftoff annotation transfer and gene extraction across multiple genomes.

Example:
    python liftoff_pipeline.py \
        --reference-fasta LSDV.fna \
        --reference-gff LSDV.gff \
        --targets-dir GENOMES \
        --output-dir results \
        --liftoff-bin liftoff \
        --threads 8
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
import shlex


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Lift annotations with Liftoff and collate genes by locus tag."
    )
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--reference-gff", required=True, type=Path)
    parser.add_argument("--targets-dir", required=True, type=Path)
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="Directory for Liftoff results."
    )
    parser.add_argument(
        "--liftoff-bin",
        default="liftoff",
        help="Path to Liftoff executable (default: liftoff in $PATH).",
    )
    parser.add_argument(
        "--threads", type=int, default=4, help="Threads to pass to Liftoff (-p)."
    )
    parser.add_argument(
        "--feature-type",
        default="gene",
        choices=["gene", "CDS"],
        help="GFF feature to extract sequences for (default: gene).",
    )
    parser.add_argument(
        "--upstream-extension",
        type=int,
        default=0,
        help="Number of nucleotides to extend upstream from Liftoff coordinates to search for start codons (default: 0).",
    )
    parser.add_argument(
        "--start-window",
        type=int,
        default=20,
        help="Number of codons (3x nucleotides) to search downstream from gene start for ATG start codon (default: 20, i.e., 60 nucleotides).",
    )
    parser.add_argument(
        "--liftoff-min-alignment",
        type=float,
        default=0.4,
        help="Value to pass to Liftoff -a (minimum aligned fraction, default 0.4).",
    )
    parser.add_argument(
        "--liftoff-min-identity",
        type=float,
        default=0.4,
        help="Value to pass to Liftoff -s (minimum sequence identity, default 0.4).",
    )
    parser.add_argument(
        "--liftoff-extra",
        default="",
        help="Additional arguments appended verbatim to the Liftoff command.",
    )
    parser.add_argument(
        "--liftoff-quiet",
        action="store_true",
        help="Add -quiet to Liftoff (requires Liftoff >=1.6).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-run Liftoff even if output GFF already exists.",
    )
    parser.add_argument(
        "--missing-report",
        type=Path,
        help="Optional path for presence/absence report. Defaults to OUTPUT_DIR/missing_annotations.tsv",
    )
    parser.add_argument(
        "--blast-rescue",
        action="store_true",
        help="Attempt to recover missing loci via BLASTn searches.",
    )
    parser.add_argument(
        "--blastn-bin",
        default="blastn",
        help="Path to BLASTn executable (default: blastn in $PATH).",
    )
    parser.add_argument(
        "--makeblastdb-bin",
        default="makeblastdb",
        help="Path to makeblastdb executable (default: makeblastdb in $PATH).",
    )
    parser.add_argument(
        "--blast-min-pident",
        type=float,
        default=40.0,
        help="Minimum percent identity for reporting BLAST rescue hits (default 40).",
    )
    parser.add_argument(
        "--blast-min-coverage",
        type=float,
        default=0.3,
        help="Minimum query coverage (alignment length / gene length) for BLAST hits (default 0.3).",
    )
    parser.add_argument(
        "--blast-report",
        type=Path,
        help="Optional path for BLAST rescue hits table. Defaults to OUTPUT_DIR/blast_rescue_hits.tsv",
    )
    parser.add_argument(
        "--annotation-status-report",
        type=Path,
        help="Optional path for per-sample locus status table. Defaults to OUTPUT_DIR/final_annotation_status.tsv",
    )
    parser.add_argument(
        "--skip-gene-status",
        action="store_true",
        help="Skip translation-based gene activity summary generation.",
    )
    parser.add_argument(
        "--run-3frame-analysis",
        action="store_true",
        help="Run parallel 3-frame translation analysis (independent of Liftoff coordinates).",
    )
    parser.add_argument(
        "--3frame-summary",
        type=Path,
        dest="frame3_summary",
        help="Path for 3-frame analysis summary (default: OUTPUT_DIR/gene_activity_3frame_summary.tsv).",
    )
    parser.add_argument(
        "--3frame-protein-dir",
        type=Path,
        dest="frame3_protein_dir",
        help="Directory for 3-frame protein FASTAs (default: OUTPUT_DIR/per_locus_proteins_3frame).",
    )
    parser.add_argument(
        "--3frame-matrix",
        type=Path,
        dest="frame3_matrix",
        help="Path for 3-frame activity matrix (default: OUTPUT_DIR/gene_activity_3frame_matrix.tsv).",
    )
    parser.add_argument(
        "--3frame-heatmap",
        type=Path,
        dest="frame3_heatmap",
        help="Path for 3-frame heatmap (default: OUTPUT_DIR/gene_activity_3frame_heatmap.png).",
    )
    parser.add_argument(
        "--gene-summary",
        type=Path,
        help="Path for gene activity summary TSV (default OUTPUT_DIR/gene_activity_summary.tsv).",
    )
    parser.add_argument(
        "--protein-fastas-dir",
        type=Path,
        help="Directory for per-locus protein FASTAs (default OUTPUT_DIR/per_locus_proteins).",
    )
    parser.add_argument(
        "--activity-matrix",
        type=Path,
        help="Path for gene activity matrix TSV (default OUTPUT_DIR/gene_activity_matrix.tsv).",
    )
    parser.add_argument(
        "--activity-heatmap",
        type=Path,
        help="Path for gene activity heatmap image (default OUTPUT_DIR/gene_activity_heatmap.png).",
    )
    return parser.parse_args()


CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def find_target_fastas(target_dir: Path) -> list[Path]:
    fastas = sorted(
        path
        for path in target_dir.iterdir()
        if path.suffix.lower() in {".fa", ".fasta", ".fna"} or path.name.endswith(".fa")
    )
    if not fastas:
        raise FileNotFoundError(f"No FASTA files found in {target_dir}")
    return fastas


def run_liftoff(
    liftoff_bin: str,
    target_fasta: Path,
    reference_fasta: Path,
    reference_gff: Path,
    output_gff: Path,
    threads: int,
    min_alignment: float,
    min_identity: float,
    quiet: bool,
    extra_args: str,
) -> None:
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        liftoff_bin,
        "-g",
        str(reference_gff),
        "-o",
        str(output_gff),
        "-p",
        str(threads),
        "-a",
        str(min_alignment),
        "-s",
        str(min_identity),
        str(target_fasta),
        str(reference_fasta),
    ]
    if quiet:
        cmd.append("-quiet")
    if extra_args:
        cmd.extend(shlex.split(extra_args))
    print(f"[liftoff] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def parse_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current_id = None
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            else:
                if current_id is None:
                    raise ValueError(f"FASTA missing header before sequence in {path}")
                sequences[current_id].append(line.upper())
    return {seq_id: "".join(parts) for seq_id, parts in sequences.items()}


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]


def parse_attributes(attr_field: str) -> dict[str, str]:
    attrs = {}
    for part in attr_field.split(";"):
        if not part.strip():
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def load_reference_loci(reference_gff: Path, feature_type: str) -> set[str]:
    loci: set[str] = set()
    with reference_gff.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            _, _, ftype, _, _, _, _, _, attrs = parts
            if ftype != feature_type:
                continue
            info = parse_attributes(attrs)
            locus_tag = info.get("locus_tag") or info.get("ID")
            if locus_tag:
                loci.add(locus_tag)
    if not loci:
        raise ValueError(f"No {feature_type} entries with locus tags found in {reference_gff}")
    return loci


def extract_reference_feature_sequences(
    reference_gff: Path, reference_sequences: dict[str, str], feature_type: str
) -> dict[str, str]:
    features: dict[str, str] = {}
    with reference_gff.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = line.strip().split("\t")
            if ftype != feature_type:
                continue
            info = parse_attributes(attrs)
            locus_tag = info.get("locus_tag") or info.get("ID")
            if not locus_tag:
                continue
            chrom_seq = reference_sequences.get(chrom)
            if chrom_seq is None:
                raise KeyError(f"Reference sequence '{chrom}' missing in FASTA.")
            subseq = chrom_seq[int(start) - 1 : int(end)]
            if strand == "-":
                subseq = reverse_complement(subseq)
            features[locus_tag] = subseq
    if not features:
        raise ValueError(
            f"No {feature_type} sequences recovered from {reference_gff}. "
            "Ensure locus_tag or ID attributes are present."
        )
    return features


def extract_genes(
    gff_path: Path,
    fasta_sequences: dict[str, str],
    feature_type: str,
    sample_name: str,
    upstream_extension: int = 0,
) -> tuple[dict[str, tuple[str, str]], dict[str, dict]]:
    """
    Extract genes from GFF file.
    
    Returns:
        (genes_dict, metadata_dict) where:
        - genes_dict: {locus_tag: (header, sequence)}
        - metadata_dict: {locus_tag: {chrom, orig_start, orig_end, extracted_start, extracted_end, strand, upstream_ext}}
    """
    genes = {}
    metadata = {}
    with gff_path.open() as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = line.strip().split(
                "\t"
            )
            if ftype != feature_type:
                continue
            info = parse_attributes(attrs)
            locus_tag = info.get("locus_tag") or info.get("ID")
            if not locus_tag:
                continue
            seq = fasta_sequences.get(chrom)
            if seq is None:
                raise KeyError(f"Sequence {chrom} missing in FASTA for {sample_name}")
            
            # Original coordinates from GFF (1-based)
            orig_start = int(start)
            orig_end = int(end)
            
            # Calculate extracted coordinates (0-based for slicing)
            start_pos = orig_start - 1
            end_pos = orig_end
            # Extend upstream: for forward strand, extend to lower coordinates;
            # for reverse strand, extend to higher coordinates (will be reverse-complemented)
            if strand == "+":
                extracted_start = max(0, start_pos - upstream_extension)
                extracted_end = end_pos
            else:  # strand == "-"
                extracted_start = start_pos
                extracted_end = min(len(seq), end_pos + upstream_extension)
            
            subseq = seq[extracted_start : extracted_end]
            if strand == "-":
                subseq = reverse_complement(subseq)
            
            header = f"{sample_name}|{locus_tag}"
            genes[locus_tag] = (header, subseq)
            
            # Store metadata
            metadata[locus_tag] = {
                "chrom": chrom,
                "orig_start": orig_start,
                "orig_end": orig_end,
                "extracted_start": extracted_start + 1,  # Convert to 1-based for reporting
                "extracted_end": extracted_end,
                "strand": strand,
                "upstream_extension": upstream_extension,
                "sequence_length": len(subseq),
            }
    return genes, metadata


def write_locus_fastas(
    genes_by_locus: dict[str, list[tuple[str, str]]], output_dir: Path
) -> None:
    fasta_dir = output_dir / "per_locus_fastas"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for locus_tag, entries in genes_by_locus.items():
        locus_file = fasta_dir / f"{locus_tag}.fa"
        with locus_file.open("w") as out:
            for header, sequence in entries:
                out.write(f">{header}\n")
                for i in range(0, len(sequence), 60):
                    out.write(sequence[i : i + 60] + "\n")


def write_raw_extracted_sequences(
    all_metadata: dict[str, dict[str, dict]], 
    genes_by_locus: dict[str, list[tuple[str, str]]],
    output_dir: Path
) -> None:
    """
    Write raw extracted sequences with full metadata for debugging.
    
    Args:
        all_metadata: {sample_name: {locus_tag: metadata_dict}}
        genes_by_locus: {locus_tag: [(header, sequence), ...]}
        output_dir: Output directory
    """
    raw_dir = output_dir / "raw_extracted_sequences"
    raw_dir.mkdir(parents=True, exist_ok=True)
    
    # Write per-sample metadata
    for sample_name, sample_metadata in all_metadata.items():
        sample_file = raw_dir / f"{sample_name}_extraction_metadata.tsv"
        with sample_file.open("w") as out:
            out.write(
                "locus_tag\tchrom\torig_start\torig_end\textracted_start\textracted_end\t"
                "strand\tupstream_extension\tsequence_length\n"
            )
            for locus_tag, meta in sorted(sample_metadata.items()):
                out.write(
                    f"{locus_tag}\t{meta['chrom']}\t{meta['orig_start']}\t{meta['orig_end']}\t"
                    f"{meta['extracted_start']}\t{meta['extracted_end']}\t{meta['strand']}\t"
                    f"{meta['upstream_extension']}\t{meta['sequence_length']}\n"
                )
    
    # Write per-locus summary
    summary_file = raw_dir / "extraction_summary.tsv"
    all_loci = set()
    for sample_metadata in all_metadata.values():
        all_loci.update(sample_metadata.keys())
    
    with summary_file.open("w") as out:
        out.write("locus_tag\tsample\tchrom\torig_start\torig_end\textracted_start\textracted_end\tstrand\tupstream_ext\tlength\n")
        for locus_tag in sorted(all_loci):
            for sample_name in sorted(all_metadata.keys()):
                if locus_tag in all_metadata[sample_name]:
                    meta = all_metadata[sample_name][locus_tag]
                    out.write(
                        f"{locus_tag}\t{sample_name}\t{meta['chrom']}\t{meta['orig_start']}\t{meta['orig_end']}\t"
                        f"{meta['extracted_start']}\t{meta['extracted_end']}\t{meta['strand']}\t"
                        f"{meta['upstream_extension']}\t{meta['sequence_length']}\n"
                    )
    
    # Write raw sequences as FASTA files (per sample)
    fasta_dir = raw_dir / "fasta_sequences"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a mapping from (sample, locus) to sequence
    sample_locus_sequences: dict[tuple[str, str], str] = {}
    for locus_tag, entries in genes_by_locus.items():
        for header, sequence in entries:
            # Parse sample name from header (format: sample_name|locus_tag)
            parts = header.split("|")
            if len(parts) >= 1:
                sample_name = parts[0]
                sample_locus_sequences[(sample_name, locus_tag)] = sequence
    
    # Write per-sample FASTA files
    for sample_name in sorted(all_metadata.keys()):
        sample_fasta = fasta_dir / f"{sample_name}_raw_sequences.fa"
        with sample_fasta.open("w") as out:
            for locus_tag in sorted(all_metadata[sample_name].keys()):
                seq = sample_locus_sequences.get((sample_name, locus_tag))
                if seq:
                    meta = all_metadata[sample_name][locus_tag]
                    # Create detailed header with metadata
                    header = (
                        f">{sample_name}|{locus_tag}|"
                        f"chrom={meta['chrom']}|"
                        f"orig_start={meta['orig_start']}|orig_end={meta['orig_end']}|"
                        f"extracted_start={meta['extracted_start']}|extracted_end={meta['extracted_end']}|"
                        f"strand={meta['strand']}|upstream_ext={meta['upstream_extension']}|"
                        f"length={meta['sequence_length']}"
                    )
                    out.write(f"{header}\n")
                    # Write sequence in 60-character lines
                    for i in range(0, len(seq), 60):
                        out.write(seq[i : i + 60] + "\n")
    
    # Also write per-locus FASTA files (all samples for each locus)
    per_locus_fasta_dir = raw_dir / "per_locus_fasta"
    per_locus_fasta_dir.mkdir(parents=True, exist_ok=True)
    
    for locus_tag in sorted(all_loci):
        locus_fasta = per_locus_fasta_dir / f"{locus_tag}_raw.fa"
        with locus_fasta.open("w") as out:
            for sample_name in sorted(all_metadata.keys()):
                if locus_tag in all_metadata[sample_name]:
                    seq = sample_locus_sequences.get((sample_name, locus_tag))
                    if seq:
                        meta = all_metadata[sample_name][locus_tag]
                        header = (
                            f">{sample_name}|{locus_tag}|"
                            f"chrom={meta['chrom']}|"
                            f"orig_start={meta['orig_start']}|orig_end={meta['orig_end']}|"
                            f"extracted_start={meta['extracted_start']}|extracted_end={meta['extracted_end']}|"
                            f"strand={meta['strand']}|upstream_ext={meta['upstream_extension']}|"
                            f"length={meta['sequence_length']}"
                        )
                        out.write(f"{header}\n")
                        for i in range(0, len(seq), 60):
                            out.write(seq[i : i + 60] + "\n")
    
    print(f"Raw extraction metadata and FASTA sequences written to {raw_dir}/")


def write_missing_report(
    missing_by_locus: dict[str, set[str]], report_path: Path
) -> None:
    remaining = {locus: sorted(samples) for locus, samples in missing_by_locus.items() if samples}
    if not remaining:
        report_path.write_text("# All reference loci recovered in every sample.\n")
        return
    with report_path.open("w") as handle:
        handle.write("locus_tag\tmissing_samples\tmissing_count\n")
        for locus_tag in sorted(remaining):
            samples = remaining[locus_tag]
            handle.write(f"{locus_tag}\t{','.join(samples)}\t{len(samples)}\n")


def write_annotation_status(
    annotation_status: dict[str, dict[str, str]], report_path: Path
) -> None:
    with report_path.open("w") as handle:
        handle.write("sample\tlocus_tag\tstatus\n")
        for sample in sorted(annotation_status):
            for locus_tag in sorted(annotation_status[sample]):
                handle.write(f"{sample}\t{locus_tag}\t{annotation_status[sample][locus_tag]}\n")


def ensure_blast_db(fasta_path: Path, db_dir: Path, makeblastdb_bin: str) -> Path:
    db_dir.mkdir(parents=True, exist_ok=True)
    prefix = db_dir / fasta_path.stem
    if prefix.with_suffix(".nhr").exists():
        return prefix
    cmd = [
        makeblastdb_bin,
        "-in",
        str(fasta_path),
        "-dbtype",
        "nucl",
        "-parse_seqids",
        "-out",
        str(prefix),
    ]
    print(f"[makeblastdb] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return prefix


def run_blastn_for_locus(
    locus_tag: str,
    gene_sequence: str,
    sample_name: str,
    db_prefix: Path,
    blastn_bin: str,
    min_pident: float,
    min_coverage: float,
    workdir: Path,
) -> dict | None:
    workdir.mkdir(parents=True, exist_ok=True)
    query_path = workdir / f"{sample_name}_{locus_tag}.fa"
    with query_path.open("w") as handle:
        handle.write(f">{locus_tag}\n")
        for i in range(0, len(gene_sequence), 60):
            handle.write(gene_sequence[i : i + 60] + "\n")
    out_path = workdir / f"{sample_name}_{locus_tag}.tsv"
    cmd = [
        blastn_bin,
        "-query",
        str(query_path),
        "-db",
        str(db_prefix),
        "-outfmt",
        "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
        "-max_target_seqs",
        "5",
        "-max_hsps",
        "1",
    ]
    print(f"[blastn] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    best_hit = None
    if out_path.exists():
        with out_path.open() as handle:
            for line in handle:
                parts = line.strip().split("\t")
                if len(parts) != 10:
                    continue
                qseqid, sseqid, pident, length, qstart, qend, sstart, send, evalue, bitscore = (
                    parts
                )
                pident = float(pident)
                length = float(length)
                coverage = length / len(gene_sequence)
                if pident >= min_pident and coverage >= min_coverage:
                    best_hit = {
                        "locus_tag": qseqid,
                        "sample": sample_name,
                        "subject": sseqid,
                        "pident": pident,
                        "coverage": coverage,
                        "length": int(length),
                        "qstart": int(float(qstart)),
                        "qend": int(float(qend)),
                        "sstart": int(float(sstart)),
                        "send": int(float(send)),
                        "evalue": float(evalue),
                        "bitscore": float(bitscore),
                    }
                    break
    query_path.unlink(missing_ok=True)
    out_path.unlink(missing_ok=True)
    return best_hit


def run_blast_rescue(
    missing_by_sample: dict[str, set[str]],
    reference_sequences: dict[str, str],
    sample_fastas: dict[str, Path],
    output_dir: Path,
    blastn_bin: str,
    makeblastdb_bin: str,
    min_pident: float,
    min_coverage: float,
) -> list[dict]:
    db_dir = output_dir / "blast_dbs"
    tmp_dir = output_dir / "blast_tmp"
    hits: list[dict] = []
    for sample_name, loci in missing_by_sample.items():
        if not loci:
            continue
        target_fasta = sample_fastas[sample_name]
        db_prefix = ensure_blast_db(target_fasta, db_dir, makeblastdb_bin)
        for locus_tag in loci:
            gene_seq = reference_sequences.get(locus_tag)
            if not gene_seq:
                continue
            hit = run_blastn_for_locus(
                locus_tag,
                gene_seq,
                sample_name,
                db_prefix,
                blastn_bin,
                min_pident,
                min_coverage,
                tmp_dir,
            )
            if hit:
                hits.append(hit)
    return hits


def extract_blast_hit_sequence(
    hit: dict, sample_sequences: dict[str, dict[str, str]]
) -> str | None:
    sample = hit["sample"]
    subject = hit["subject"]
    sequences = sample_sequences.get(sample)
    if sequences is None:
        return None
    contig = sequences.get(subject)
    if contig is None:
        return None
    start = hit["sstart"]
    end = hit["send"]
    if start <= end:
        subseq = contig[start - 1 : end]
    else:
        subseq = reverse_complement(contig[end - 1 : start])
    return subseq


def write_blast_hits(hits: list[dict], path: Path) -> None:
    if not hits:
        path.write_text("# No BLAST rescue hits met the thresholds.\n")
        return
    with path.open("w") as handle:
        handle.write(
            "locus_tag\tsample\tsubject\tpident\tcoverage\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
        )
        for hit in hits:
            handle.write(
                "{locus_tag}\t{sample}\t{subject}\t{pident:.2f}\t{coverage:.2f}\t{length}\t"
                "{qstart}\t{qend}\t{sstart}\t{send}\t{evalue:.2e}\t{bitscore:.1f}\n".format(
                    **hit
                )
            )


def parse_per_locus_fastas(locus_dir: Path) -> dict[tuple[str, str], str]:
    sequences: dict[tuple[str, str], str] = {}
    if not locus_dir.exists():
        return sequences
    for fasta_path in sorted(locus_dir.glob("*.fa*")):
        locus_tag = fasta_path.stem
        with fasta_path.open() as handle:
            current_sample = None
            seq_chunks: list[str] = []
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_sample is not None:
                        sequences[(current_sample, locus_tag)] = "".join(seq_chunks).upper()
                    header = line[1:]
                    parts = header.split("|")
                    current_sample = parts[0]
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
            if current_sample is not None:
                sequences[(current_sample, locus_tag)] = "".join(seq_chunks).upper()
    return sequences


def clean_sequence(seq: str) -> str:
    return seq.upper().replace("U", "T").replace(" ", "")


def evaluate_cds(sequence: str, start_window: int = 20, upstream_extension: int = 0) -> tuple[str, str, str, str]:
    """
    Evaluate CDS and return (status, protein, reason, trimmed_sequence).
    trimmed_sequence is the sequence actually used for translation (after trimming to ATG if found upstream).
    
    Process:
    1. First, search for ATG in the first start_window codons (60 nucleotides) starting from upstream_extension
       (This searches the actual gene region, not the upstream extension)
    2. If no ATG found, search 30 nucleotides upstream from position 0
       (This searches the upstream extension region)
    """
    seq = clean_sequence(sequence)
    original_seq = seq
    
    # STEP 1: First, search for ATG in the first 60 nt of the ACTUAL GENE (not upstream region)
    # If upstream_extension was used, the actual gene starts at position upstream_extension
    # So we search from upstream_extension to upstream_extension + 60
    atg_offset = None
    gene_start_pos = upstream_extension  # Where the actual gene starts (after upstream extension)
    search_end_pos = gene_start_pos + (start_window * 3)  # 60 nt downstream from gene start
    
    # Search in the actual gene region first (positions gene_start_pos to search_end_pos)
    for pos in range(gene_start_pos, min(search_end_pos, len(seq) - 2), 3):  # Step by 3 (codon boundaries)
        codon = seq[pos : pos + 3]
        if codon == "ATG":
            atg_offset = pos
            break
    
    # STEP 2: If no ATG found in the gene region, search 30 nt upstream (positions 0 to upstream_extension)
    if atg_offset is None and upstream_extension > 0:
        upstream_search_limit = min(30, upstream_extension)  # Search up to 30 nt upstream, or less if extension is smaller
        for pos in range(0, upstream_search_limit - 2, 3):  # Step by 3 (codon boundaries)
            codon = seq[pos : pos + 3]
            if codon == "ATG":
                atg_offset = pos
                break
        
        # If upstream_extension > 30, also check the rest of the upstream region
        if atg_offset is None and upstream_extension > 30:
            for pos in range(30, min(upstream_extension, len(seq)) - 2, 3):
                codon = seq[pos : pos + 3]
                if codon == "ATG":
                    atg_offset = pos
                    break
    
    # STEP 3: If still no ATG and no upstream_extension was used, search from position 0
    # (This handles cases where sequence wasn't extended but we still want to check the beginning)
    if atg_offset is None and upstream_extension == 0:
        # Search first 30 nt as fallback
        for pos in range(0, min(30, len(seq) - 2), 3):
            codon = seq[pos : pos + 3]
            if codon == "ATG":
                atg_offset = pos
                break
    
    # Trim sequence to start from ATG if found
    trimmed_seq = seq
    if atg_offset is not None and atg_offset > 0:
        trimmed_seq = seq[atg_offset:]
    elif atg_offset is None:
        # No ATG found - use original sequence
        trimmed_seq = seq
    
    codon_count = len(trimmed_seq) // 3
    if codon_count == 0:
        return "unknown", "", "too_short_for_translation", trimmed_seq
    has_start = (trimmed_seq[:3] == "ATG")
    protein: list[str] = []
    stop_index = None
    ambiguous_codons = []
    for idx in range(codon_count):
        codon = trimmed_seq[idx * 3 : idx * 3 + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            # Handle ambiguous nucleotides (N, etc.) by translating as X
            if "N" in codon or any(base not in "ACGT" for base in codon):
                aa = "X"
                ambiguous_codons.append(idx)
            else:
                return "unknown", "", f"ambiguous_codon_{codon}", trimmed_seq
        if aa == "*":
            stop_index = idx
            break
        protein.append(aa)
    if not has_start:
        reason = "missing_start_codon"
        if ambiguous_codons:
            reason += f"_with_{len(ambiguous_codons)}_ambiguous_codons"
        return "inactivated", "".join(protein), reason, trimmed_seq
    if stop_index is None:
        # No stop codon detected
        if ambiguous_codons:
            reason = f"no_stop_detected_with_{len(ambiguous_codons)}_ambiguous_codons"
            return "likely-active", "".join(protein), reason, trimmed_seq
        else:
            return "active", "".join(protein), "no_stop_detected", trimmed_seq
    # Stop codon found - check if it's within last 5 AA (then it's fine)
    remaining = codon_count - stop_index - 1
    if remaining > 5:
        # Early stop (more than 5 AA before end)
        reason = f"early_stop_{remaining}_aa_before_end"
        if ambiguous_codons:
            reason += f"_with_{len(ambiguous_codons)}_ambiguous_codons"
        return "inactivated", "".join(protein), reason, trimmed_seq
    # Stop is within last 5 AA or exactly at position for 5 AA - this is fine (active)
    if ambiguous_codons:
        reason = f"stop_within_terminal_window_with_{len(ambiguous_codons)}_ambiguous_codons"
        return "likely-active", "".join(protein), reason, trimmed_seq
    else:
        return "active", "".join(protein), "stop_within_terminal_window", trimmed_seq


def find_best_reading_frame(sequence: str, start_window: int = 20, upstream_extension: int = 0) -> tuple[int, str, str, str]:
    """
    Try all 3 reading frames and find the best one.
    Returns: (frame_offset, translated_seq, status, reason)
    Frame 0 = no offset, Frame 1 = offset by 1, Frame 2 = offset by 2
    
    For each frame, searches for ATG in first 60 nt, then 30 nt upstream if not found.
    """
    seq = clean_sequence(sequence)
    best_frame = 0
    best_status = "unknown"
    best_protein = ""
    best_reason = "no_valid_frame"
    best_length = 0
    
    for frame_offset in range(3):
        frame_seq = seq[frame_offset:]
        if len(frame_seq) < 3:
            continue
        
        # Pass upstream_extension so evaluate_cds knows where the actual gene starts
        # (sequences in per_locus_fastas already include upstream extension if it was used)
        # Adjust upstream_extension for frame offset
        effective_upstream = max(0, upstream_extension - frame_offset) if upstream_extension > 0 else 0
        status, protein, reason, _ = evaluate_cds(frame_seq, start_window, upstream_extension=effective_upstream)
        
        # Score: prefer active > likely-active > inactivated > unknown
        # Among same status, prefer longer proteins
        score = 0
        if status == "active":
            score = 1000 + len(protein)
        elif status == "likely-active":
            score = 750 + len(protein)
        elif status == "inactivated":
            score = 500 + len(protein)
        elif status == "unknown":
            score = len(protein)
        
        if score > best_length or (score == best_length and status in ("active", "likely-active")):
            best_frame = frame_offset
            best_status = status
            best_protein = protein
            best_reason = f"frame_{frame_offset}_{reason}"
            best_length = score
    
    return best_frame, best_protein, best_status, best_reason


def evaluate_cds_3frame(sequence: str | None, start_window: int = 20, upstream_extension: int = 0) -> tuple[str, str, str]:
    """
    Evaluate CDS by trying all 3 reading frames and picking the best.
    This is independent of Liftoff coordinates.
    """
    if sequence is None or len(sequence) < 3:
        return "not_found", "", "not_found"
    
    frame_offset, protein, status, reason = find_best_reading_frame(sequence, start_window, upstream_extension)
    return status, protein, reason


def determine_gene_activity(
    base_status: str, sequence: str | None, start_window: int, upstream_extension: int = 0
) -> tuple[str, str, str, str]:
    """
    Determine gene activity and return (status, protein, reason, trimmed_sequence).
    """
    if base_status == "missing" or sequence is None:
        return "not_found", "", "not_found", ""
    return evaluate_cds(sequence, start_window, upstream_extension)


def status_to_matrix_symbol(status: str) -> str:
    if status == "active":
        return "1"
    if status == "likely-active":
        return "3"
    if status == "inactivated":
        return "0"
    if status == "not_found":
        return "Not found"
    return "NA"


def status_to_numeric(status: str) -> float:
    mapping = {"active": 1.0, "likely-active": 2.0, "inactivated": 0.0, "not_found": -1.0}
    return mapping.get(status, float("nan"))


def write_protein_fastas(
    protein_dir: Path, protein_entries: dict[str, list[tuple[str, str]]]
) -> None:
    protein_dir.mkdir(parents=True, exist_ok=True)
    for locus_tag, entries in protein_entries.items():
        if not entries:
            continue
        out_path = protein_dir / f"{locus_tag}.faa"
        with out_path.open("w") as handle:
            for header, seq in entries:
                handle.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    handle.write(seq[i : i + 60] + "\n")


def write_activity_matrix(
    matrix_path: Path,
    samples: list[str],
    loci: list[str],
    matrix_status: dict[str, dict[str, str]],
) -> None:
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    with matrix_path.open("w") as handle:
        handle.write("sample\t" + "\t".join(loci) + "\n")
        for sample in samples:
            values = [
                status_to_matrix_symbol(matrix_status.get(sample, {}).get(locus, "not_found"))
                for locus in loci
            ]
            handle.write(sample + "\t" + "\t".join(values) + "\n")


def render_activity_heatmap(
    heatmap_path: Path,
    samples: list[str],
    loci: list[str],
    matrix_status: dict[str, dict[str, str]],
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    data = np.array(
        [
            [status_to_numeric(matrix_status.get(sample, {}).get(locus, "not_found")) for locus in loci]
            for sample in samples
        ]
    )
    # Colors: red (inactivated), yellow (likely-active), green (active)
    cmap = ListedColormap(["#d73027", "#fee08b", "#abdda4", "#1a9850"])
    cmap.set_bad("#bdbdbd")
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(figsize=(0.25 * len(loci) + 2, 0.3 * len(samples) + 2))
    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(loci)))
    ax.set_xticklabels(loci, rotation=90, fontsize=6)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(samples, fontsize=6)
    cbar = plt.colorbar(im, ticks=[-1, 0, 1, 2], fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(["Not found", "Inactivated", "Active", "Likely-active"])
    cbar.ax.tick_params(labelsize=8)
    plt.tight_layout()
    heatmap_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(heatmap_path, dpi=300)
    plt.close(fig)


def run_gene_activity_analysis(
    per_locus_dir: Path,
    annotation_status: dict[str, dict[str, str]],
    loci: list[str],
    summary_path: Path,
    protein_dir: Path,
    matrix_path: Path,
    heatmap_path: Path,
    start_window: int,
    upstream_extension: int = 0,
) -> None:
    sequences = parse_per_locus_fastas(per_locus_dir)
    samples = sorted(annotation_status)
    loci_sorted = sorted(loci)
    protein_entries: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)
    matrix_status: defaultdict[str, dict[str, str]] = defaultdict(dict)
    # Store trimmed nucleotide sequences used for translation
    trimmed_nucleotide_entries: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w") as out:
        out.write(
            "sample\tlocus_tag\tgene_status\tprotein_length\tprotein_sequence\treason\n"
        )
        for sample in samples:
            for locus in loci_sorted:
                base_state = annotation_status[sample].get(locus, "missing")
                seq = sequences.get((sample, locus))
                status, protein_seq, reason, trimmed_seq = determine_gene_activity(
                    base_state, seq, start_window, upstream_extension
                )
                matrix_status[sample][locus] = status
                protein_len = len(protein_seq)
                out.write(
                    f"{sample}\t{locus}\t{status}\t{protein_len}\t{protein_seq}\t{reason}\n"
                )
                if protein_seq:
                    header = f"{sample}|{locus}|{status}"
                    protein_entries[locus].append((header, protein_seq))
                # Store trimmed nucleotide sequence if available
                if trimmed_seq:
                    header = f"{sample}|{locus}|{status}"
                    trimmed_nucleotide_entries[locus].append((header, trimmed_seq))

    write_protein_fastas(protein_dir, protein_entries)
    # Write trimmed nucleotide sequences used for translation
    trimmed_nuc_dir = protein_dir.parent / "per_locus_trimmed_sequences"
    write_protein_fastas(trimmed_nuc_dir, trimmed_nucleotide_entries)
    write_activity_matrix(matrix_path, samples, loci_sorted, matrix_status)
    render_activity_heatmap(heatmap_path, samples, loci_sorted, matrix_status)
    print(f"Gene activity summary written to {summary_path}")
    print(f"Protein FASTAs written to {protein_dir}")
    print(f"Trimmed nucleotide sequences written to {trimmed_nuc_dir}")
    print(f"Activity matrix written to {matrix_path}")
    print(f"Heatmap saved to {heatmap_path}")


def run_3frame_analysis(
    per_locus_dir: Path,
    loci: list[str],
    target_fastas: list[Path],
    summary_path: Path,
    protein_dir: Path,
    matrix_path: Path,
    heatmap_path: Path,
    start_window: int,
    upstream_extension: int = 0,
) -> None:
    """
    Run 3-frame translation analysis independently of Liftoff coordinates.
    This analyzes nucleotide sequences directly by trying all 3 reading frames.
    """
    sequences = parse_per_locus_fastas(per_locus_dir)
    samples = sorted({fasta.stem for fasta in target_fastas})
    loci_sorted = sorted(loci)
    protein_entries: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)
    matrix_status: defaultdict[str, dict[str, str]] = defaultdict(dict)
    trimmed_nucleotide_entries: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)
    
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w") as out:
        out.write(
            "sample\tlocus_tag\tgene_status\tprotein_length\tprotein_sequence\treason\tframe_offset\n"
        )
        for sample in samples:
            for locus in loci_sorted:
                seq = sequences.get((sample, locus))
                if seq is None:
                    status, protein_seq, reason = "not_found", "", "not_found"
                    frame_offset = -1
                    trimmed_seq = ""
                else:
                    frame_offset, protein_seq, status, reason = find_best_reading_frame(seq, start_window, upstream_extension)
                    # Get the trimmed sequence from the best frame
                    if frame_offset >= 0:
                        frame_seq = clean_sequence(seq)[frame_offset:]
                        # Evaluate to get trimmed sequence
                        # Note: upstream_extension applies to the original sequence, so we need to adjust
                        # for the frame offset. If we're in frame 1 or 2, the upstream region shifts.
                        effective_upstream = max(0, upstream_extension - frame_offset) if upstream_extension > 0 else 0
                        _, _, _, trimmed_seq = evaluate_cds(frame_seq, start_window, upstream_extension=effective_upstream)
                    else:
                        trimmed_seq = ""
                
                matrix_status[sample][locus] = status
                protein_len = len(protein_seq)
                out.write(
                    f"{sample}\t{locus}\t{status}\t{protein_len}\t{protein_seq}\t{reason}\t{frame_offset}\n"
                )
                if protein_seq:
                    header = f"{sample}|{locus}|{status}|frame{frame_offset}"
                    protein_entries[locus].append((header, protein_seq))
                if trimmed_seq:
                    header = f"{sample}|{locus}|{status}|frame{frame_offset}"
                    trimmed_nucleotide_entries[locus].append((header, trimmed_seq))
    
    write_protein_fastas(protein_dir, protein_entries)
    # Write trimmed nucleotide sequences used for translation
    trimmed_nuc_dir = protein_dir.parent / "per_locus_trimmed_sequences_3frame"
    write_protein_fastas(trimmed_nuc_dir, trimmed_nucleotide_entries)
    write_activity_matrix(matrix_path, samples, loci_sorted, matrix_status)
    render_activity_heatmap(heatmap_path, samples, loci_sorted, matrix_status)
    print(f"3-frame analysis summary written to {summary_path}")
    print(f"3-frame protein FASTAs written to {protein_dir}")
    print(f"3-frame trimmed nucleotide sequences written to {trimmed_nuc_dir}")
    print(f"3-frame activity matrix written to {matrix_path}")
    print(f"3-frame heatmap saved to {heatmap_path}")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    liftoff_out_dir = args.output_dir / "liftoff_gffs"
    liftoff_out_dir.mkdir(exist_ok=True)
    report_path = args.missing_report or args.output_dir / "missing_annotations.tsv"
    blast_report = args.blast_report or args.output_dir / "blast_rescue_hits.tsv"
    status_report = (
        args.annotation_status_report
        or args.output_dir / "final_annotation_status.tsv"
    )
    gene_summary_path = args.gene_summary or args.output_dir / "gene_activity_summary.tsv"
    protein_dir = args.protein_fastas_dir or args.output_dir / "per_locus_proteins"
    matrix_path = args.activity_matrix or args.output_dir / "gene_activity_matrix.tsv"
    heatmap_path = args.activity_heatmap or args.output_dir / "gene_activity_heatmap.png"
    target_fastas = find_target_fastas(args.targets_dir)
    target_fasta_map = {fasta.stem: fasta for fasta in target_fastas}
    reference_sequences = parse_fasta(args.reference_fasta)
    reference_loci = load_reference_loci(args.reference_gff, args.feature_type)
    reference_feature_sequences = extract_reference_feature_sequences(
        args.reference_gff, reference_sequences, args.feature_type
    )
    annotation_status = {
        sample: {locus: "missing" for locus in reference_loci}
        for sample in target_fasta_map
    }
    sample_sequences_cache: dict[str, dict[str, str]] = {}
    genes_by_locus: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)
    missing_by_locus: defaultdict[str, set[str]] = defaultdict(set)
    missing_by_sample: defaultdict[str, set[str]] = defaultdict(set)
    all_extraction_metadata: dict[str, dict[str, dict]] = {}  # Store metadata for debugging

    for target_fasta in target_fastas:
        sample_name = target_fasta.stem
        sample_gff = liftoff_out_dir / f"{sample_name}.gff3"
        if args.overwrite or not sample_gff.exists():
            run_liftoff(
                args.liftoff_bin,
                target_fasta,
                args.reference_fasta,
                args.reference_gff,
                sample_gff,
                args.threads,
                args.liftoff_min_alignment,
                args.liftoff_min_identity,
                args.liftoff_quiet,
                args.liftoff_extra,
            )
        else:
            print(f"[liftoff] Skipping {sample_name}, existing annotation found.")

        sequences = parse_fasta(target_fasta)
        sample_sequences_cache[sample_name] = sequences
        sample_genes, sample_metadata = extract_genes(
            sample_gff, sequences, args.feature_type, sample_name, args.upstream_extension
        )
        
        # Store metadata for debugging
        all_extraction_metadata[sample_name] = sample_metadata
        
        for locus_tag in reference_loci:
            gene_entry = sample_genes.get(locus_tag)
            if gene_entry:
                genes_by_locus[locus_tag].append(gene_entry)
                annotation_status[sample_name][locus_tag] = "liftoff"
            else:
                missing_by_locus[locus_tag].add(sample_name)
                missing_by_sample[sample_name].add(locus_tag)

    print(
        f"Completed Liftoff transfer for {len(genes_by_locus)} locus tags across "
        f"{len(target_fastas)} genomes."
    )
    
    # Write raw extraction metadata and FASTA sequences for debugging
    write_raw_extracted_sequences(all_extraction_metadata, genes_by_locus, args.output_dir)
    
    if args.blast_rescue:
        blast_hits = run_blast_rescue(
            missing_by_sample,
            reference_feature_sequences,
            target_fasta_map,
            args.output_dir,
            args.blastn_bin,
            args.makeblastdb_bin,
            args.blast_min_pident,
            args.blast_min_coverage,
        )
        write_blast_hits(blast_hits, blast_report)
        for hit in blast_hits:
            sample = hit["sample"]
            locus = hit["locus_tag"]
            if annotation_status.get(sample, {}).get(locus) == "liftoff":
                continue
            seq = extract_blast_hit_sequence(hit, sample_sequences_cache)
            if seq:
                header = f"{sample}|{locus}|blast"
                genes_by_locus[locus].append((header, seq))
            annotation_status[sample][locus] = "blast_rescued"
            missing_by_locus[locus].discard(sample)
            missing_by_sample[sample].discard(locus)
        print(
            f"BLAST rescue complete: {len(blast_hits)} candidate hits written to {blast_report}"
        )
    
    write_locus_fastas(genes_by_locus, args.output_dir)
    write_missing_report(missing_by_locus, report_path)
    missing_loci = sum(1 for samples in missing_by_locus.values() if samples)
    if args.blast_rescue:
        print(
            f"Missing gene report written to {report_path} "
            f"({missing_loci} loci absent after Liftoff and BLAST)."
        )
    else:
        print(
            f"Missing gene report written to {report_path} "
            f"({missing_loci} loci absent after Liftoff)."
        )
    write_annotation_status(annotation_status, status_report)
    print(f"Final annotation status table written to {status_report}")
    if not args.skip_gene_status:
        run_gene_activity_analysis(
            args.output_dir / "per_locus_fastas",
            annotation_status,
            list(reference_loci),
            gene_summary_path,
            protein_dir,
            matrix_path,
            heatmap_path,
            start_window=args.start_window,
            upstream_extension=args.upstream_extension,
        )
    
    if args.run_3frame_analysis:
        print("\n" + "="*70)
        print("Running parallel 3-frame translation analysis...")
        print("="*70)
        frame3_summary = args.frame3_summary or args.output_dir / "gene_activity_3frame_summary.tsv"
        frame3_protein_dir = args.frame3_protein_dir or args.output_dir / "per_locus_proteins_3frame"
        frame3_matrix = args.frame3_matrix or args.output_dir / "gene_activity_3frame_matrix.tsv"
        frame3_heatmap = args.frame3_heatmap or args.output_dir / "gene_activity_3frame_heatmap.png"
        
        run_3frame_analysis(
            args.output_dir / "per_locus_fastas",
            list(reference_loci),
            target_fastas,
            frame3_summary,
            frame3_protein_dir,
            frame3_matrix,
            frame3_heatmap,
            start_window=args.start_window,
            upstream_extension=args.upstream_extension,
        )


if __name__ == "__main__":
    main()

