import pysam
import sys

# Usage: python bam_ani_editdist.py <bam_file> <reference_fasta>
bam_path = sys.argv[1]
ref_path = sys.argv[2]

bamfile = pysam.AlignmentFile(bam_path, "rb")
ref_fasta = pysam.FastaFile(ref_path)

matches = 0
aligned_bases = 0
total_edit_distance = 0
total_aligned_reads = 0

for read in bamfile.fetch(until_eof=True):
    if read.is_unmapped:
        continue
    
    total_aligned_reads += 1

    # Get reference sequence for this alignment
    ref_seq = ref_fasta.fetch(bamfile.get_reference_name(read.reference_id),
                              read.reference_start, read.reference_end)

    query_seq = read.query_sequence
    aligned_pairs = read.get_aligned_pairs(matches_only=True)

    # Count matches
    for qpos, rpos in aligned_pairs:
        aligned_bases += 1
        if qpos is not None and rpos is not None:
            if query_seq[qpos] == ref_seq[rpos - read.reference_start]:
                matches += 1

    # Calculate edit distance using CIGAR
    # get_cigar_stats returns tuple with counts of operations: (matches, mismatches, insertions, deletions)
    # [1] is the number of mismatches, [2] insertions, [3] deletions
    # This requires pysam >= 0.16
    cigar_stats = read.get_cigar_stats()[0]
    mismatches = cigar_stats[1] if len(cigar_stats) > 1 else 0
    insertions = cigar_stats[2] if len(cigar_stats) > 2 else 0
    deletions = cigar_stats[3] if len(cigar_stats) > 3 else 0

    edit_distance = mismatches + insertions + deletions
    total_edit_distance += edit_distance

bamfile.close()
ref_fasta.close()

ani = (matches / aligned_bases) * 100 if aligned_bases > 0 else 0
mean_edit_distance_per_read = (total_edit_distance / total_aligned_reads) if total_aligned_reads > 0 else 0

print(f"Aligned bases: {aligned_bases}")
print(f"Matches: {matches}")
print(f"ANI: {ani:.4f}%")
print(f"Total aligned reads: {total_aligned_reads}")
print(f"Total edit distance: {total_edit_distance}")
print(f"Mean edit distance per aligned read: {mean_edit_distance_per_read:.2f}")
