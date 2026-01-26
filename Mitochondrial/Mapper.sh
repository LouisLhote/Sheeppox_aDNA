#!/bin/bash
set -euo pipefail

############################################
# FUNCTIONS
############################################

# Convert relative paths to absolute paths
get_abs_path() {
    if [[ -d "$1" ]]; then
        (cd "$1" && pwd)
    elif [[ -f "$1" ]]; then
        local dir=$(dirname "$1")
        local file=$(basename "$1")
        local abs_dir=$(cd "$dir" && pwd)
        echo "$abs_dir/$file"
    else
        # For files/dirs that don't exist yet, resolve the parent directory
        local parent=$(dirname "$1")
        local name=$(basename "$1")
        if [[ -d "$parent" ]]; then
            local abs_parent=$(cd "$parent" && pwd)
            echo "$abs_parent/$name"
        elif [[ "$parent" == "." ]]; then
            local cur_dir=$(pwd)
            echo "$cur_dir/$name"
        else
            echo "$1"  # Return as-is if can't resolve
        fi
    fi
}

############################################
# ARGUMENTS
############################################

while getopts s:r:o: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;     # FASTQ
        r) reference=${OPTARG};;  # reference fasta
        o) outdir=${OPTARG};;     # output directory
    esac
done

# Convert to absolute paths
sample=$(get_abs_path "$sample")
reference=$(get_abs_path "$reference")
outdir=$(get_abs_path "$outdir")

echo "Sample:    $sample"
echo "Reference: $reference"
echo "Output:    $outdir"

############################################
# SAMPLE & REFERENCE NAMES
############################################

name=$(basename "$sample")
SAMPLE=$(echo "$name" | cut -f1 -d'.')

refn=$(basename "$reference")
ref_name=$(echo "$refn" | cut -f1 -d'.')

sample_name=$(echo "$SAMPLE" | cut -f1 -d'-')

read_group="@RG\tID:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}\tSM:${sample_name}"

############################################
# DIRECTORY SETUP
############################################

# Create output directory if it doesn't exist
mkdir -p "$outdir"
cd "$outdir"

mkdir -p sai_file/bam_file/q30/rmdup
mkdir -p final_results/{mapDamage,edit_distance/{tables,plots}}

############################################
# BWA ALIGNMENT (aDNA PARAMETERS)
############################################

echo "Running bwa aln..."

bwa aln \
  -l 1024 \
  -n 0.01 \
  -o 2 \
  -t 8 \
  "$reference" \
  "$sample" \
  > sai_file/${SAMPLE}_${ref_name}.sai

############################################
# SAM → BAM (FILTER UNMAPPED)
############################################

echo "Generating BAM..."

bwa samse \
  -r "$read_group" \
  "$reference" \
  sai_file/${SAMPLE}_${ref_name}.sai \
  "$sample" | \
samtools view -F 4 -Sb - \
> sai_file/bam_file/${SAMPLE}_${ref_name}_F4.bam

############################################
# SORT BAM
############################################

samtools sort \
  sai_file/bam_file/${SAMPLE}_${ref_name}_F4.bam \
  -o sai_file/bam_file/${SAMPLE}_${ref_name}_F4.bam

############################################
# MAPQ ≥ 30
############################################

samtools view -q 30 \
  sai_file/bam_file/${SAMPLE}_${ref_name}_F4.bam \
  -o sai_file/bam_file/q30/${SAMPLE}_${ref_name}_F4_q30.bam

samtools sort \
  sai_file/bam_file/q30/${SAMPLE}_${ref_name}_F4_q30.bam \
  -o sai_file/bam_file/q30/${SAMPLE}_${ref_name}_F4_q30_sort.bam

############################################
# REMOVE DUPLICATES
############################################

java -jar /raid_md0/Software/picard.jar MarkDuplicates \
  OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
  VALIDATION_STRINGENCY=SILENT \
  REMOVE_DUPLICATES=true \
  I=sai_file/bam_file/q30/${SAMPLE}_${ref_name}_F4_q30_sort.bam \
  O=sai_file/bam_file/q30/rmdup/${SAMPLE}_${ref_name}_F4_q30_rmdup.bam \
  M=sai_file/bam_file/q30/rmdup/${SAMPLE}_${ref_name}_F4_q30_rmdup.metrics

############################################
# FINAL BAM
############################################

FINAL_BAM=sai_file/bam_file/q30/rmdup/${SAMPLE}_${ref_name}_F4_q30_rmdup.bam

############################################
# QUALIMAP (FIXED)
############################################

QUALIMAP_DIR=sai_file/bam_file/q30/rmdup/qualimap
mkdir -p ${QUALIMAP_DIR}

qualimap bamqc -bam "$FINAL_BAM"

mv sai_file/bam_file/q30/rmdup/${SAMPLE}_${ref_name}_F4_q30_rmdup_stats \
   ${QUALIMAP_DIR}/
############################################
# MAPDAMAGE
############################################

echo "Running mapDamage..."
samtools index "$FINAL_BAM"
mapDamage \
  -i "$FINAL_BAM" \
  -r "$reference" \
  -d final_results/mapDamage

############################################
# EDIT DISTANCE (NM TAG)
############################################

echo "Extracting edit distances..."

samtools view "$FINAL_BAM" | \
awk '
{
  for (i=12; i<=NF; i++) {
    if ($i ~ /^NM:i:/) {
      split($i,a,":");
      print a[3]
    }
  }
}
' > final_results/edit_distance/tables/NM_values.txt

############################################
# NM COUNTS (0–5)
############################################

awk '$1 <= 5' final_results/edit_distance/tables/NM_values.txt | \
sort -n | uniq -c | \
awk '{print $2"\t"$1}' \
> final_results/edit_distance/tables/edit_distance_0_5.tsv

############################################
# PAIRWISE DISTANCE TO REFERENCE
############################################

samtools view "$FINAL_BAM" | \
awk '
{
  nm="NA"
  for (i=12; i<=NF; i++) {
    if ($i ~ /^NM:i:/) {
      split($i,a,":");
      nm=a[3]
    }
  }
  if (nm != "NA") {
    len=length($10)
    if (len > 0)
      print nm/len
  }
}
' > final_results/edit_distance/tables/pairwise_distance.txt

############################################
# PLOTTING (PYTHON)
############################################

echo "Plotting results.."

############################################
# PLOTTING WITH R (BASE R ONLY)
############################################

# Create R script file
R_SCRIPT="final_results/edit_distance/plots/plot_script.R"

cat > "$R_SCRIPT" << 'EOF'

# ---------- Edit distance plot (0–5) ----------
ed <- read.table(
  "final_results/edit_distance/tables/edit_distance_0_5.tsv",
  header = FALSE
)

png(
  "final_results/edit_distance/plots/edit_distance_0_5.png",
  width = 1000,
  height = 800
)

barplot(
  height = ed$V2,
  names.arg = ed$V1,
  xlab = "Number of mismatches (NM)",
  ylab = "Number of reads",
  main = "Edit distance distribution (0-5 mismatches)",
  col = "steelblue"
)

dev.off()

# ---------- Pairwise distance plot ----------
dist <- scan(
  "final_results/edit_distance/tables/pairwise_distance.txt",
  quiet = TRUE
)

png(
  "final_results/edit_distance/plots/pairwise_distance.png",
  width = 1000,
  height = 800
)

hist(
  dist,
  breaks = 50,
  xlab = "Pairwise distance (NM / aligned length)",
  ylab = "Number of reads",
  main = "Pairwise distance to reference",
  col = "darkseagreen"
)

dev.off()

EOF

Rscript "$R_SCRIPT"

# Clean up R script file
rm -f "$R_SCRIPT"

############################################
# SUMMARY STATISTICS
############################################

echo "Generating summary statistics..."

SUMMARY_FILE="final_results/summary_statistics.tsv"

# 1. Number of reads mapped
NUM_READS=$(samtools view -c "$FINAL_BAM")

# 2. Damage at first base for 5' to T (from mapDamage)
# mapDamage creates 5pCtoT_freq.txt with position and 5pC>T frequency
# Format: pos     5pC>T
# We want position 1 (first base, 5' end)
DAMAGE_5T="NA"

# First, try to parse from 5pCtoT_freq.txt (most direct method)
if [[ -f "final_results/mapDamage/5pCtoT_freq.txt" ]]; then
    # Extract value at position 1 (first data row after header)
    DAMAGE_5T=$(awk 'NR==2 && NF>=2 {print $2; exit}' final_results/mapDamage/5pCtoT_freq.txt 2>/dev/null)
    # Validate it's a number
    if [[ -n "$DAMAGE_5T" ]]; then
        if ! [[ "$DAMAGE_5T" =~ ^[0-9]*\.?[0-9]+$ ]]; then
            DAMAGE_5T="NA"
        fi
    fi
fi

# Fallback: try to parse from misincorporation.txt
if [[ -z "$DAMAGE_5T" || "$DAMAGE_5T" == "NA" ]]; then
    if [[ -f "final_results/mapDamage/misincorporation.txt" ]]; then
    # First, check the header to find the C->T column
    C_TO_T_COL=$(head -1 final_results/mapDamage/misincorporation.txt | tr '\t' '\n' | grep -n -i "C.*T\|C->T" | head -1 | cut -d: -f1)
    
    if [[ -n "$C_TO_T_COL" ]]; then
        # Extract C->T at position 1 (first data row after header)
        DAMAGE_5T=$(awk -v col="$C_TO_T_COL" 'NR==2 && NF>=col {print $col; exit}' final_results/mapDamage/misincorporation.txt 2>/dev/null)
    else
        # Fallback: assume C->T is column 7 (common format)
        DAMAGE_5T=$(awk 'NR==2 && NF>=7 {print $7; exit}' final_results/mapDamage/misincorporation.txt 2>/dev/null)
    fi
    
    # Validate it's a number (0-1 range for damage) and not a path
    if [[ -n "$DAMAGE_5T" ]]; then
        if ! [[ "$DAMAGE_5T" =~ ^[0-9]*\.?[0-9]+$ ]] || [[ "$DAMAGE_5T" =~ "/" ]] || [[ "$DAMAGE_5T" =~ "home" ]]; then
            DAMAGE_5T="NA"
        fi
    fi
    fi
fi

# Alternative: parse using R (more robust)
if [[ -z "$DAMAGE_5T" || "$DAMAGE_5T" == "NA" ]]; then
    DAMAGE_5T=$(Rscript -e "
        options(warn=-1)
        damage_val <- 'NA'
        # First try 5pCtoT_freq.txt (most direct)
        if (file.exists('final_results/mapDamage/5pCtoT_freq.txt')) {
            freq <- tryCatch({
                read.table('final_results/mapDamage/5pCtoT_freq.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE, comment.char='')
            }, error=function(e) NULL)
            if (!is.null(freq) && nrow(freq) > 0) {
                # Get value at position 1
                pos1_row <- freq[freq[,1] == 1, ]
                if (nrow(pos1_row) > 0 && ncol(pos1_row) >= 2) {
                    val <- pos1_row[1, 2]
                    if (!is.na(val) && is.numeric(val)) {
                        damage_val <- as.character(val)
                    }
                }
            }
        }
        # Fallback to misincorporation.txt
        if (damage_val == 'NA' && file.exists('final_results/mapDamage/misincorporation.txt')) {
            mis <- tryCatch({
                read.table('final_results/mapDamage/misincorporation.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE, comment.char='', fill=TRUE)
            }, error=function(e) NULL)
            if (!is.null(mis) && nrow(mis) > 0) {
                # Find C->T column by name
                col_names <- colnames(mis)
                c_to_t_col <- grep('C.*T|C->T', col_names, ignore.case=TRUE, value=FALSE)[1]
                if (length(c_to_t_col) > 0 && !is.na(c_to_t_col)) {
                    val <- mis[1, c_to_t_col]
                } else if (ncol(mis) >= 7) {
                    # Fallback to column 7
                    val <- mis[1, 7]
                } else {
                    val <- NA
                }
                if (!is.na(val) && is.numeric(val)) {
                    damage_val <- as.character(val)
                }
            }
        }
        cat(damage_val)
    " 2>/dev/null)
    
    # Final validation - ensure it's a number and not a path
    if [[ -n "$DAMAGE_5T" ]]; then
        if ! [[ "$DAMAGE_5T" =~ ^[0-9]*\.?[0-9]+$ ]] || [[ "$DAMAGE_5T" =~ "/" ]] || [[ "$DAMAGE_5T" =~ "home" ]]; then
            DAMAGE_5T="NA"
        fi
    fi
fi

# Final fallback: set to NA if still invalid or contains path characters
if [[ -z "$DAMAGE_5T" ]] || [[ "$DAMAGE_5T" =~ "/" ]] || [[ "$DAMAGE_5T" =~ "home" ]] || ! [[ "$DAMAGE_5T" =~ ^[0-9]*\.?[0-9]+$|^NA$ ]]; then
    DAMAGE_5T="NA"
fi

# 3. Evenness of coverage (coefficient of variation of coverage)
# Calculate coverage using samtools depth
COVERAGE_FILE="final_results/coverage_depth.txt"
# Ensure BAM is indexed
samtools index "$FINAL_BAM" 2>/dev/null
# Calculate depth - try with reference chromosome name first, then fallback
REF_CHR=$(samtools view -H "$FINAL_BAM" | grep '^@SQ' | head -1 | cut -f2 | sed 's/SN://')
if [[ -n "$REF_CHR" ]]; then
    samtools depth -a -r "$REF_CHR" "$FINAL_BAM" > "$COVERAGE_FILE" 2>/dev/null || \
    samtools depth -a "$FINAL_BAM" > "$COVERAGE_FILE" 2>/dev/null || touch "$COVERAGE_FILE"
else
    samtools depth -a "$FINAL_BAM" > "$COVERAGE_FILE" 2>/dev/null || touch "$COVERAGE_FILE"
fi

EVENNESS="NA"
if [[ -s "$COVERAGE_FILE" ]]; then
    EVENNESS=$(awk '
    {
        if (NF >= 3 && $3 > 0) {
            sum += $3
            sum_sq += $3 * $3
            n++
        }
    }
    END {
        if (n > 0) {
            mean = sum / n
            if (mean > 0) {
                variance = (sum_sq / n) - (mean * mean)
                stddev = sqrt(variance)
                cv = stddev / mean
                # Evenness is 1 / (1 + CV) - higher is more even
                evenness = 1 / (1 + cv)
                printf "%.6f", evenness
            } else {
                print "NA"
            }
        } else {
            print "NA"
        }
    }' "$COVERAGE_FILE")
fi

# 4. Depth of coverage (mean depth)
# Method 1: Calculate as total mapped bases / reference length
DEPTH="NA"

# Get reference length from BAM header (try multiple methods)
REF_LENGTH=$(samtools view -H "$FINAL_BAM" 2>/dev/null | grep '^@SQ' | head -1 | awk -F'\t' '{for(i=1;i<=NF;i++){if($i~/^LN:/){sub(/LN:/,"",$i); print $i; exit}}}')

# Alternative: get from reference file if available
if [[ -z "$REF_LENGTH" || "$REF_LENGTH" -le 0 ]]; then
    if [[ -f "$reference" ]]; then
        REF_LENGTH=$(grep -v '^>' "$reference" 2>/dev/null | tr -d '\n' | wc -c)
    fi
fi

# Calculate total mapped bases (sum of read sequence lengths)
TOTAL_BASES=$(samtools view "$FINAL_BAM" 2>/dev/null | awk '{sum += length($10)} END {print sum+0}')

# Debug: ensure we have values (uncomment for debugging)
# echo "DEBUG: REF_LENGTH=$REF_LENGTH, TOTAL_BASES=$TOTAL_BASES" >&2

# Calculate depth if we have both values
if [[ -n "$REF_LENGTH" && "$REF_LENGTH" -gt 0 && -n "$TOTAL_BASES" && "$TOTAL_BASES" -gt 0 ]]; then
    # Use awk for calculation: Depth = Total Mapped Bases / Reference Length
    # Use 6 decimal places to capture very small depth values
    DEPTH=$(awk -v bases="$TOTAL_BASES" -v len="$REF_LENGTH" 'BEGIN {
        if (len > 0 && bases > 0) {
            depth = bases / len
            printf "%.6f", depth
        } else {
            print "NA"
        }
    }')
    
    # Validate the result is a number
    if [[ ! "$DEPTH" =~ ^[0-9]+\.?[0-9]*$ ]]; then
        DEPTH="NA"
    fi
fi

# Method 2: Fallback - use samtools depth output
if [[ -z "$DEPTH" || "$DEPTH" == "NA" || "$DEPTH" == "0.00" ]]; then
    if [[ -s "$COVERAGE_FILE" ]]; then
        # Calculate mean depth from samtools depth output (mean across all positions)
        DEPTH=$(awk 'NF >= 3 && $3 != "" {sum += $3; n++} END {if (n > 0) {avg = sum/n; printf "%.2f", avg} else print "NA"}' "$COVERAGE_FILE")
    fi
fi

# Method 3: Use samtools stats as final fallback
if [[ -z "$DEPTH" || "$DEPTH" == "NA" || "$DEPTH" == "0.00" ]]; then
    STATS_FILE="final_results/temp_stats.txt"
    samtools stats "$FINAL_BAM" > "$STATS_FILE" 2>/dev/null
    if [[ -f "$STATS_FILE" ]]; then
        # Extract total bases from stats
        STATS_BASES=$(grep "^SN" "$STATS_FILE" | grep "bases mapped" | awk '{print $4}')
        if [[ -n "$STATS_BASES" && "$STATS_BASES" -gt 0 && -n "$REF_LENGTH" && "$REF_LENGTH" -gt 0 ]]; then
            DEPTH=$(awk -v bases="$STATS_BASES" -v len="$REF_LENGTH" 'BEGIN {if (len > 0 && bases > 0) printf "%.2f", bases/len; else print "NA"}')
        fi
        rm -f "$STATS_FILE"
    fi
fi

# 5. 1 - mean pairwise distance
MEAN_PAIRWISE_DIST="NA"
if [[ -f "final_results/edit_distance/tables/pairwise_distance.txt" && -s "final_results/edit_distance/tables/pairwise_distance.txt" ]]; then
    MEAN_PAIRWISE_DIST=$(awk '{sum += $1; n++} END {if (n > 0) {mean = sum/n; result = 1 - mean; printf "%.6f", result} else print "NA"}' final_results/edit_distance/tables/pairwise_distance.txt)
fi

# Final safety check: ensure DAMAGE_5T is valid before writing
if [[ -z "$DAMAGE_5T" ]] || [[ "$DAMAGE_5T" =~ "/" ]] || [[ "$DAMAGE_5T" =~ "home" ]] || [[ "$DAMAGE_5T" =~ "TUTORIAL" ]] || ! [[ "$DAMAGE_5T" =~ ^[0-9]*\.?[0-9]+$|^NA$ ]]; then
    DAMAGE_5T="NA"
fi

# Write summary TSV file
{
    echo -e "Sample\tReference\tNum_reads_mapped\tDamage_5prime_to_T\tEvenness_of_coverage\tDepth_of_coverage\tOne_minus_mean_pairwise_distance"
    echo -e "${SAMPLE}\t${ref_name}\t${NUM_READS}\t${DAMAGE_5T}\t${EVENNESS}\t${DEPTH}\t${MEAN_PAIRWISE_DIST}"
} > "$SUMMARY_FILE"

# Clean up temporary coverage file
rm -f "$COVERAGE_FILE"

echo "Summary statistics written to: $SUMMARY_FILE"
echo "Pipeline finished successfully."

