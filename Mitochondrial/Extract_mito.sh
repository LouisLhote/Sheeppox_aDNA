#!/bin/bash

# Reference contigs per host species
declare -A REF
REF["Goat"]="gi|316926505|ref|NC_005044.2|"
REF["Sheep"]="NC_001941.1"
REF["Cattle"]="V00654.1"

# Sample → species assignments
declare -A SPECIES
SPECIES["16R"]="Goat Cattle"
SPECIES["48R"]="Goat"
SPECIES["79RR"]="Goat Sheep"
SPECIES["80R"]="Goat"
SPECIES["ANC255"]="Sheep"
SPECIES["CL299"]="Cattle Sheep"
SPECIES["MYRZ180"]="Sheep"
SPECIES["PEM01"]="Cattle Sheep"
SPECIES["PEM02"]="Cattle"
SPECIES["PEM03"]="Cattle Sheep"
SPECIES["PEM04"]="Cattle"
SPECIES["PEM05"]="Cattle Sheep"
SPECIES["RF002"]="Sheep"
SPECIES["RF003"]="Sheep"
SPECIES["RF009"]="Sheep"
SPECIES["RUB002"]="Sheep"
SPECIES["TCC02"]="Cattle"
SPECIES["TP01"]="Cattle Sheep"
SPECIES["TP02"]="Cattle Sheep"
SPECIES["TP04"]="Sheep"
SPECIES["YG01"]="Sheep Goat"
SPECIES["YG15"]="Cattle"
SPECIES["YG16"]="Cattle"
SPECIES["BA27"]="Cattle Sheep"

# I/O directories
BAMDIR="."
OUTDIR="./filtered_species_bams"
mkdir -p "$OUTDIR"

for BAM in *_Ruminant_mito_merged_merged.bam; do

    SAMPLE=$(echo "$BAM" | awk -F"_" '{print $1}')
    HOST_LIST=${SPECIES[$SAMPLE]}

    if [ -z "$HOST_LIST" ]; then
        echo "⚠️  No species mapping for sample $SAMPLE — skipping."
        continue
    fi

    echo "Processing $SAMPLE → species: $HOST_LIST"

    for HOST in $HOST_LIST; do
        CONTIG=${REF[$HOST]}

        OUT="$OUTDIR/${SAMPLE}_${HOST}.bam"

        echo "  Extracting $HOST → contig $CONTIG"
        samtools view -b "$BAM" "$CONTIG" > "$OUT"
        samtools index "$OUT"
    done
done
