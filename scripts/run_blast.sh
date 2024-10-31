#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 -r <reference> -q <query> -o <output_prefix> [-t <threads>]"
    echo "  -r: Reference genome (can be gzipped)"
    echo "  -q: Query file (can be gzipped)"
    echo "  -o: Output prefix"
    echo "  -t: Number of threads (optional, default: 40)"
    exit 1
}

# Parse command line arguments
while getopts "r:q:o:t:" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        q) QUERY="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if required arguments are provided
if [ -z "$REFERENCE" ] || [ -z "$QUERY" ] || [ -z "$OUTPUT_PREFIX" ]; then
    usage
fi

# Set default threads if not provided
THREADS=${THREADS:-40}

# Load BLAST module (adjust as necessary for your system)
module load blast

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_PREFIX")"

# Uncompress files if they're gzipped
if [[ "$REFERENCE" == *.gz ]]; then
    REFERENCE_UNZIPPED="${REFERENCE%.gz}"
    gunzip -c "$REFERENCE" > "$REFERENCE_UNZIPPED"
    REFERENCE="$REFERENCE_UNZIPPED"
fi

if [[ "$QUERY" == *.gz ]]; then
    QUERY_UNZIPPED="${QUERY%.gz}"
    gunzip -c "$QUERY" > "$QUERY_UNZIPPED"
    QUERY="$QUERY_UNZIPPED"
fi

# Create BLAST database from the reference
makeblastdb -in "$REFERENCE" -dbtype nucl

# Run BLAST
blastn -task blastn-short \
    -query "$QUERY" \
    -db "$REFERENCE" \
    -evalue 0.01 \
    -out "${OUTPUT_PREFIX}.tsv" \
    -outfmt 6 \
    -max_target_seqs 5 \
    -num_threads "$THREADS"

# Filter results
awk -F"\t" '{if ($3 >99 && $4 >30 ) print $0,$2,$9}' "${OUTPUT_PREFIX}.tsv" > "${OUTPUT_PREFIX}_filtered.tsv"

# Clean up
if [[ "$REFERENCE" == "$REFERENCE_UNZIPPED" ]]; then
    rm "$REFERENCE_UNZIPPED"
fi
if [[ "$QUERY" == "$QUERY_UNZIPPED" ]]; then
    rm "$QUERY_UNZIPPED"
fi

echo "BLAST mapping completed. Results stored in ${OUTPUT_PREFIX}.tsv"
echo "Filtered results stored in ${OUTPUT_PREFIX}_filtered.tsv"