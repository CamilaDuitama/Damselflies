#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status.

# Error handling function
error_handler() {
    local line_number=$1
    local command=$2
    local error_code=$3
    echo "Error on line $line_number: Command '$command' exited with status $error_code"
    if [ "$command" != "((current_task++))" ]; then
        exit $error_code
    fi
}

# Set up error trapping
trap 'error_handler ${LINENO} "$BASH_COMMAND" $?' ERR

# Function to display a progress bar
show_progress() {
    local current=$1
    local total=$2
    local width=50
    local percentage=$((current * 100 / total))
    local completed=$((width * current / total))
    local remaining=$((width - completed))
    printf "\rProgress: [%s%s] %d%%" "$(printf '#%.0s' $(seq 1 $completed))" "$(printf ' %.0s' $(seq 1 $remaining))" "$percentage"
}

# Parse command line arguments
while getopts "r:" opt; do
  case $opt in
    r) REF_GENOME="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# Check if reference genome is provided
if [ -z "$REF_GENOME" ]; then
    echo "Error: Reference genome not provided. Use -r flag to specify the reference genome."
    exit 1
fi

# Convert to absolute path
REF_GENOME=$(realpath "$REF_GENOME")

# Verify if the reference genome file exists
if [ ! -f "$REF_GENOME" ]; then
    echo "Error: The specified reference genome file does not exist: $REF_GENOME"
    exit 1
fi

# Get the number of available CPU cores
THREADS=$(nproc)

# Store the base directory
BASE_DIR=$(pwd)

# Function to convert BAM to BLAST-like format
convert_bam_to_blast() {
    local bam_file=$1
    local output_file=$2

    samtools view "$bam_file" | awk -v OFS='\t' '
    {
        qname=$1;  # Query name (kmer ID)
        flag=$2;
        rname=$3;  # Reference sequence name
        pos=$4;    # 1-based leftmost mapping position
        mapq=$5;   # Mapping quality
        cigar=$6;
        seq=$10;   # Sequence

        # Calculate alignment length from CIGAR string
        aln_len=0;
        while (match(cigar, /[0-9]+[MIDNSHP]/)) {
            num = substr(cigar, RSTART, RLENGTH-1) + 0;
            op = substr(cigar, RSTART+RLENGTH-1, 1);
            if (op ~ /[MID]/) aln_len += num;
            cigar = substr(cigar, RSTART+RLENGTH);
        }

        # Calculate end position
        end = pos + aln_len - 1;

        # Placeholder values
        pident = 100;  # Percent identity
        mismatch = 0;
        gapopen = 0;
        qstart = 1;
        qend = length(seq);
        evalue = "1e-10";
        bitscore = 60;

        print qname, rname, pident, aln_len, mismatch, gapopen, qstart, qend, pos, end, evalue, bitscore, rname, pos;
    }' > "$output_file"
}

# Function to process each file
process_file() {
    local folder=$1
    local comparison=$2
    
    echo "Processing $comparison in folder: $folder"
    
    cd "$folder"

    # Create output directory for BAM and BAI files
    local output_dir="mapping_output"
    mkdir -p "$output_dir"

    local fasta="significant_unitigs_${comparison}.fasta.gz"
    
    if [ ! -f "$fasta" ]; then
        echo "Error: $fasta not found in $folder. Skipping."
        cd "$BASE_DIR"
        return 1
    fi

    local base=$(basename "$fasta" .fasta.gz)
    
    echo "Processing $base..."

    # Mapping with minimap2 using all available threads
    echo "Running minimap2..."
    minimap2 -t $THREADS -ax sr "$REF_GENOME" "$fasta" > "$output_dir/${base}_mapped.sam"
    
    # Convert SAM to BAM, sort, and index using all available threads
    echo "Converting SAM to BAM and sorting..."
    samtools view -@ $THREADS -bS "$output_dir/${base}_mapped.sam" | samtools sort -@ $THREADS -o "$output_dir/${base}_mapped_sorted.bam"

    echo "Indexing BAM file..."
    samtools index -@ $THREADS "$output_dir/${base}_mapped_sorted.bam"

    # Convert BAM to BLAST-like format
    echo "Converting BAM to BLAST-like format..."
    convert_bam_to_blast "$output_dir/${base}_mapped_sorted.bam" "$output_dir/${base}_blast_like.tsv"

    # Count total reads
    echo "Counting total reads..."
    local total_reads=$(zcat "$fasta" | grep -c '^>')

    # Find the correct contig name
    echo "Finding contig name..."
    local contig_name=$(samtools view -H "$output_dir/${base}_mapped_sorted.bam" | grep -oP 'SN:SUPER_13_unloc_2\S*' | cut -d':' -f2)

    if [ -z "$contig_name" ]; then
        echo "Warning: Could not find SUPER_13_unloc_2 contig in ${base}_mapped_sorted.bam"
        cd "$BASE_DIR"
        return 1
    fi

    # Count reads mapping to the found contig
    echo "Counting mapped reads..."
    local mapped_reads=$(samtools view -@ $THREADS "$output_dir/${base}_mapped_sorted.bam" "$contig_name" | wc -l)

    # Calculate percentage
    local percentage=$(echo "scale=2; $mapped_reads / $total_reads * 100" | bc)

    echo "$base: $percentage% of reads map to $contig_name" >> "$output_dir/mapping_results.txt"

    # Clean up intermediate files
    rm "$output_dir/${base}_mapped.sam"

    cd "$BASE_DIR"
}

# Main script
folders=(
    "minimap/filter_0.1_bin_thr_0.8_chi2"
    "minimap/filter_0.1_bin_thr_0_chi2"
    "minimap/filter_0.2_bin_thr_0.8_chi2"
    "minimap/filter_0.2_bin_thr_0_chi2"
)

comparisons=("AvsO" "AvsI" "OvsI")

total_tasks=$((${#folders[@]} * ${#comparisons[@]}))
current_task=0

echo "Found ${#folders[@]} folders to process with ${#comparisons[@]} comparisons each."

for folder in "${folders[@]}"; do
    echo "Processing folder: $folder"
    for comparison in "${comparisons[@]}"; do
        process_file "$folder" "$comparison"
        current_task=$((current_task + 1))
        show_progress $current_task $total_tasks
    done
done

echo "Processing complete. Results are stored in mapping_output/mapping_results.txt in each processed folder."