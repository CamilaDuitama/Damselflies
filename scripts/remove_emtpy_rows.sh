#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_file> -o <output_file> -t <threshold> [-j <num_jobs>]"
    echo "  -i: Input file path"
    echo "  -o: Output file path"
    echo "  -t: Threshold percentage (0-100)"
    echo "  -j: Number of parallel jobs (optional, default: number of CPU cores)"
    exit 1
}

# Parse command line arguments
while getopts "i:o:t:j:" opt; do
    case $opt in
        i) input="$OPTARG" ;;
        o) output="$OPTARG" ;;
        t) threshold="$OPTARG" ;;
        j) jobs="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if all required parameters are provided
if [ -z "$input" ] || [ -z "$output" ] || [ -z "$threshold" ]; then
    usage
fi

# Set default number of jobs if not specified
if [ -z "$jobs" ]; then
    jobs=$(nproc)
fi

# Check if input file exists
if [ ! -f "$input" ]; then
    echo "Error: Input file '$input' does not exist."
    exit 1
fi

# Check if output file already exists
if [ -f "$output" ]; then
    echo "Error: Output file '$output' already exists."
    exit 1
fi

# Create a temporary directory for chunks and a progress file
temp_dir=$(mktemp -d)
progress_file="$temp_dir/progress"

# Initialize progress file
echo "0" > "$progress_file"

# Get the number of lines in the input file
total_lines=$(wc -l < "$input")
lines_per_job=$(( (total_lines + jobs - 1) / jobs ))

echo "Total lines: $total_lines"
echo "Number of jobs: $jobs"
echo "Lines per job: $lines_per_job"
echo "Processing..."

# AWK script for processing
awk_script='
BEGIN { FS = OFS = "," }
{
    if (NR == 1) {
        print
        next
    }
    total = NF - 1
    positive = 0
    for (i = 2; i <= NF; i++) {
        if ($i > 0) positive++
    }
    percentage = (positive / total) * 100
    if (percentage >= threshold) {
        print
    }
    if (NR % 1000 == 0) {
        system("echo " NR " >> '"$progress_file"'")
    }
}
END {
    system("echo " NR " >> '"$progress_file"'")
}
'

# Function to update progress
update_progress() {
    local current_progress=$(cat "$progress_file" | sort -n | tail -n 1)
    local progress=$(( (current_progress * 100) / total_lines ))
    echo -ne "Progress: $progress%\r"
}

# Process the file in chunks
for ((i=1; i<=total_lines; i+=lines_per_job)); do
    end=$((i + lines_per_job - 1))
    sed -n "${i},${end}p" "$input" | awk -F, -v threshold="$threshold" "$awk_script" > "$temp_dir/chunk_$i" &
done

# Monitor progress while jobs are running
while [ $(jobs -r | wc -l) -gt 0 ]; do
    update_progress
    sleep 1
done

# Final progress update
update_progress

# Combine results
cat "$temp_dir"/chunk_* > "$output"

# Clean up
rm -rf "$temp_dir"

echo -e "\nProcessing complete. Results written to '$output'."
