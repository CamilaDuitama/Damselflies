#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status.

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_matrix_file>"
    exit 1
fi

input_file="$1"
num_cores=$(nproc)  # Get the number of CPU cores

# Define output files for A and f
output_a="A.mat"
output_f="f.mat"

# Create a temporary directory for chunks
temp_dir=$(mktemp -d)
echo "Temporary directory: $temp_dir"

# Function to process a chunk of the file
process_chunk() {
    local chunk_filename="$1"
    local chunk_number="$2"
    local a_output="${temp_dir}/A.mat.${chunk_number}.tmp"
    local f_output="${temp_dir}/f.mat.${chunk_number}.tmp"

    awk -F " " '
    {
        a_line = "";
        f_line = "";
        for (i = 1; i <= NF; i++) {
            if (i == 1) {
                a_line = a_line $i;
                f_line = f_line $i;
            } else {
                split($i, parts, ";");
                a_line = a_line " " (parts[1] == "" ? "NA" : parts[1]);
                f_line = f_line " " (parts[2] == "" ? "NA" : parts[2]);
            }
        }
        print a_line > "'"$a_output"'";
        print f_line > "'"$f_output"'";
    }' "${chunk_filename}"

    echo "Chunk ${chunk_filename}: A=$(wc -l < "$a_output"), F=$(wc -l < "$f_output")"
}

export -f process_chunk
export temp_dir

# Get the number of lines in the input file
total_lines=$(wc -l < "${input_file}")
lines_per_chunk=$((total_lines / num_cores + 1))

echo "Total lines: $total_lines"
echo "Number of cores: $num_cores"
echo "Lines per chunk: $lines_per_chunk"

# Split the file for parallel processing
split -l "$lines_per_chunk" --numeric-suffixes=1 --additional-suffix=".split" "${input_file}" "${temp_dir}/input_chunk_"

# Process chunks and show progress bar using 'parallel'
ls "${temp_dir}"/input_chunk_*.split | parallel --bar -j "${num_cores}" process_chunk {} {#}

# Combine the results from each chunk
> "${output_a}"  # Clear the file if it exists
> "${output_f}"  # Clear the file if it exists

for tmp_file in "${temp_dir}"/A.mat.*.tmp; do
    cat "${tmp_file}" >> "${output_a}"
done

for tmp_file in "${temp_dir}"/f.mat.*.tmp; do
    cat "${tmp_file}" >> "${output_f}"
done

# Clean up
rm -rf "${temp_dir}"

echo "Files have been successfully split into ${output_a} and ${output_f}"

# Verify row counts
input_rows=$(wc -l < "${input_file}")
a_rows=$(wc -l < "${output_a}")
f_rows=$(wc -l < "${output_f}")

echo "Input file rows: $input_rows"
echo "A.mat rows: $a_rows"
echo "f.mat rows: $f_rows"

if [ "$input_rows" -eq "$a_rows" ] && [ "$input_rows" -eq "$f_rows" ]; then
    echo "Row counts match."
else
    echo "Warning: Row counts do not match."
fi

# Additional verification
echo "First few lines of A.mat:"
head -n 5 "${output_a}"
echo "First few lines of f.mat:"
head -n 5 "${output_f}"
echo "Last few lines of A.mat:"
tail -n 5 "${output_a}"
echo "Last few lines of f.mat:"
tail -n 5 "${output_f}"
