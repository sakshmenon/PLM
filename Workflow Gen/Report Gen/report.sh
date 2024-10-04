#! /usr/bin/bash

module load cuda/12.2
module load anaconda3
source activate bioembeddings

# Check if the user provided at least three arguments
if [ $# -lt 5 ]; then
    echo "Usage: $0 <reference_file> <query_file> <format> <top_matches> <method> [output_file]"
    exit 1
fi

reference_file=$1
query_file=$2
format=$3
top_matches=$4
method=$5

# Optional output file (if provided)
output_file=$6

# Run the first command and generate reference_file.csv
# python cluster_run.py -ff "${reference_file}.faa" -e prott5 -of "${reference_file}" -l protein

# Run the second command and generate query_file.csv
# python cluster_run.py -ff "${query_file}.faa" -e prott5 -of "${query_file}" -l protein

# Run the third command, using the .csv outputs from the previous commands as input
if [ -z "$output_file" ]; then
    # If output_file is not provided, run without saving to a specific file
    python search_similar_proteins.py -rf "${reference_file}.csv" -qf "${query_file}.csv" -m "${method}" -k "${top_matches}" -f "${format}"
else
    # If output_file is provided, run and save results to the specified location
    python search_similar_proteins.py -rf "${reference_file}.csv" -qf "${query_file}.csv" -m "${method}" -k "${top_matches}" -f "${format}" -of "${output_file}"
    echo "Search complete. Results saved to ${output_file}"
fi