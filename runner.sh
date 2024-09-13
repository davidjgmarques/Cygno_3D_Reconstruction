#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 start_run end_run data_folder output_folder"
  exit 1
fi

START_RUN=$1
END_RUN=$2
DATA_FOLDER=$3
OUTPUT_FOLDER=$4

# Create the output folder if it doesn't exist
if [ ! -d "$OUTPUT_FOLDER" ]; then
  mkdir -p "$OUTPUT_FOLDER"
fi

# Function to process a single file
process_file() {
  local input_file=$1
  local output_file=$2
  echo "Processing: $input_file -> $output_file"
  ./alpha_3d.out full "$input_file" "$output_file" >> "log/log_$(basename "$output_file").txt" 2>&1
}

# Check if start_run and end_run are both -1
if [ "$START_RUN" -eq -1 ] && [ "$END_RUN" -eq -1 ]; then
  # Process all files in the data folder
  for input_file in "$DATA_FOLDER"/*.root; do
    if [ -f "$input_file" ]; then
      output_file="${OUTPUT_FOLDER}/out_3D_$(basename "$input_file" .root)"
      process_file "$input_file" "$output_file"
    else
      echo "No .root files found in $DATA_FOLDER"
    fi
  done
else
  # Loop through the specified range
  for (( run=$START_RUN; run<=$END_RUN; run++ )); do
    input_file="${DATA_FOLDER}/reco_run${run}_3D.root"
    
    # Check if the input file exists
    if [ -f "$input_file" ]; then
      output_file="${OUTPUT_FOLDER}/out_3D_$(basename "$input_file" .root)"
      process_file "$input_file" "$output_file"
    else
      echo "File not found: $input_file"
    fi
  done
fi