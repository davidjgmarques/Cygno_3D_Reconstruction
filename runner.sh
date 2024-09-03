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

# Loop through the specified range
for (( run=$START_RUN; run<=$END_RUN; run++ )); do
  INPUT_FILE="${DATA_FOLDER}/reco_run${run}_3D.root"
  
  # Check if the input file exists
  if [ -f "$INPUT_FILE" ]; then
    OUTPUT_FILE="${OUTPUT_FOLDER}/out_3D_${run}"
    echo "Processing: $INPUT_FILE -> $OUTPUT_FILE"
    ./alpha_3d.out full "$INPUT_FILE" "$OUTPUT_FILE" >> "log/log_${run}.txt" 2>&1
  else
    echo "File not found: $INPUT_FILE"
  fi
done