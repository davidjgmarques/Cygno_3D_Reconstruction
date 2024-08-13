#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 start_run end_run data_folder"
  exit 1
fi

START_RUN=$1
END_RUN=$2
DATA_FOLDER=$3

# Loop through the specified range
for (( run=$START_RUN; run<=$END_RUN; run++ )); do
  INPUT_FILE="${DATA_FOLDER}/reco_run${run}_3D.root"
  
  # Check if the input file exists
  if [ -f "$INPUT_FILE" ]; then
    OUTPUT_FILE="bat_match_only_clusters/out3D_${run}"
    echo "Processing: $INPUT_FILE -> $OUTPUT_FILE"
    # ./alpha_3d.out full "$INPUT_FILE" "$INPUT_FILE" "$OUTPUT_FILE" || true
    ./alpha_3d.out full "$INPUT_FILE" "$INPUT_FILE" "$OUTPUT_FILE"
  else
    echo "File not found: $INPUT_FILE"
  fi
done