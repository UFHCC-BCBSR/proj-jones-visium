#!/bin/bash

# Output file
output_file="../data/combined_metrics_summary.csv"

# Flag to track if header is written
header_written=0

# Loop through all sample directories
for metrics_file in /orange/cancercenter-dept/JONES/AM_10x_Visium_v2/SpaceRanger_Output_and_Scripts/*/outs/metrics_summary.csv; do
  if [ -f "$metrics_file" ]; then
    if [ $header_written -eq 0 ]; then
      # Write header from first file
      head -n 1 "$metrics_file" > "$output_file"
      header_written=1
    fi
    # Append data (skip header)
    tail -n +2 "$metrics_file" >> "$output_file"
  fi
done

echo "Combined metrics saved to: $output_file"
echo "Total samples found: $(tail -n +2 "$output_file" | wc -l)"

# Display the file
cat "$output_file"
