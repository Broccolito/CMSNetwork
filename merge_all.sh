#!/bin/bash

# Output file
output="v2g.csv"

# List of all CSVs to merge
files=( csv_files/*.csv )
total=${#files[@]}

# Write header from the first file
head -n 1 "${files[0]}" > "$output"

# Progress variables
count=0

# Loop through each file
for file in "${files[@]}"; do
    ((count++))
    # Skip header for all but the first file
    tail -n +2 "$file" >> "$output"

    # Display progress bar
    progress=$((count * 100 / total))
    bar=$(printf "%-${progress}s" "#" | tr ' ' '#')
    echo -ne "Merging: [${bar:0:50}] $progress% ($count/$total)\r"
done

echo -e "\nDone: merged into $output"

