#!/bin/bash

# This script compares all the files named the same within a given path (as argument)
# If files with the same name have different contents, it prints the paths of those files
# If the files are identical, nothing is printed
# README files are excluded

search_dir="$1"

files_to_ignore="README.md materials.xml elements.xml CMakeLists.txt"

# Create an associative array to store file paths with the same names
declare -A file_names

# Iterate through files in the directory
find "$search_dir" -not -path '*/.*' -type f | while read -r file; do

    # Get the base name of the file (without the path)
    file_name=$(basename "$file")

    # Ignore README files
    if [[ $files_to_ignore == *"$file_name"* ]]; then
        continue
    fi

    # Check if there is already a file with the same name
    if [[ -n "${file_names[$file_name]}" ]]; then
        # Compare the contents of the files
        if ! cmp -s "$file" "${file_names[$file_name]}"; then
            # Print the paths of the files if they differ
            echo "Error: files with the same name but different contents:"
            echo "$file"
            echo "${file_names[$file_name]}"
            echo
        fi
    else
        # Store the file path in the array
        file_names["$file_name"]="$file"
    fi

done

