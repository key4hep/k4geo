#!/bin/bash

# This script compares all the files named the same within a given path (as argument)
# If files with the same name have different contents, it prints the paths of those files
# If the files are identical, nothing is printed
# Files listed in an ignore file are excluded

# Usage: ./script.sh <search_dir> <ignore_file>

search_dir="$1"
ignore_file="$2"

# Read filenames to ignore from the ignore file
files_to_ignore=$(<"$ignore_file")

# Create an associative array to store file paths with the same names
declare -A file_names

# Initialize status code variable
status_code=0

# Iterate through files in the directory
while IFS= read -r -d '' file; do
    # Get the base name of the file (without the path)
    file_name=$(basename "$file")

    # Ignore files listed in the ignore file
    if grep -qx "$file_name" <<< "$files_to_ignore"; then
        continue
    fi

    # Check if there is already a file with the same name
    if [[ -n "${file_names[$file_name]}" ]]; then
        # Compare the contents of the files
        if ! cmp -s "$file" "${file_names[$file_name]}"; then
            # Print the paths of the files if they differ
            echo "Error. Files with the same name but different contents:"
            echo "$file"
            echo "${file_names[$file_name]}"
            echo
            status_code=1
        fi
    else
        # Store the file path in the array
        file_names["$file_name"]="$file"
    fi
done < <(find "$search_dir" -not -path '*/.*' -type f -print0)

exit $status_code
