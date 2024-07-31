#!/bin/bash

# This script compares all the files named the same within a given path (as argument)
# The dfiference between two files the with same name is printed out
# If the files are identical, nothing is printed
# README files are excluded

search_dir="$1"

files_to_ignore="README.md materials.xml elements.xml CMakeLists.txt"

# Create an associative array to store file paths with the same names
declare -A file_names

# Iterate through files in the directory
find "$search_dir"  -not -path '*/.*' -type f | while read -r file; do
    
    # Get the base name of the file (without the path)
    file_name=$(basename "$file")
    
    # Ignore README files 
    if [[ $files_to_ignore == *"$file_name"*  ]]; then
	continue
    fi

    # Check if there is already a file with the same name
    if [[ -n "${file_names[$file_name]}" ]]; then
        # Compare the contents of the files
	# use flag -q instead of -c to print just the filenames
        diff -c "$file" "${file_names[$file_name]}"
    else
        # Store the file path in the array
        file_names["$file_name"]="$file"
    fi

done

