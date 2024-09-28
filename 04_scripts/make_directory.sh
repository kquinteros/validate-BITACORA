#!/bin/bash
# Author: Kevin Quinteros
# Date: December 13, 2023
# Purpose: This script creates nested directories for InterProScan and validation steps based on user-provided input.
#          It takes an input directory and number ranges for 's' and 'p' as arguments and creates the necessary directory structure.
#          The script also checks if the directories already exist before creating them.

# Check if the input directory and ranges are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_directory> <s_range_start> <s_range_end> <p_range_start> <p_range_end>"
    exit 1
fi

input_directory=$1
s_range_start=$2
s_range_end=$3
p_range_start=$4
p_range_end=$5

# Create nested directories
for s in $(seq "$s_range_start" "$s_range_end"); do
    for p in $(seq -w "$p_range_start" "$p_range_end"); do
        dir0="$input_directory/S00${s}/P${p}/interproscan/"
        dir1="$input_directory/S00${s}/P${p}/validated/"

        # Create introproscan directories
        # Check if the directory already exists
        if [ -d "$dir0" ]; then
            echo "Directory '$dir0' already exists."
        else
            # Create the directory
            mkdir -p "$dir0"
            echo "Directory '$dir0' created."
        fi

        # Create validated directory
        # Check if the directory already exists
        if [ -d "$dir1" ]; then
            echo "Directory '$dir1' already exists."
        else
            # Create the directory
            mkdir -p "$dir1"
            echo "Directory '$dir1' created."
        fi
    done
done

# results a message of completion and awkward sock drawer comment. 
echo "Congratulations, script execution complete. Your directories are now more organized than my sock drawer."

