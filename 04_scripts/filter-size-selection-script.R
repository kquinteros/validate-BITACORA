#!/usr/bin/env Rscript

# filter-size-selection-script.R
# Author: Kevin Quinteros
# Date: 2023-12-31
# Description: Script that executes a function to filter FASTA sequences within a specified size range.

#Source function
source("04_scripts/filter-size-selection-function.R")

#commandline input
# Load library
library(argparser)

# Create an argument parser
parser <- arg_parser(description = "Protein sequence size selection")

# Define command-line arguments
parser <- add_argument(parser, "--input", help = "Input FASTA file path", type = "character")
parser <- add_argument(parser, "--output", help = "Output filtered FASTA file path", type = "character")
parser <- add_argument(parser, "--min", help = "Minimum length for Protein Sequence", type = "character")
parser <- add_argument(parser, "--max", help = "Maximum length for Protein Sequence", type = "character")

# Parse command-line arguments
args <- parse_args(parser)

# Call the filter_protein_domain function with the parsed arguments
filter_size_selection(
  input_fasta_path = args$input,
  output_fasta_path = args$output,
  min_length = args$min,
  max_length = args$max
)