#!/usr/bin/env Rscript

# filter-protein-scriot.R
# Author: Kevin Quinteros
# Date: 2023-12-18
# Description: Script that runs function for checking protein domains and filtering FASTA sequences.

#commandline input

# filter-protein-function.R
source("04_scripts/filter-protein-function.R")

# Load library
library(argparser)

# Create an argument parser
parser <- arg_parser(description = "Filter Protein Domain")

# Define command-line arguments
parser <- add_argument(parser, "--input", help = "Input FASTA file path", type = "character")
parser <- add_argument(parser, "--output", help = "Output filtered FASTA file path", type = "character")
parser <- add_argument(parser, "--interproscan", help = "InterProScan TSV file path", type = "character")
parser <- add_argument(parser, "--acc_ids", help = "Comma-separated ACC_ID values", type = "character")

# Parse command-line arguments
args <- parse_args(parser)

# Call the filter_protein_domain function with the parsed arguments
filter_protein_domain(
  input_fasta_path = args$input,
  output_fasta_path = args$output,
  interproscan_path = args$interproscan,
  ACC_ID = args$acc_ids
)