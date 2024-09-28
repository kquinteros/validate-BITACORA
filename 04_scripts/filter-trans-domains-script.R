#!/usr/bin/env Rscript

# filter-trans-domains-script.R
# Author: Kevin Quinteros
# Date: Sys.Date()
# Description: Script that executes a function to filter protein FASTA sequences with a target number of transmembrane protein domains

# Load libraries
library(argparser)
source("04_scripts/filter-trans-domains-function.R")

# Command-line input
# Create an argument parser
parser <- arg_parser(description = "Filter protein sequences by the number of transmembrane proteins")

# Define command-line arguments
parser <- add_argument(parser, "--fasta", help = "Path to the input protein FASTA file", type = "character")
parser <- add_argument(parser, "--pred", help = "Path to the file containing transmembrane domain counts", type = "character")
parser <- add_argument(parser, "--target", help = "Target output FASTA file path", type = "character")
parser <- add_argument(parser, "--number", help = "The desired number of transmembrane domains", type = "character")

# Parse command-line arguments
args <- parse_args(parser)

# Call the filter_trans_domain function with the parsed arguments
filter_trans_domain(
  fasta_file = args$fasta,
  tmbed_counts_file = args$pred,
  target_output_file = args$target,
  target_tm_count = args$number
)
