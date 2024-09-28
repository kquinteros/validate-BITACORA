# filter-size-selection-function.R
# Author: Kevin Quinteros
# Date: 2023-12-31
# Description: A function to filter FASTA sequences within a specified size range.

#' Check If Protein Sequences are within a specified size range
#'
#' This function reads in protein fasta sequences 
#' and filters FASTA sequences based on specified size range
#'
#' @param input_fasta_path Path to the input FASTA file.
#' @param output_fasta_path Path to save the filtered FASTA file.
#' @param interproscan_path Path to the InterProScan TSV file.
#' @param ACC_ID Vector of ACC_ID values to filter the sequences.
#' @param results_header CSV file specifying column names for the InterProScan file.
#'
#' @return A filtered FASTA file containing sequences with the correct size range

# Define function

# Helper function to check if the sequence length is within the specified range
is_within_range <- function(seq, min_length, max_length) {
  seq_length <- length(seq)
  return(all(seq_length >= min_length & seq_length <= max_length))
}

# Perform size selection between the specified range
filter_size_selection <- function(input_fasta_path, output_fasta_path, min_length, max_length) {
  # Check if the input FASTA file exists
  if (!file.exists(input_fasta_path)) {
    stop("ERROR: Input FASTA file not found.")
  }

  # Read in the FASTA file
  fasta_sequences <- Biostrings::readAAStringSet(input_fasta_path)

  # Perform size selection between the specified range
  selected_sequences <- fasta_sequences[sapply(fasta_sequences, is_within_range, min_length = min_length, max_length = max_length)]

  # Make sure that all sequences are unique
  unique_sequences <- unique(selected_sequences)

  # Print a message with Run information
  print(paste("Before selection:", length(fasta_sequences), "After selection:", length(selected_sequences), "Number of Unique Sequences after selection:", length(unique_sequences)))


  # Check if the selected sequences are empty
  if (length(unique_sequences) == 0) {
    stop("ERROR: No sequences found within the specified size range.")
  }

  # Write out the selected FASTA sequences
  Biostrings::writeXStringSet(unique_sequences, file = output_fasta_path, format = "fasta", width = 80)

  # Check if the output file is empty
  if (file.size(output_fasta_path) == 0) {
    stop("ERROR: Output FASTA file is empty. No valid sequences found.")
  }
}



