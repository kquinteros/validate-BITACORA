# filter-tran-domains.R
# Author: Kevin Quinteros
# Date: 2024-01-03
# Description: A function to filter protein FASTA sequences based on transmembrane domain counts.

#' Filter Protein Sequences Based on Transmembrane Domain Counts
#'
#' This function reads in protein FASTA sequences and filters them based on transmembrane domain counts.
#' Sequences are divided into two groups: one with a specified number of transmembrane domains
#' and the other with a different number of transmembrane domains.
#'
#' @param fasta_file Path to the input protein FASTA file.
#' @param tmbed_counts_file Path to the file containing transmembrane domain counts.
#' @param equal_output_file Path to save the FASTA file with the specified number of transmembrane domains.
#' @param less_output_file Path to save the FASTA file with a different number of transmembrane domains.
#' @param target_tm_count The desired number of transmembrane domains (default is 7).
#'
#' @return Two filtered FASTA files based on the transmembrane domain counts.

# Example Usage:
#' \code{filter_sequences(
#'   fasta_file = "path/to/input.fasta",
#'   tmbed_counts_file = "path/to/tmbed_counts.txt",
#'   equal_output_file = "path/to/equal_output.fasta",
#'   less_output_file = "path/to/less_output.fasta"
#' )}

filter_trans_domain <- function(fasta_file, tmbed_counts_file, target_output_file, target_tm_count = 7) {
  # Read fasta file
  fasta <- seqinr::read.fasta(fasta_file)
  
  # Read tmbed helix counts
  pred <- read.table(tmbed_counts_file, sep = ",", header = FALSE)
  
  # Check if files were read successfully
  if (is.null(fasta) || ncol(pred) != 2) {
    stop("Error: Failed to read input files.")
  }
  
  # Select proteins that have the target number of transmembrane domains
  index.equal <- which(pred$V2 == target_tm_count)
  
  # Select proteins that have less or more than the target number of transmembrane domains
  index.diff <- which(pred$V2 != target_tm_count)
  
  # Subset the fasta files
  fasta.equal <- fasta[index.equal]
  fasta.diff <- fasta[index.diff]

  #create path for non target fasta file
  non_target_output_file <- gsub("target.fasta", "non_target.fasta", target_output_file)
    
  # Check if both subsets are empty
  if (length(fasta.equal) == 0 && length(fasta.diff) == 0) {
    stop("Error: No sequences matching the specified criteria.")
  }
  
  # Check if fasta.equal is not empty
  if (length(fasta.equal) > 0) {
    seqinr::write.fasta(fasta.equal,
                        names = names(fasta.equal),
                        nbchar = 80,
                        file.out = target_output_file)
  }
  
  # Check if fasta.diff is not empty
  if (length(fasta.diff) > 0) {
    seqinr::write.fasta(fasta.diff,
                        names = names(fasta.diff),
                        nbchar = 80,
                        file.out = non_target_output_file)
    
    # Check if fasta.equal is empty
    if (length(fasta.equal) == 0) {
      stop("Error: No sequences matching the specified criteria.")
    }
  }
  
  # Check if the output files are empty
  if (file.size(target_output_file) == 0 ) {
    stop("Error: Failed to write target fasta output files.")
  }
  
  cat("Filtering completed successfully.\n")
}


