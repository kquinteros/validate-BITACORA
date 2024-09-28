# filter-protein-function.R
# Author: Kevin Quinteros
# Date: 2023-12-18
# Description: Function for checking protein domains and filtering FASTA sequences.

#' Check Protein Domain and Filter FASTA Sequences
#'
#' This function checks protein domain information from an InterProScan file
#' and filters FASTA sequences based on specified ACC_ID values.
#'
#' @param input_fasta_path Path to the input FASTA file from BITACORA output.
#' @param output_fasta_path Path to save the filtered FASTA file.
#' @param interproscan_path Path to the InterProScan TSV file.
#' @param ACC_ID Vector of ACC_ID values to filter the sequences.
#' @param results_header CSV file specifying column names for the InterProScan file.
#'
#' @return A filtered FASTA file containing sequences with the correct ACC_ID.
#'
#' @examples
#' \dontrun{
#' Rscript your_script.R input.fasta output_filtered.fasta interproscan.tsv PF123,PF456
#' }
#'
#' @seealso
#' \code{\link{read.table}}, \code{\link{seqinr::read.fasta}}, \code{\link{seqinr::write.fasta}}
#'
#' @importFrom tools file_ext
#' @importFrom seqinr read.fasta write.fasta
#'
#' @export

filter_protein_domain <- function(input_fasta_path, output_fasta_path, interproscan_path, ACC_ID) {
    # Check if the required packages are installed
    if (!requireNamespace("tools", quietly = TRUE) || !requireNamespace("Biostrings", quietly = TRUE)) {
        stop("Error: The required package(s) is/are not installed. Please install them using BiocManager::install().")
    }

    # Default value for results.header
    results_header <- c("ID", "md5", "seq length", "Library", "ACC", "family", "hmmer3-location.start", "hmmer3-location.end", "hmmer3.Evalue", "Post.processed", "Date-run", "ACC.entry", "Type", "V1", "V2")

    # Validate interproscan file is a tsv file
    if (!tools::file_ext(interproscan_path) %in% c("tsv")) {
        stop("ERROR: Failed to read interproscan file. Check extension.")
    }

    # Read interproscan file
    interproscan <- read.table(interproscan_path, header = FALSE, sep = "\t", col.names = results_header)
    print("Successfully read interproscan file.")

    # Validate output file extension is fasta
    if (!tools::file_ext(output_fasta_path) %in% c("fasta", "fas", "fa")) {
        stop("ERROR: Output fasta filename does not include fasta extension.")
    }

    # Split ACC_ID by comma and convert to vector
    ACC_ID <- unlist(strsplit(ACC_ID, ","))
    print(ACC_ID)
    
    # Filter interproscan data based on ACC_ID
    filtered_interproscan <- interproscan[interproscan$ACC %in% ACC_ID, ]

    # Read in fasta from BITACORA output
    print(paste("Reading fasta file", input_fasta_path))
    fasta <- Biostrings::readAAStringSet(input_fasta_path, format = "fasta")

    # Remove predicted protein sequences with incorrect protein domain signature
    invalid_entry_indices <- which(!(names(fasta) %in% filtered_interproscan$ID))

    # Remove these sequences
    if (length(invalid_entry_indices) == 0){
        filtered_sequences <- fasta
    } else{
        filtered_sequences <- fasta[-invalid_entry_indices]
    }

    # Check if the output file is empty
    if (length(filtered_sequences) == 0) {
        stop("ERROR: Filtered sequences are empty. No valid sequences found.")
    }

    # Print the number of unique protein domains
    print(paste("Multiple protein domains detected. Removing non-matching protein domains", unique(filtered_interproscan$ACC)))

    # Write to file
    Biostrings::writeXStringSet(filtered_sequences, file = output_fasta_path, format = "fasta", width=80)

    # Check if the output file is empty
    if (file.size(output_fasta_path) == 0) {
        stop("ERROR: Output fasta file is empty. No valid sequences found.")
    } else {
        print("Fasta sequences successfully written to file.")
    }
}