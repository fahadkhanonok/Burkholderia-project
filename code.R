library(readxl)
library(Biostrings)

# Set working directory
setwd("D:/rashidul sir/nasir vai/nuceotide works/n")

# Create output folder for lipasecluster results
output_dir <- "lipase_matched_sequences"
if(!dir.exists(output_dir)) dir.create(output_dir)

# List all FFN files
ffn_files <- list.files(pattern = "\\.ffn$")

# Read lipasecluster Excel sheet
lipasecluster <- read_excel("lipasecluster.xlsx", col_names = FALSE)
lip_matrix <- as.matrix(lipasecluster)

# Extract row names from first column (e.g., lipA, lipB, etc.)
row_names <- trimws(as.character(lip_matrix[,1]))

# Clean Excel cells (remove special characters except letters, numbers, |, _)
lip_clean <- apply(lip_matrix, c(1,2), function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[^[:alnum:]|_]", "", x)
  return(x)
})

# Loop over all FFN files
for(ffn_file in ffn_files){
  
  ffn_name <- sub("\\.ffn$", "", ffn_file)
  cat("\nProcessing FFN file:", ffn_file, "\n")
  
  # Read FFN sequences
  ffn_seq <- readDNAStringSet(ffn_file)
  
  all_matched_seq <- DNAStringSet()  # To collect sequences for all rows
  
  # Loop over rows in Excel (lipA, lipB, etc.)
  for(i in 1:nrow(lip_clean)){
    
    # Get all cells in this row
    row_cells <- lip_clean[i, ]
    
    # Find matches for this FFN name in this row
    matches <- row_cells[grepl(ffn_name, row_cells, ignore.case = TRUE)]
    
    if(length(matches) == 0) next  # skip if no match
    
    # Extract IDs after "|"
    ids <- sapply(strsplit(matches, "\\|"), function(x) trimws(x[2]))
    
    # Retrieve sequences matching IDs (case-insensitive, handles extra header text)
    matched_seq <- ffn_seq[sapply(names(ffn_seq), function(h) any(tolower(ids) %in% tolower(strsplit(h, " ")[[1]])))]
    
    if(length(matched_seq) == 0) next
    
    # Modify headers: >RowName|FFN_name|ID
    new_names <- sapply(names(matched_seq), function(h) {
      matched_id <- ids[tolower(ids) %in% tolower(strsplit(h, " ")[[1]])]
      paste0(row_names[i], "|", ffn_name, "|", matched_id)
    })
    names(matched_seq) <- new_names
    
    # Add to combined matched sequences
    all_matched_seq <- c(all_matched_seq, matched_seq)
  }
  
  if(length(all_matched_seq) == 0){
    cat("No sequences found in", ffn_file, "for any row\n")
    next
  }
  
  # Save matched sequences to new folder
  output_file <- file.path(output_dir, paste0(ffn_name, "_matched_sequences.ffn"))
  writeXStringSet(all_matched_seq, output_file)
  
  cat("Saved", length(all_matched_seq), "sequences with row-specific headers to", output_file, "\n")
}
