library(readxl)
library(stringr)

# -----------------------------
# Set working directory
# -----------------------------
setwd("D:/rashidul sir/nasir vai/nuceotide works/n")

# -----------------------------
# Input files and parameters
# -----------------------------
excel_file <- "polycluster.xlsx"
fna_file <- "Bglu49.fna"
search_term <- "Bglu49"
output_fasta <- paste0(search_term, "Bglu_1output-review.fasta")

# -----------------------------
# Step 1: Read the Excel file without header
# -----------------------------
df <- read_excel(excel_file, col_names = FALSE)
df <- as.data.frame(df, stringsAsFactors = FALSE)

# -----------------------------
# Step 2: Find all matches in the whole sheet
# -----------------------------
matches <- list()

for (i in 1:nrow(df)) {
  for (j in 1:ncol(df)) {
    cell_value <- df[i, j]
    if (!is.na(cell_value) && str_detect(cell_value, fixed(search_term))) {
      gene_name <- df[i, 1]  # First column = gene name
      matches[[length(matches) + 1]] <- list(
        gene = gene_name,
        cell = cell_value
      )
    }
  }
}

cat("Total matches found:", length(matches), "\n")  # Debug

# -----------------------------
# Step 3: Extract accession numbers
# -----------------------------
accessions <- data.frame(gene = character(), acc = character(), stringsAsFactors = FALSE)

for (m in matches) {
  gene <- m$gene
  items <- unlist(str_split(m$cell, "[;,\\s]+"))
  for (item in items) {
    if (str_detect(item, "\\|")) {
      acc <- str_split(item, "\\|")[[1]][2]
      if (!is.na(acc) && acc != "") {
        accessions <- rbind(accessions, data.frame(gene = gene, acc = acc, stringsAsFactors = FALSE))
      }
    }
  }
}

# Remove duplicates
accessions <- unique(accessions)
cat("Total unique accession numbers extracted:", nrow(accessions), "\n")  # Debug

# -----------------------------
# Step 4: Read the .fna file
# -----------------------------
fna_lines <- readLines(fna_file)

# -----------------------------
# Step 5: Extract sequences for matching accession numbers
# -----------------------------
fasta_out <- character()
for (i in 1:nrow(accessions)) {
  gene <- accessions$gene[i]
  acc <- accessions$acc[i]
  
  # Pattern to match in header: protein_id=ACCESSION
  pattern <- paste0("protein_id=", acc)
  
  header_index <- grep(pattern, fna_lines, ignore.case = TRUE)
  
  if (length(header_index) > 0) {
    header <- fna_lines[header_index[1]]
    seq_lines <- character()
    idx <- header_index[1] + 1
    while (idx <= length(fna_lines) && !startsWith(fna_lines[idx], ">")) {
      seq_lines <- c(seq_lines, fna_lines[idx])
      idx <- idx + 1
    }
    sequence <- paste(seq_lines, collapse = "")
    
    # Custom header
    custom_header <- paste0(">", gene, "|", search_term, "|", acc)
    
    
    fasta_out <- c(fasta_out, custom_header, sequence)
  } else {
    cat("WARNING: Accession", acc, "not found in FNA file\n")  # Debug
  }
}

# -----------------------------
# Step 6: Write output FASTA
# -----------------------------
if (length(fasta_out) > 0) {
  writeLines(fasta_out, output_fasta)
  cat("FASTA file created:", output_fasta, "\n")
} else {
  cat("No sequences were extracted. Please check accession names and FNA headers.\n")
}
