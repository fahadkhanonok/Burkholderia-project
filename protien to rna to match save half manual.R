# -------------------------
# Set working directory
# -------------------------
setwd("D:/rashidul sir/nasir vai/nuceotide works/n")

# -------------------------
# Input and output files
# -------------------------
genome_fasta <- "wholeg21g.fna"       # Whole genome FASTA file
output_fasta <- "nucleotide_output21g.fasta"

# -------------------------
# CDS info for peh genes
# -------------------------
cds_info <- data.frame(
  gene_name = c(
    "pehA", "pehA", "pehB"
  ),
  protein_id = c(
    "WP_012733522.1", "WP_017922174.1", "WP_017423921.1"
  ),
  start = c(
    1126375, 1602190, 3212278
  ),
  end = c(
    1127763, 1603569, 3213996
  ),
  strand = c(
    "-", "+", "-"
  ),
  stringsAsFactors = FALSE
)

# -------------------------
# Function to read genome
# -------------------------
read_genome <- function(fasta_file) {
  lines <- readLines(fasta_file)
  lines <- lines[!grepl("^>", lines)]  # remove header lines
  genome_seq <- toupper(paste(lines, collapse = ""))
  return(genome_seq)
}

# -------------------------
# Reverse complement function
# -------------------------
revcomp <- function(seq) {
  seq <- strsplit(seq, "")[[1]]
  seq <- rev(seq)
  seq <- chartr("ACGTacgt", "TGCATGCA", seq)
  paste(seq, collapse = "")
}

# -------------------------
# Read genome
# -------------------------
genome_seq <- read_genome(genome_fasta)

# -------------------------
# Extract sequences
# -------------------------
nuc_seqs <- list()
for (i in 1:nrow(cds_info)) {
  start_pos <- cds_info$start[i]
  end_pos <- cds_info$end[i]
  strand <- cds_info$strand[i]
  seq <- substr(genome_seq, start_pos, end_pos)
  if (strand == "-") seq <- revcomp(seq)
  nuc_seqs[[i]] <- seq
}

# -------------------------
# Write FASTA with modified header
# Format: >gene_name|BD21G|protein_id
# -------------------------
write_fasta <- function(seqs, cds_info, file_out) {
  con <- file(file_out, "w")
  for (i in 1:length(seqs)) {
    header <- paste0(">", cds_info$gene_name[i], "|BD21G|", cds_info$protein_id[i])
    cat(header, "\n", file = con)
    seq_lines <- strwrap(seqs[[i]], width = 60)
    for (line in seq_lines) {
      cat(line, "\n", file = con)
    }
  }
  close(con)
}

# -------------------------
# Run FASTA writing
# -------------------------
write_fasta(nuc_seqs, cds_info, output_fasta)
cat("✅ Nucleotide FASTA saved as:", output_fasta, "\n")
