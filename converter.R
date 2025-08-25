# -------------------------
# Input files
# -------------------------
genome_fasta <- "wholeg21g.fna"        # Your whole genome FASTA
output_fasta <- "nucleotide_output21g.fasta"

# -------------------------
# CDS info for multiple proteins
# start and end are genome coordinates
# strand: "+" forward, "-" complement
# -------------------------
cds_info <- data.frame(
  protein_id = c("WP_012733585.1", "WP_251107590.1"),
  start = c(1229215, 1230291),
  end = c(1230291, 1231349),
  strand = c("+", "+"),  # "+" because no complement() in coordinates
  stringsAsFactors = FALSE
)

# -------------------------
# Function to read genome FASTA
# -------------------------
read_genome <- function(fasta_file){
  lines <- readLines(fasta_file)
  lines <- lines[!grepl("^>", lines)]   # remove header lines
  genome_seq <- toupper(paste(lines, collapse = ""))  # concatenate and uppercase
  return(genome_seq)
}

# -------------------------
# Function to reverse-complement
# -------------------------
revcomp <- function(seq){
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
# Extract nucleotide sequences
# -------------------------
nuc_seqs <- list()
for(i in 1:nrow(cds_info)){
  pid <- cds_info$protein_id[i]
  start_pos <- cds_info$start[i]
  end_pos <- cds_info$end[i]
  strand <- cds_info$strand[i]
  
  seq <- substr(genome_seq, start_pos, end_pos)
  if(strand == "-") seq <- revcomp(seq)
  
  nuc_seqs[[pid]] <- seq
}

# -------------------------
# Write to FASTA
# -------------------------
write_fasta <- function(seqs, file_out){
  con <- file(file_out, "w")
  for(name in names(seqs)){
    cat(paste0(">", name, "\n"), file = con)
    # wrap sequence at 60 characters per line
    seq_lines <- strwrap(seqs[[name]], width = 60)
    for(line in seq_lines){
      cat(line, "\n", file = con)
    }
  }
  close(con)
}

write_fasta(nuc_seqs, output_fasta)
cat("✅ Nucleotide FASTA saved as:", output_fasta, "\n")
