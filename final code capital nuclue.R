library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

setwd("D:/rashidul sir/nasir vai/nuceotide works/n")

# Inputs
excel_file <- "flagellacluster.xlsx"
gbff_file <- "Bglu.gbff"
isolate_name <- "BD_21G"
protein_fasta <- "Bglu_proteins1.fasta"
nucleotide_fasta <- "Bglu_nucleotides1.fasta"

# 1️⃣ Read Excel and clean
tox_data <- read_excel(excel_file) %>%
  mutate_all(~str_trim(as.character(.)))

# 2️⃣ Flatten all cells and extract protein IDs
all_cells <- tox_data %>%
  mutate(row_gene = tox_data[[1]]) %>%
  pivot_longer(-row_gene, names_to="column", values_to="cell") %>%
  filter(!is.na(cell)) %>%
  rowwise() %>%
  mutate(parts = list(unlist(str_split(cell, "[;,\t\n ]")))) %>%
  unnest(cols = c(parts)) %>%
  mutate(parts = str_trim(parts)) %>%
  filter(parts != "")

protein_gene_map <- all_cells %>%
  filter(str_detect(parts, regex("Bglu\\|WP_[0-9]+", ignore_case = TRUE))) %>%
  mutate(protein_id = str_remove(parts, "Bglu\\|")) %>%
  select(protein_id, gene_name = row_gene) %>%
  distinct()

protein_ids <- protein_gene_map$protein_id
cat("Total unique protein IDs found:", length(protein_ids), "\n")

# 3️⃣ Read GBFF file
gb_lines <- readLines(gbff_file)

# Extract ORIGIN sequence
origin_start <- grep("^ORIGIN", gb_lines)
origin_seq <- paste(gb_lines[(origin_start+1):length(gb_lines)], collapse = "")
origin_seq <- gsub("[^acgt]", "", origin_seq)

# ✅ Function: Wrap sequence into lines of fixed width
wrap_fasta <- function(seq, width = 60){
  return(paste(strwrap(seq, width = width), collapse = "\n"))
}

# ✅ Function: Extract protein sequence
extract_protein_sequence <- function(pid, gb_lines){
  match_line <- grep(paste0('/protein_id="', pid, '"'), gb_lines)
  if(length(match_line) == 0) return(NULL)
  
  trans_start <- grep('/translation="', gb_lines)
  trans_start <- trans_start[trans_start > match_line[1]][1]
  if(is.na(trans_start)) return(NULL)
  
  seq_lines <- c()
  for(i in trans_start:length(gb_lines)){
    line <- gb_lines[i]
    seq_lines <- c(seq_lines, line)
    if(grepl('"$', line)) break
  }
  
  full_text <- paste(seq_lines, collapse = "")
  prot_seq <- gsub('.*?/translation="|"$', '', full_text)
  prot_seq <- gsub("\\s+", "", prot_seq)
  return(prot_seq)
}

# ✅ Function: Extract CDS nucleotide sequence
extract_cds_sequence <- function(pid, gb_lines, origin_seq){
  start_line <- grep(paste0('/protein_id="', pid, '"'), gb_lines)
  if(length(start_line) == 0) return(NULL)
  
  # Find CDS location
  cds_start <- max(grep("^ {5}CDS", gb_lines[1:start_line[1]]))
  cds_line <- gb_lines[cds_start]
  
  location <- gsub(".*CDS +", "", cds_line)
  
  # Extract numeric positions
  positions <- as.numeric(unlist(str_extract_all(location, "\\d+")))
  if(length(positions) < 2) return(NULL)
  
  seq <- ""
  for(i in seq(1, length(positions), by=2)){
    seq <- paste0(seq, substr(origin_seq, positions[i], positions[i+1]))
  }
  
  # Reverse complement if complement()
  if(grepl("complement", location)){
    seq <- paste(rev(strsplit(seq, "")[[1]]), collapse = "")
    seq <- chartr("acgt", "tgca", seq)
  }
  
  seq <- toupper(seq)  # ✅ Convert to uppercase
  return(seq)
}

# ✅ Save FASTA files
prot_conn <- file(protein_fasta, "w")
nuc_conn <- file(nucleotide_fasta, "w")

for(pid in protein_ids){
  prot_seq <- extract_protein_sequence(pid, gb_lines)
  nuc_seq <- extract_cds_sequence(pid, gb_lines, origin_seq)
  gene_name <- protein_gene_map$gene_name[protein_gene_map$protein_id == pid]
  
  if(!is.null(prot_seq) & !is.null(nuc_seq)){
    header <- paste0(">", gene_name, "|", isolate_name, "|", pid)
    
    writeLines(header, prot_conn)
    writeLines(wrap_fasta(prot_seq), prot_conn)
    
    writeLines(header, nuc_conn)
    writeLines(wrap_fasta(nuc_seq), nuc_conn)
  } else {
    cat("Missing sequence for:", pid, "\n")
  }
}

close(prot_conn)
close(nuc_conn)
cat("Protein FASTA:", protein_fasta, "\n")
cat("Nucleotide FASTA:", nucleotide_fasta, "\n")
