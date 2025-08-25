library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(Biostrings)

setwd("D:/rashidul sir/nasir vai/nuceotide works/n")

# Read Excel
tox_data <- read_excel("toxcluster.xlsx") %>% mutate_all(~str_trim(as.character(.)))

# Flatten all cells
all_cells <- tox_data %>%
  pivot_longer(everything(), names_to="column", values_to="cell") %>%
  filter(!is.na(cell)) %>%
  rowwise() %>%
  mutate(parts = list(unlist(str_split(cell, "[;,\t\n ]")))) %>%
  unnest(cols = c(parts)) %>%
  mutate(parts = str_trim(parts)) %>%
  filter(parts != "")

# Extract protein IDs
matches <- all_cells %>%
  filter(str_detect(parts, regex("Bglu\\|WP_[0-9]+", ignore_case = TRUE))) %>%
  mutate(protein_id = str_remove(parts, "Bglu\\|"))

protein_ids <- unique(matches$protein_id)
# Print unique protein IDs
print(protein_ids)

# Or, for a nicer format
cat("Unique protein IDs found:\n")
cat(protein_ids, sep = "\n")
####
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

setwd("D:/rashidul sir/nasir vai/nuceotide works/n")

# List of protein IDs
protein_ids <- c(
  "WP_042967738.1","WP_012733473.1","WP_012733474.1",
  "WP_012734873.1","WP_039200734.1","WP_012733469.1",
  "WP_012733468.1","WP_012733014.1","WP_012733467.1",
  "WP_012733015.1","WP_015876878.1","WP_230674341.1",
  "WP_012733464.1","WP_012733470.1"
)

isolate_name <- "BD_21G"

# Read Excel file
tox_data <- read_excel("toxcluster.xlsx") %>%
  mutate_all(~str_trim(as.character(.)))

# Flatten all cells with column info
all_cells <- tox_data %>%
  mutate(row_gene = tox_data[[1]]) %>%  # First column = gene name
  pivot_longer(-row_gene, names_to = "column", values_to = "cell") %>%
  filter(!is.na(cell)) %>%
  rowwise() %>%
  mutate(parts = list(unlist(str_split(cell, "[;,\t\n ]")))) %>%
  unnest(cols = c(parts)) %>%
  mutate(parts = str_trim(parts)) %>%
  filter(parts != "")

# Map protein_id to gene name
protein_gene_map <- all_cells %>%
  filter(str_detect(parts, paste(protein_ids, collapse = "|"))) %>%
  mutate(protein_id = str_remove(parts, "Bglu\\|")) %>%
  select(protein_id, gene_name = row_gene)

print(protein_gene_map)  # check mapping

# Read GBFF file
gb_lines <- readLines("Bglu.gbff")

# Function to extract protein sequence
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

# Save sequences to FASTA
fasta_file <- "Bglu_proteins_BD21G_gene_first.fasta"
fileConn <- file(fasta_file, "w")

for(pid in protein_ids){
  seq <- extract_protein_sequence(pid, gb_lines)
  gene_name <- protein_gene_map$gene_name[protein_gene_map$protein_id == pid]
  if(!is.null(seq) & length(gene_name) > 0){
    # Header: gene_name|BD_21G|protein_id
    writeLines(paste0(">", gene_name, "|", isolate_name, "|", pid), fileConn)
    writeLines(seq, fileConn)
  } else {
    cat("Protein ID not found or gene name missing:", pid, "\n")
  }
}

close(fileConn)
cat("FASTA file saved as", fasta_file, "\n")

