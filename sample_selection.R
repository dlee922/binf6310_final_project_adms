library(tidyverse)

# Load BOTH files
cat("Select File 1: CSV with assembly accessions (MOESM7)\n")
assembly_file <- read.csv(file.choose())

cat("\nSelect File 2: TXT with full metadata (MOESM4)\n")
metadata_file <- read.delim(file.choose(), sep = "\t")

# Filter metadata to chromosomes
metadata_clean <- metadata_file %>%
  filter(contig.type == "chromosome") %>%
  filter(!is.na(ezclermont.phylogroup)) %>%
  filter(!is.na(mlst.st))

# Merge
combined <- assembly_file %>%
  inner_join(metadata_clean, by = "isolate.id") %>%
  filter(!is.na(genome.assembly.accession))

cat("\nCombined isolates:", nrow(combined), "\n\n")

# === SELECT MANUALLY ===
set.seed(789)

# Step 1: Get 2 from each phylogroup
phylo_reps <- data.frame()
for (pg in unique(combined$ezclermont.phylogroup)) {
  samples_from_pg <- combined %>% filter(ezclermont.phylogroup == pg)
  n_to_take <- min(2, nrow(samples_from_pg))
  sampled <- samples_from_pg %>% slice_sample(n = n_to_take)
  phylo_reps <- bind_rows(phylo_reps, sampled)
}

cat("Phylogroup reps:", nrow(phylo_reps), "\n")

# Step 2: Get key STs
key_sts_pool <- combined %>%
  filter(mlst.st %in% c(12, 69, 127)) %>%
  filter(!isolate.id %in% phylo_reps$isolate.id)

n_key_available <- nrow(key_sts_pool)
n_key_to_take <- min(20, n_key_available)

key_sts_selected <- key_sts_pool %>% slice_sample(n = n_key_to_take)

cat("Key STs added:", nrow(key_sts_selected), "\n")

# Step 3: Fill remainder
already_selected <- c(phylo_reps$isolate.id, key_sts_selected$isolate.id)
others_pool <- combined %>% filter(!isolate.id %in% already_selected)

n_remaining <- 50 - nrow(phylo_reps) - nrow(key_sts_selected)
others_selected <- others_pool %>% slice_sample(n = n_remaining)

# Combine all
final <- bind_rows(phylo_reps, key_sts_selected, others_selected)

# === CHECK ===
cat("\n=== FINAL ===\n")
cat("Total:", nrow(final), "\n\n")

cat("Phylogroups:\n")
print(table(final$ezclermont.phylogroup))

cat("\nKey STs:\n")
cat("ST12:", sum(final$mlst.st == 12), "\n")
cat("ST69:", sum(final$mlst.st == 69), "\n")
cat("ST127:", sum(final$mlst.st == 127), "\n")

# === SAVE ===
write.csv(final, "FINAL_50_genomes_v3.csv", row.names = FALSE)

write.table(final$genome.assembly.accession,
            "FINAL_accession_list_v3.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("\n✅ SAVED!\n")
cat("Location:", getwd(), "\n")

# Show first 10 accessions
cat("\nFirst 10 genome accessions:\n")
print(head(final$genome.assembly.accession, 10))