#!/usr/bin/env Rscript

cat("\n========================================\n")
cat("TEM-1 Expression Model Only\n")
cat("========================================\n\n")

setwd("/courses/BINF6310.202610/students/lo.mich/final_group_project")
if (!dir.exists("output")) { dir.create("output") }

# Install packages if needed
pkgs <- c("readr", "readxl", "dplyr", "tidyr", "ape", "phytools", 
          "MCMCglmm", "coda", "parameters", "bayestestR", "reshape2", 
          "ggplot2", "RColorBrewer", "forcats")

for (p in pkgs) {
  if (!require(p, character.only = TRUE)) {
    cat(sprintf("Installing %s...\n", p))
    install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

cat("\nRunning modelExpression.R...\n")
source("modelExpression.R")

cat("\n========================================\n")
cat("Expression model complete!\n")
cat("========================================\n\n")
