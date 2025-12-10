#!/usr/bin/env Rscript
###############################################################################
# TEM-1 Beta-Lactamase Expression & MIC Hierarchical Bayesian Models
# Reproducibility study: E. coli phylogeny and antibiotic resistance
# Nature Communications - Simplified for Tables 2 & 4
###############################################################################

# Load libraries - MINIMAL SET
library(ape)
library(readr)
library(MCMCglmm)
library(phytools)
library(dplyr)
library(tidyr)
library(coda)
library(readxl)

###############################################################################
# SECTION 1: Read in data
###############################################################################

cat("Reading data files...\n")

data_dir <- "./data"
output_dir <- "./output"
dir.create(output_dir, showWarnings = FALSE)

df <- read_delim(file.path(data_dir, "supplementary_data_1.tsv"), 
                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)

tem1.report <- read_delim(file.path(data_dir, "supplementary_data_2.tsv"), 
                          delim = "\t", escape_double = FALSE, trim_ws = TRUE)

phylo <- read.tree(file.path(data_dir, "ml_tree.txt"))

exp <- read_excel(file.path(data_dir, "supplementary_data_3.xlsx"))

cat("Data loaded successfully.\n")

###############################################################################
# SECTION 2: Prepare phylogeny
###############################################################################

cat("Preparing phylogeny...\n")

phylo$tip.label <- as.numeric(phylo$tip.label)
phylo$node.label <- NULL

phylo <- midpoint.root(phylo)
phylo.u <- phylo


inv.phylo <- inverseA(phylo.u, nodes = "TIPS", scale = FALSE)

cat("Phylogeny prepared.\n")

###############################################################################
# SECTION 3: Prepare modeling data
###############################################################################

cat("Preparing modeling data...\n")

# Filter for TEM-1 isolate and replicon
df.model <- df %>%
  filter(tem1.isolate == 1 & tem1.replicon == 1) %>%
  select(accession, isolate.assembly, contig.assembly, 
         contig.copy.number, contig.type)

tem1.report.less <- tem1.report[c("isolate.assembly", "contig.assembly", 
                                   "promoter.snv", "tem1.snv")]

df.model <- merge(df.model, tem1.report.less, 
                  by = c("isolate.assembly", "contig.assembly"),
                  all.x = TRUE)

# Promoter type from SNV
df.model$promoter.type <- ifelse(df.model$promoter.snv == "tggcgg", "TG",
                                 ifelse(df.model$promoter.snv == "tggcga", "TA",
                                 ifelse(df.model$promoter.snv == "cggcgg", "CG",
                                 ifelse(df.model$promoter.snv == "cggcga", "CA", NA))))

df.model <- df.model %>% filter(!is.na(promoter.type))

# Expression data
exp$Isolate <- sapply(strsplit(as.character(exp$Isolate), " "), function(x) x[1])
exp <- exp[c("Isolate", "delta Ct (control)")]
colnames(exp) <- c("accession", "exp")

df.model <- merge(df.model, exp, by = "accession", all.y = TRUE)
df.model <- df.model %>% 
  filter(!is.na(exp)) %>% 
  filter(!is.na(isolate.assembly)) %>% 
  distinct() %>%
  group_by(accession) %>% 
  mutate(run = row_number())

df.model$contig.copy.number <- as.numeric(df.model$contig.copy.number)

# Standardize
mean.exp <- mean(df.model$exp)
sd.exp <- sd(df.model$exp)
mean.cp <- mean(df.model$contig.copy.number)
sd.cp <- sd(df.model$contig.copy.number)

df.model <- df.model %>%
  mutate(exp.scaled = (exp - mean.exp) / sd.exp,
         contig.copy.number.scaled = (contig.copy.number - mean.cp) / sd.cp)

# Cap at 95% quantile
up.exp <- quantile(df.model$exp.scaled, 0.95)
df.model$exp.scaled[df.model$exp.scaled > up.exp] <- up.exp

# Promoter position indicators
new.cols <- data.frame(t(data.frame(strsplit(df.model$promoter.type, ""))))
df.model$pos1 <- new.cols$X1
df.model$pos2 <- new.cols$X2
df.model$pos1.bool <- ifelse(df.model$pos1 == "T", TRUE, FALSE)
df.model$pos2.bool <- ifelse(df.model$pos2 == "A", TRUE, FALSE)

# Phylogeny
df.model$phylo <- as.character(df.model$isolate.assembly)
phylo.trim <- keep.tip(phylo.u, df.model$phylo)
inv.phylo.trim <- inverseA(phylo.trim, nodes = "TIPS", scale = FALSE)

df.model$isolate.assembly <- as.factor(df.model$isolate.assembly)
df.model$tem1.snv <- as.factor(df.model$tem1.snv)
df.model <- as.data.frame(df.model)

cat("Data preparation complete.\n")
cat("Sample size:", nrow(df.model), "\n")
cat("Number of unique isolates:", n_distinct(df.model$isolate.assembly), "\n")

###############################################################################
# TABLE 2: MCMC Expression Model
###############################################################################

cat("\n=== Running Expression Model for Table 2 ===\n")

prior.1 <- list(
  G = list(
    G1 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 1e+3),
    G2 = list(V = 1, nu = 0.02)
  ),
  R = list(V = 1, nu = 0.02)
)

cat("Running chain 1 for expression...\n")
set.seed(1)
chain.1 <- MCMCglmm(
  exp.scaled ~ pos1.bool * pos2.bool + contig.copy.number.scaled,
  random = ~ phylo + isolate.assembly,
  family = "gaussian",
  ginverse = list(phylo = inv.phylo.trim$Ainv),
  data = df.model,
  prior = prior.1,
  nitt = 10000000,
  burnin = 1000000,
  thin = 100,
  DIC = FALSE,
  pr = TRUE
)

cat("Running chain 2 for expression...\n")
set.seed(2)
chain.2 <- MCMCglmm(
  exp.scaled ~ pos1.bool * pos2.bool + contig.copy.number.scaled,
  random = ~ phylo + isolate.assembly,
  family = "gaussian",
  ginverse = list(phylo = inv.phylo.trim$Ainv),
  data = df.model,
  prior = prior.1,
  nitt = 10000000,
  burnin = 1000000,
  thin = 100,
  DIC = FALSE,
  pr = TRUE
)

cat("Convergence checks:\n")
print(summary(chain.1))
mclist <- mcmc.list(chain.1$Sol, chain.2$Sol)
cat("\nGelman-Rubin diagnostic:\n")
print(gelman.diag(mclist))

# Extract Table 2
model <- chain.1
solutions <- as.data.frame(summary(model)$solutions)
solutions$Variable <- rownames(solutions)

table2 <- data.frame(
  Variable = solutions$Variable,
  Beta_coefficient_posterior_mean = solutions$post.mean,
  HPD_lower = solutions$l.95..CI,
  HPD_upper = solutions$u.95..CI,
  pMCMC = solutions$pMCMC
)

write.csv(table2, file.path(output_dir, "TABLE2_expression_model.csv"), row.names = FALSE)
cat("\nTable 2 saved to", file.path(output_dir, "TABLE2_expression_model.csv"), "\n")

# Save model
saveRDS(model, file.path(output_dir, "expression_model_fit.rds"))

###############################################################################
# TABLE 4: MCMC MIC Model
###############################################################################

cat("\n=== Running MIC Model for Table 4 ===\n")

# Filter data for MIC model
df.model$isolate.assembly <- as.character(df.model$isolate.assembly)

# Prepare phylogeny for all 377 isolates
df.mic <- df %>%
  select(isolate.assembly, contig.assembly, contig.copy.number, 
         contig.type, coamox.mic, tem1.isolate, tem1.replicon)

df.mic <- df.mic %>%
  filter(tem1.isolate == 1 & tem1.replicon == 1) %>%
  distinct()

tem1.report.for.mic <- tem1.report[c("isolate.assembly", "contig.assembly", 
                                      "promoter.snv", "tem1.snv")]

df.mic <- merge(df.mic, tem1.report.for.mic,
                by = c("isolate.assembly", "contig.assembly"),
                all.x = TRUE)

# Promoter type
df.mic$promoter.type <- ifelse(df.mic$promoter.snv == "tggcgg", "TG",
                               ifelse(df.mic$promoter.snv == "tggcga", "TA",
                               ifelse(df.mic$promoter.snv == "cggcgg", "CG",
                               ifelse(df.mic$promoter.snv == "cggcga", "CA", NA))))

df.mic <- df.mic %>% filter(!is.na(promoter.type))

df.mic$contig.copy.number <- as.numeric(df.mic$contig.copy.number)

# Standardization
mean.cp.mic <- mean(df.mic$contig.copy.number, na.rm = TRUE)
sd.cp.mic <- sd(df.mic$contig.copy.number, na.rm = TRUE)

df.mic <- df.mic %>%
  mutate(contig.copy.number.scaled = (contig.copy.number - mean.cp.mic) / sd.cp.mic)

# Binary indicators
new.cols.mic <- data.frame(t(data.frame(strsplit(df.mic$promoter.type, ""))))
df.mic$pos1 <- new.cols.mic$X1
df.mic$pos2 <- new.cols.mic$X2
df.mic$pos1.bool <- ifelse(df.mic$pos1 == "T", TRUE, FALSE)
df.mic$pos2.bool <- ifelse(df.mic$pos2 == "A", TRUE, FALSE)

# Genome copy number indicator
df.mic$genome.copy.indicator <- ifelse(df.mic$contig.copy.number > 1, 1, 0)

# Phylogeny for all isolates
df.mic$phylo <- as.character(df.mic$isolate.assembly)
phylo.full <- keep.tip(phylo.u, df.mic$phylo)
inv.phylo.full <- inverseA(phylo.full, nodes = "TIPS", scale = FALSE)

df.mic$isolate.assembly <- as.factor(df.mic$isolate.assembly)
df.mic <- as.data.frame(df.mic)

cat("MIC model sample size:", nrow(df.mic), "\n")

# Convert MIC to ordinal
df.mic$coamox.mic.ordinal <- factor(df.mic$coamox.mic, 
                                     levels = c("<=2.2", "4.2", "8.2", "16.2", "32.2", ">32.2"),
                                     ordered = TRUE)

cat("Running MIC model...\n")
set.seed(1)
chain.mic <- MCMCglmm(
  coamox.mic.ordinal ~ contig.copy.number.scaled + genome.copy.indicator + 
                       pos1.bool + pos2.bool + pos1.bool:pos2.bool,
  random = ~ phylo + isolate.assembly,
  family = "ordinal",
  ginverse = list(phylo = inv.phylo.full$Ainv),
  data = df.mic,
  prior = prior.1,
  nitt = 10000000,
  burnin = 1000000,
  thin = 100,
  DIC = FALSE,
  pr = TRUE
)

cat("MIC Model convergence check:\n")
print(summary(chain.mic))

# Extract Table 4
solutions.mic <- as.data.frame(summary(chain.mic)$solutions)
solutions.mic$Variable <- rownames(solutions.mic)

table4 <- data.frame(
  Variable = solutions.mic$Variable,
  Beta_coefficient_posterior_mean = solutions.mic$post.mean,
  HPD_lower = solutions.mic$l.95..CI,
  HPD_upper = solutions.mic$u.95..CI,
  pMCMC = solutions.mic$pMCMC
)

write.csv(table4, file.path(output_dir, "TABLE4_mic_model.csv"), row.names = FALSE)
cat("\nTable 4 saved to", file.path(output_dir, "TABLE4_mic_model.csv"), "\n")

# Save model
saveRDS(chain.mic, file.path(output_dir, "mic_model_fit.rds"))

###############################################################################
# SUMMARY
###############################################################################

cat("\n=== Analysis Complete ===\n")
cat("Output files:\n")
cat("  TABLE2_expression_model.csv\n")
cat("  TABLE4_mic_model.csv\n")
cat("  expression_model_fit.rds\n")
cat("  mic_model_fit.rds\n")
cat("\nSession info:\n")
print(sessionInfo())
