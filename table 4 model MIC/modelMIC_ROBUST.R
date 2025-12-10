############################################################
# ROBUST MIC MODEL - WITH AGGRESSIVE CHECKPOINTING
# Saves data at every processing stage
############################################################

library(ape)
library(readr)
library(MCMCglmm)
library(phytools)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(parameters)
library(coda)
library(stringr)

cat("\n===== MIC MODEL (WITH CHECKPOINTING) STARTING =====\n")

# -------------------------------------------------------------------
# CHECKPOINT & OUTPUT SETTINGS
# -------------------------------------------------------------------

checkpoint_dir <- "output"
if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir)

checkpoint_file_1   <- file.path(checkpoint_dir, "mic_chain1_checkpoint.rds")
checkpoint_file_2   <- file.path(checkpoint_dir, "mic_chain2_checkpoint.rds")
final_chains_file   <- file.path(checkpoint_dir, "mic_chains_final.rds")
final_results_file  <- file.path(checkpoint_dir, "mic_model_results.csv")

# -------------------------------------------------------------------
# READ IN DATA
# -------------------------------------------------------------------

cat("Loading data files...\n")

df <- read_delim("data/supplementary_data_1.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)

tem1.report <- read_delim("data/supplementary_data_2.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

phylo <- read.tree("data/ml_tree.txt")

cat("✓ Data loaded successfully.\n")

# -------------------------------------------------------------------
# PREPARE PHYLOGENY
# -------------------------------------------------------------------

cat("Preparing phylogeny...\n")

phylo$tip.label <- as.numeric(phylo$tip.label)
phylo$node.label <- NULL
phylo <- midpoint.root(phylo)
phylo.u <- chronos(phylo, lambda=1, model="correlated")

cat("  Rooted:", is.rooted(phylo.u), "\n")
cat("  Ultrametric:", is.ultrametric(phylo.u), "\n")

inv.phylo <- inverseA(phylo.u, nodes="TIPS", scale=TRUE)

cat("✓ Phylogeny prepared.\n\n")

# -------------------------------------------------------------------
# PREPARE METADATA - WITH CHECKPOINTING AT EACH STAGE
# -------------------------------------------------------------------

cat("Preparing metadata with stage-wise checkpoints...\n")

df.model <- df %>% 
  mutate(tem1.contig.copy.number = contig.copy.number * tem1.replicon) %>%
  group_by(isolate.id) %>%
  mutate(tem1.isolate.copy.number = sum(tem1.contig.copy.number)) %>%
  ungroup() %>%
  mutate(
    tem1.isolate.copy.number.scaled = 
      (tem1.isolate.copy.number - mean(tem1.isolate.copy.number)) / sd(tem1.isolate.copy.number),
    tem1.isolate.scaled = ifelse(tem1.isolate > 1, TRUE, FALSE)
  ) %>%
  mutate(tem1.isolate.scaled = as.factor(tem1.isolate.scaled)) %>%
  filter(contig.type == "chromosome") %>%
  ungroup() %>%
  select(isolate.id, isolate.assembly, coamox.mic, tem1.isolate.scaled, 
         tem1.isolate.copy.number.scaled, ampc.promoter.snv) %>%
  mutate(
    coamox.mic = factor(
      coamox.mic, 
      ordered = TRUE, 
      levels = c("<=2.2", "4.2", "8.2", "16.2", "32.2", ">32.2")
    )
  )

# Cap at 95th percentile
upper_percentile_value <- quantile(df.model$tem1.isolate.copy.number.scaled, 0.95)
df.model$tem1.isolate.copy.number.scaled[
  df.model$tem1.isolate.copy.number.scaled > upper_percentile_value
] <- upper_percentile_value

cat("  ✓ Checkpoint 1: Initial data (n=", nrow(df.model), ")\n")
saveRDS(df.model, file.path(checkpoint_dir, "01_df_initial.rds"))
write.csv(df.model, file.path(checkpoint_dir, "01_df_initial.csv"))

# Promoter SNV filtering
tem1.model <- tem1.report %>%
  group_by(isolate.assembly) %>%
  filter(!is.na(promoter.start)) %>%
  mutate(n.linked = n()) %>%
  filter(n.linked == tem1.isolate) %>%
  group_by(isolate.assembly) %>%
  mutate(promoter.snv.types = n_distinct(promoter.snv)) %>%
  filter(promoter.snv.types == 1) %>%
  group_by(promoter.snv) %>%
  mutate(promoter.snv.n = n()) %>%
  filter(promoter.snv.n >= 10) %>%
  mutate(promoter.snv = as.factor(toupper(promoter.snv))) %>%
  select(isolate.assembly, promoter.snv) %>%
  distinct()

df.model <- merge(df.model, tem1.model, by="isolate.assembly", all.x=TRUE) %>%
  filter(!is.na(promoter.snv))

cat("  ✓ Checkpoint 2: After promoter SNV (n=", nrow(df.model), ")\n")
saveRDS(df.model, file.path(checkpoint_dir, "02_df_after_promoter_snv.rds"))
write.csv(df.model, file.path(checkpoint_dir, "02_df_after_promoter_snv.csv"))

# AmpC promoter filtering
df.model <- df.model %>%
  group_by(ampc.promoter.snv) %>%
  filter(n() >= 10)

cat("  ✓ Checkpoint 3: After AmpC filter (n=", nrow(df.model), ")\n")
saveRDS(df.model, file.path(checkpoint_dir, "03_df_after_ampc_filter.rds"))
write.csv(df.model, file.path(checkpoint_dir, "03_df_after_ampc_filter.csv"))

# Factor levels
df.model$promoter.snv <- factor(
  df.model$promoter.snv, 
  ordered = FALSE, 
  levels = c("CGGCGG", "CGGCGA", "TGGCGA", "TGGCGG")
)

df.model$ampc.promoter.snv <- factor(
  df.model$ampc.promoter.snv, 
  ordered = FALSE, 
  levels = c("GGCTCCTAGGG", "AGCTTCTAGGG", "AGCTCCTAGGG", "GATTCCTAGGG")
)

df.model <- as.data.frame(df.model)
df.model$phylo <- as.factor(df.model$isolate.assembly)

cat("  ✓ Checkpoint 4: FINAL PREPARED DATA (n=", nrow(df.model), ")\n")
saveRDS(df.model, file.path(checkpoint_dir, "04_df_model_FINAL.rds"))
write.csv(df.model, file.path(checkpoint_dir, "04_df_model_FINAL.csv"))

cat("  MIC levels:", nlevels(df.model$coamox.mic), "\n")
cat("  Ready for modeling.\n\n")

# -------------------------------------------------------------------
# PRIOR
# -------------------------------------------------------------------

cat("Setting prior...\n")

prior <- list(
  G = list(G1 = list(V = 1, nu = 0.02, alpha.mu = 0, alpha.V = 1e+3)),
  R = list(fix = 1, V = 1, nu = 0.02)
)

cat("✓ Prior set.\n\n")

# -------------------------------------------------------------------
# RUN CHAIN 1 - WITH ERROR HANDLING
# -------------------------------------------------------------------

cat("===== RUNNING CHAIN 1 =====\n")

chain.1 <- NULL
chain1_error <- FALSE

tryCatch({
  set.seed(1)
  chain.1 <- MCMCglmm(
    coamox.mic ~ tem1.isolate.copy.number.scaled + tem1.isolate.scaled + 
      ampc.promoter.snv + promoter.snv,
    random = ~ phylo,
    family = "ordinal",
    ginverse = list(phylo = inv.phylo$Ainv),
    prior = prior,
    data = df.model,
    nitt = 10000000,
    burnin = 1000000,
    thin = 100,
    DIC = FALSE,
    pr = TRUE,
    trunc = TRUE
  )
  
  cat("✓ Chain 1 completed successfully.\n")
  saveRDS(chain.1, checkpoint_file_1)
  cat("  Saved:", checkpoint_file_1, "\n")
  
  # Save summary
  chain1_summary <- summary(chain.1)$solutions
  write.csv(chain1_summary, file.path(checkpoint_dir, "05_chain1_summary.csv"))
  cat("  Saved: 05_chain1_summary.csv\n")
  
}, error = function(e) {
  cat("\n!!!! CHAIN 1 FAILED !!!!\n")
  cat("Error:", e$message, "\n")
  cat("All prepared data has been saved in checkpoints 01-04.\n")
  chain1_error <<- TRUE
})

# -------------------------------------------------------------------
# RUN CHAIN 2 - WITH ERROR HANDLING (only if chain 1 succeeded)
# -------------------------------------------------------------------

chain.2 <- NULL
chain2_error <- FALSE

if (!chain1_error) {
  cat("\n===== RUNNING CHAIN 2 =====\n")
  
  tryCatch({
    set.seed(2)
    chain.2 <- MCMCglmm(
      coamox.mic ~ tem1.isolate.copy.number.scaled + tem1.isolate.scaled + 
        ampc.promoter.snv + promoter.snv,
      random = ~ phylo,
      family = "ordinal",
      ginverse = list(phylo = inv.phylo$Ainv),
      prior = prior,
      data = df.model,
      nitt = 10000000,
      burnin = 1000000,
      thin = 100,
      DIC = FALSE,
      pr = TRUE,
      trunc = TRUE
    )
    
    cat("✓ Chain 2 completed successfully.\n")
    saveRDS(chain.2, checkpoint_file_2)
    cat("  Saved:", checkpoint_file_2, "\n")
    
    # Save summary
    chain2_summary <- summary(chain.2)$solutions
    write.csv(chain2_summary, file.path(checkpoint_dir, "06_chain2_summary.csv"))
    cat("  Saved: 06_chain2_summary.csv\n")
    
  }, error = function(e) {
    cat("\n!!!! CHAIN 2 FAILED !!!!\n")
    cat("Error:", e$message, "\n")
    cat("Chain 1 results are still available.\n")
    chain2_error <<- TRUE
  })
} else {
  cat("\nSkipping Chain 2 because Chain 1 failed.\n")
}

# -------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS (if both chains succeeded)
# -------------------------------------------------------------------

if (!chain1_error && !chain2_error) {
  cat("\n===== CONVERGENCE DIAGNOSTICS =====\n")
  
  cat("\n--- Chain 1 Summary ---\n")
  print(summary(chain.1))
  
  cat("\n--- Chain 1 Autocorrelation ---\n")
  print(autocorr.diag(chain.1$VCV))
  
  cat("\n--- Chain 2 Summary ---\n")
  print(summary(chain.2))
  
  cat("\n--- Chain 2 Autocorrelation ---\n")
  print(autocorr.diag(chain.2$VCV))
  
  # Gelman-Rubin diagnostic
  cat("\n--- Gelman-Rubin PSRF (should be <1.1) ---\n")
  mclist <- mcmc.list(chain.1$Sol, chain.2$Sol)
  gr <- gelman.diag(mclist)
  print(gr)
  
  # Use chain 1 for downstream analysis
  model <- chain.1
  
  cat("\n===== MODEL PARAMETERS =====\n")
  params <- model_parameters(
    model, 
    centrality = "mean", 
    ci = 0.95, 
    ci_method = "hdi",
    component = "all"
  )
  print(params)
  
  # Save results
  cat("\n===== SAVING FINAL RESULTS =====\n")
  saveRDS(list(chain.1 = chain.1, chain.2 = chain.2), final_chains_file)
  cat("✓ Saved:", final_chains_file, "\n")
  
  write.csv(summary(model)$solutions, final_results_file)
  cat("✓ Saved:", final_results_file, "\n")
  
} else if (!chain1_error) {
  cat("\n===== PARTIAL COMPLETION: CHAIN 1 ONLY =====\n")
  model <- chain.1
  print(summary(model))
  write.csv(summary(model)$solutions, final_results_file)
  cat("✓ Saved:", final_results_file, "\n")
  
} else {
  cat("\n===== FAILURE: DATA SAVED =====\n")
  cat("Both chains failed, but all prepared data is in checkpoints 01-04.\n")
  cat("Load checkpoint files to continue analysis or debugging.\n")
  model <- NULL
}

# -------------------------------------------------------------------
# EXTRACT PHYLOGENETIC EFFECTS
# -------------------------------------------------------------------

if (!is.null(model)) {
  cat("\n===== EXTRACTING PHYLOGENETIC EFFECTS =====\n")
  
  phylo.effects <- model$Sol[, grep("phylo", colnames(model$Sol))]
  phylo.effects <- as.data.frame(phylo.effects)
  phylo.effects.long <- tidyr::pivot_longer(
    phylo.effects, 
    cols = everything(), 
    names_to = "tip", 
    values_to = "value"
  )
  phylo.effects.long$tip <- gsub("phylo.", "", phylo.effects.long$tip)
  
  write.csv(phylo.effects.long, file.path(checkpoint_dir, "07_phylo_effects_full.csv"))
  cat("✓ Saved: 07_phylo_effects_full.csv\n")
  
  # Summary stats
  phylo_summary <- phylo.effects.long %>%
    group_by(tip) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      q025 = quantile(value, 0.025),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    )
  
  write.csv(phylo_summary, file.path(checkpoint_dir, "08_phylo_effects_summary.csv"))
  cat("✓ Saved: 08_phylo_effects_summary.csv\n")
}

cat("\n===== MIC MODEL COMPLETE =====\n")
cat("\nCheckpoint files saved in output/:\n")
cat("  01_df_initial.rds/.csv\n")
cat("  02_df_after_promoter_snv.rds/.csv\n")
cat("  03_df_after_ampc_filter.rds/.csv\n")
cat("  04_df_model_FINAL.rds/.csv\n")
cat("  05_chain1_summary.csv (if chain 1 ran)\n")
cat("  06_chain2_summary.csv (if chain 2 ran)\n")
cat("  07_phylo_effects_full.csv (if model ran)\n")
cat("  08_phylo_effects_summary.csv (if model ran)\n")
cat("  mic_chains_final.rds (if both chains succeeded)\n")
cat("  mic_model_results.csv (if model ran)\n")
cat("===== END =====\n\n")
EOF
