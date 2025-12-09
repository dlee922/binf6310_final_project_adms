library(ape)
library(readr)
library(MCMCglmm)
library(phytools)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
# library(ggtext)  # not used
library(tidyr)
library(parameters)
library(coda)
library(readxl)

## --------- Setup results dir and logging ----------
results_dir <- "/home/desai.ara/binf6310_final_project_adms/results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", ...)

msg("Process Check: Finished loading libraries")

## --------- read in data ----------
msg("Reading data")
df <- read_delim("/home/desai.ara/binf6310_final_project_adms/data/supplementary_data_1.tsv",
                 delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
tem1.report <- read_delim("/home/desai.ara/binf6310_final_project_adms/data/supplementary_data_2.tsv",
                          delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
phylo <- read.tree("/home/desai.ara/binf6310_final_project_adms/data/ml_tree.txt")
exp <- read_excel("/home/desai.ara/binf6310_final_project_adms/data/supplementary_data_3.xlsx")
msg("Finished Reading data")

## --------- prepare MIC data ----------
msg("Preparing MIC data")
df.mic <- df %>%
  mutate(tem1.contig.copy.number = contig.copy.number * tem1.replicon) %>%
  group_by(isolate.id) %>%
  mutate(tem1.isolate.copy.number = sum(tem1.contig.copy.number)) %>%
  ungroup() %>%
  mutate(tem1.isolate.copy.number.scaled = (tem1.isolate.copy.number - mean(tem1.isolate.copy.number)) / sd(tem1.isolate.copy.number),
         tem1.isolate.scaled = ifelse(tem1.isolate>1, TRUE, FALSE)) %>%
  mutate(tem1.isolate.scaled = as.factor(tem1.isolate.scaled)) %>%
  filter(contig.type=="chromosome") %>%
  ungroup() %>%
  select(isolate.id, isolate.assembly, coamox.mic, tem1.isolate.scaled, tem1.isolate.copy.number.scaled) %>%
  mutate(coamox.mic = factor(coamox.mic, ordered = TRUE, levels = c("<=2.2", "4.2", "8.2",
                                                                    "16.2", "32.2", ">32.2")))

upper_percentile_value <- quantile(df.mic$tem1.isolate.copy.number.scaled, 0.95, na.rm = TRUE)
df.mic$tem1.isolate.copy.number.scaled[df.mic$tem1.isolate.copy.number.scaled > upper_percentile_value] <- upper_percentile_value

tem1.model <- tem1.report %>%
  group_by(isolate.assembly) %>%
  filter(!is.na(promoter.start)) %>%
  mutate(n.linked = n()) %>%
  filter(n.linked == tem1.isolate) %>%
  group_by(isolate.assembly) %>%
  mutate(promoter.snv.types = n_distinct(promoter.snv)) %>%
  filter(promoter.snv.types==1) %>%
  group_by(promoter.snv) %>%
  mutate(promoter.snv.n = n()) %>%
  filter(promoter.snv.n >= 10) %>%
  mutate(promoter.snv = as.factor(toupper(promoter.snv))) %>%
  select(isolate.assembly, promoter.snv) %>%
  distinct()

df.mic <- merge(df.mic, tem1.model, by="isolate.assembly", all.x=TRUE) %>%
  filter(!is.na(promoter.snv))

df.mic$promoter.snv <- factor(df.mic$promoter.snv, ordered = FALSE,
                              levels = c("CGGCGG", "CGGCGA", "TGGCGA", "TGGCGG"))

df.mic$trait <- "mic"
df.mic$family <- "gaussian"
df.mic <- as.data.frame(df.mic)
df.mic$phylo <- as.factor(df.mic$isolate.assembly)

msg("Finished MIC data")

## --------- prepare expression data ----------
msg("Preparing expression data")
df.exp <- df %>%
  filter(tem1.isolate==1 & tem1.replicon==1) %>%
  select(accession, isolate.assembly, contig.assembly)

exp$Isolate <- sapply(strsplit(as.character(exp$Isolate), " "), function(x) x[1])
exp <- exp[c("Isolate", "delta Ct (control)")]
colnames(exp) <- c("accession", "exp")

df.exp <- merge(df.exp, exp, by="accession", all.y=TRUE)
df.exp <- df.exp %>% filter(!is.na(exp)) %>% filter(!is.na(isolate.assembly)) %>% distinct()
df.exp <- df.exp %>% group_by(accession) %>% mutate(run = row_number())
mean.exp <- mean(df.exp$exp, na.rm = TRUE)
sd.exp <- sd(df.exp$exp, na.rm = TRUE)
df.exp <- df.exp %>%
  mutate(exp.scaled = (exp - mean.exp) / sd.exp)
up.exp <- quantile(df.exp$exp.scaled, 0.95, na.rm = TRUE)
df.exp$exp.scaled[df.exp$exp.scaled > up.exp] <- up.exp
df.exp$trait <- "exp"
df.exp$family <- "gaussian"
msg("Finished expression data")

## --------- wrangle ----------
msg("Wrangling/model frame")
df.model <- merge(df.mic, df.exp, all=TRUE, by=c("isolate.assembly", "trait", "family"))
df.model <- df.model %>%
  group_by(isolate.assembly) %>%
  filter(!(!("mic" %in% trait) & "exp" %in% trait)) %>%
  fill(isolate.id, .direction = "downup") %>%
  fill(tem1.isolate.scaled, .direction = "downup") %>%
  fill(tem1.isolate.copy.number.scaled, .direction = "downup") %>%
  fill(promoter.snv, .direction = "downup") %>%
  fill(phylo, .direction = "downup") %>%
  ungroup() %>%
  select(isolate.id, phylo, trait, family, tem1.isolate.scaled, tem1.isolate.copy.number.scaled,
         promoter.snv, coamox.mic, exp.scaled, run)
df.model$y <- ifelse(!is.na(df.model$coamox.mic), as.factor(df.model$coamox.mic), as.numeric(df.model$exp.scaled))
df.model$trait <- as.factor(df.model$trait)
msg("Finished wrangle data")

## --------- phylogeny ----------
msg("Preparing phylogeny")
df.model$phylo <- as.character(df.model$phylo)
phylo$tip.label <- as.character(phylo$tip.label)
phylo$node.label <- NULL
phylo <- keep.tip(phylo, df.model$phylo)
phylo <- midpoint.root(phylo)
# chronos might warn; that's OK. you can adjust control if needed.
phylo.u <- tryCatch(
  chronos(phylo, lambda=1, model="correlated"),
  error = function(e){
    msg("chronos() error: ", conditionMessage(e))
    stop(e)
  }
)
msg("is rooted? ", is.rooted(phylo.u), "; is ultrametric? ", is.ultrametric(phylo.u))
inv.phylo <- inverseA(phylo.u, nodes="TIPS", scale=TRUE)
msg("Finished phylogeny data")

## --------- prior and model run ----------
msg("Setting prior and starting MCMC runs")
df.model <- as.data.frame(df.model)
prior <- list(
  G=list(
    G1=list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1e+3),
    G2=list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1e+3),
    G3=list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1e+3)
  ),
  R = list(V=diag(2), nu=0.002),
  S = list(mu=0, V=1e+3)
)

## Wrapper to run and capture errors
run_chain <- function(seed, outfile_prefix){
  set.seed(seed)
  tryCatch({
    ch <- MCMCglmm(
      y ~ -1 + trait:(1 + tem1.isolate.scaled + tem1.isolate.copy.number.scaled + promoter.snv),
      random = ~ phylo + isolate.id + us(at.level(trait, "mic")):phylo,
      rcov = ~ idh(trait):units,
      family = NULL,
      theta_scale = list(factor="trait", level="mic", random=1:2),
      ginverse = list(phylo = inv.phylo$Ainv),
      data = df.model,
      prior = prior,
      nitt = 40000000,
      burnin = 5000000,
      thin = 200,
      DIC = FALSE,
      pr = TRUE
    )
    # Save chain object immediately
    saveRDS(ch, file = file.path(results_dir, paste0(outfile_prefix, "_chain.rds")))
    return(ch)
  }, error = function(e){
    msg("ERROR running chain (seed=", seed, "): ", conditionMessage(e))
    writeLines(conditionMessage(e), con = file.path(results_dir, paste0(outfile_prefix, "_error.txt")))
    return(NULL)
  })
}

msg("Running chain 1")
chain.1 <- run_chain(1, "chain1")
msg("Running chain 2")
chain.2 <- run_chain(2, "chain2")

if (is.null(chain.1) || is.null(chain.2)) {
  msg("One or both chains failed. Check error files in: ", results_dir)
  quit(save = "no", status = 1)
}

msg("Model finished")

## --------- Save summaries, plots, diagnostics ----------
msg("Saving outputs to ", results_dir)

# summaries to text
summary_chain_1 <- summary(chain.1)
capture.output(summary_chain_1, file = file.path(results_dir, "chain1_summary.txt"))

summary_chain_2 <- summary(chain.2)
capture.output(summary_chain_2, file = file.path(results_dir, "chain2_summary.txt"))

# plots
pdf(file.path(results_dir, "chain1_plot.pdf"))
plot(chain.1)
dev.off()

pdf(file.path(results_dir, "chain2_plot.pdf"))
plot(chain.2)
dev.off()

# save R objects
save(chain.1, chain.2, file = file.path(results_dir, "chains.RData"))
saveRDS(df.model, file = file.path(results_dir, "df_model.rds"))

# Gelman diagnostic: create mcmc.list from Sol (posterior solutions)
mclist <- mcmc.list(chain.1$Sol, chain.2$Sol)
gelman_results <- gelman.diag(mclist)
capture.output(gelman_results, file = file.path(results_dir, "gelman_diag.txt"))

# Save gelman numeric PSRF as CSV for easier parsing
if (!is.null(gelman_results$psrf)) {
  write.csv(gelman_results$psrf, file = file.path(results_dir, "gelman_diag_psrf.csv"), row.names = TRUE)
}

msg("All outputs saved. Done.")