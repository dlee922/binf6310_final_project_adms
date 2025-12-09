###
### Import libraries and data
###

library(readr)
library(vcfR)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ordinal)
library(lme4)
library(clevr)
library(cowplot)
library(ggtext)
library(tidyverse)

# Main metadata file
df <- read_delim("data/supplementary_data_1.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)

# blaTEM-1 and linked promoter annotations
tem1.report <- read_delim("data/supplementary_data_2.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

###
### A curated sample of E. coli isolates with hybrid assemblies and co-amoxiclav MICs
###

# Sample size
length(unique(as.factor(df$isolate.id)))
# Sampling date range
summary(as.Date(df[df$contig.type=="chromosome",]$date, format = "%d/%m/%Y"))
# Complete assemblies
summary(as.factor(df[df$contig.type=="chromosome",]$is.complete))
# Plasmid stats
df %>% 
  filter(contig.type != "chromosome") %>%
  group_by(isolate.id) %>%
  dplyr::summarise(n.contigs=n()) %>%
  dplyr::summarise(median.plasmids = median(n.contigs),
                   lq.plasmids = quantile(n.contigs, 0.25),
                   uq.plasmids = quantile(n.contigs, 0.75))

# blaTEM-1 isolate copies
summary(df[df$tem1.isolate>0,]$tem1.isolate)
# blaTEM-1-positive contigs
nrow(df[df$tem1.replicon>0,])
nrow(df[df$tem1.replicon>1,])
# blaTEM-1 contig types
summary(as.factor(df[df$tem1.replicon>0,]$contig.type))
# blaTEM-1 SNV profiles
summary(as.factor(tem1.report$tem1.snv))
# blaTEM-1 repeated SNV profiles
tem1.report %>%
  group_by(isolate.assembly, contig.assembly) %>%
  dplyr::summarise(n.tem1 = n(),
            n.tem1.snv = n_distinct(tem1.snv)) %>%
  filter(n.tem1 > 1)

# Sequence length vs blaTEM-1 carriage
sum(df[df$contig.type=="plasmid",]$contig.length)
sum(df[df$contig.type=="chromosome",]$contig.length)
sum(df[df$contig.type=="plasmid",]$tem1.replicon)
sum(df[df$contig.type=="chromosome",]$tem1.replicon)

# Phylogroup stats
summary(as.factor(df[df$contig.type=="chromosome",]$ezclermont.phylogroup))

# Plasmid stats
nrow(df[df$contig.type=="plasmid",])
nrow(df[df$contig.type=="plasmid" & df$contig.topology=="circular",])
summary(as.factor(df[df$contig.type=="plasmid" & df$contig.topology=="circular",]$plasmidfinder))

# Promoter stats
nrow(tem1.report[!is.na(tem1.report$promoter.snv),])
summary(as.factor(tem1.report$promoter.snv))

###
### Figure S3 ###
###

tem1.promoter.df <- tem1.report %>%
  select(isolate.assembly, promoter.snv) %>%
  na.omit

df.figS3 <- merge(tem1.promoter.df, df, by="isolate.assembly", all.x=TRUE) %>%
  filter(contig.type=="chromosome") %>%
  select(ezclermont.phylogroup, promoter.snv)

color_palette <- brewer.pal(n = 10, name = "Set3") 


plot_S3a <- ggplot(df.figS3, aes(x = ezclermont.phylogroup , fill = promoter.snv)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = color_palette) +
  labs(y = "Combination", x = "", fill = "blaTEM-1 promoter SNV") +
  labs(title = "Raw counts", x = "Phylogroup", y = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom")

df.figS3 <- df.figS3 %>%
  count(ezclermont.phylogroup, promoter.snv) %>%
  group_by(ezclermont.phylogroup) %>%
  mutate(prop = n / sum(n))

plot_S3b <- ggplot(df.figS3, aes(x = ezclermont.phylogroup, y = prop, fill = promoter.snv)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = color_palette) +
  labs(title = "Proportions", x = "Phylogroup", y = "Proportion", fill = "blaTEM-1 promoter SNV") +
  theme_minimal() +
  theme(legend.position = "none")

plot_grid(plot_S3a, plot_S3b, labels = c("a", "b"), ncol = 2, align = "h")

# MIC stats
summary(as.factor(df[df$contig.type=="chromosome",]$coamox.mic))

###
### Figure S4 ###
###

df.figS4 <- df %>%
  filter(contig.type=="chromosome") %>%
  select(ezclermont.phylogroup, coamox.mic)

df.figS4$coamox.mic <- factor(df.figS4$coamox.mic, levels = c("<=2.2", "4.2", "8.2",
                                                            "16.2", "32.2", ">32.2"))

color_palette <- brewer.pal(n = 6, name = "PuBu") 

plot_S4a <- ggplot(df.figS4, aes(x = ezclermont.phylogroup , fill = coamox.mic)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = color_palette) +
  labs(y = "Combination", x = "", fill = "Co-amoxiclav MIC") +
  labs(title = "Raw counts", x = "Phylogroup", y = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom")

df.figS4 <- df.figS4 %>%
  count(ezclermont.phylogroup, coamox.mic) %>%
  group_by(ezclermont.phylogroup) %>%
  mutate(prop = n / sum(n))

plot_S4b <- ggplot(df.figS4, aes(x = ezclermont.phylogroup, y = prop, fill = coamox.mic)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = color_palette) +
  labs(title = "Proportions", x = "Phylogroup", y = "Proportion", fill = "Co-amoxiclav MIC") +
  theme_minimal() +
  theme(legend.position = "none")

plot_grid(plot_S4a, plot_S4b, labels = c("a", "b"), ncol = 2, align = "h")

###
### Plot Figure 1
###

add_na_to_list <- function(lst) {
  if (length(lst) == 0 || !grepl("^chromosome-", lst[1])) {
    c(NA, lst)
  } else {
    lst
  }
}

df.fig1 <- df %>%
  filter(tem1.replicon>0) %>%
  select(isolate.id, contig.type, tem1.replicon) %>%
  group_by(isolate.id, contig.type) %>%
  mutate(contig.copies = paste(contig.type, tem1.replicon, sep="-")) %>%
  arrange(desc(tem1.replicon), contig.type) %>%
  group_by(isolate.id) %>%
  mutate(distribution = list(contig.copies)) %>%
  distinct(isolate.id, distribution) %>%
  mutate(distribution = lapply(distribution, add_na_to_list)) %>%
  tidyr::unnest_wider(distribution, names_sep = "_") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate_all(~str_replace_all(., "chromosome-|plasmid-", "")) %>%
  rename(chromosome = distribution_1,  plasmid.1 = distribution_2,
           plasmid.2 = distribution_3, plasmid.3 = distribution_4)

df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  dplyr::summarise(count = n()) %>%
  arrange(desc(count))

df.fig1 <- merge(df.fig1, df[c("isolate.id", "ezclermont.phylogroup", "coamox.mic")], by = "isolate.id", all.x= TRUE)
df.fig1 <- unique(df.fig1[!is.na(df.fig1$ezclermont.phylogroup),])
    
df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

df.fig1.summary <- df.fig1.summary[,-c(5)]
df.fig1.summary$arrangement.id <- nrow(df.fig1.summary):1
df.fig1.summary.melt <- melt(df.fig1.summary, id.vars = "arrangement.id")

# Plot TEM-1 arrangements for Figure 1a

p.fig1.a <- ggplot(df.fig1.summary.melt, aes(x = variable, y = arrangement.id, fill = value)) +
  geom_tile(color = "black",  width = 0.9, height = 0.9, size=0.5) +
  scale_fill_manual(values = c("white", rep("black", 4))) +
  geom_text(aes(label = ifelse(value != 0, value, "")), color = "white") +
  theme_minimal() +
  theme(
    legend.position = "none",  
    panel.grid = element_blank(), 
    axis.title.y = element_blank(),  
    axis.title.x = element_blank(), 
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1.5, color = "black")) +
  scale_x_discrete(labels = c("chromosome" = "Chromosome", 
                              "plasmid.1" = "Plasmid 1", "plasmid.2" = "Plasmid 2", "plasmid.3" = "Plasmid 3")) +  # Rename x-axis labels
  coord_fixed(ratio = 1, ylim = c(0.5, 12.5)) 

df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

df.fig1.summary$arrangement.id <- nrow(df.fig1.summary):1
df.fig1.summary <- df.fig1.summary[c("count", "arrangement.id")]
df.fig1.summary$count <- as.numeric(df.fig1.summary$count)

# Plot TEM-1 arrangement counts for Figure 1b

p.fig1.b <- ggplot(df.fig1.summary, aes(x=1, y = arrangement.id, fill = count)) +
  geom_tile(color = "black",  width = 0.9, height = 0.9, size=0.5) +
  geom_text(aes(label = ifelse(count != 0, count, "")), color = "Black") +
  theme_minimal() +
  scale_fill_distiller(palette = "OrRd", direction=1) +
  theme(
    legend.position = "none",  
    panel.grid = element_blank(), 
    axis.title.y = element_blank(),  
    axis.title.x = element_blank(), 
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank()) +
  coord_fixed(ratio = 1, ylim = c(0.5, 12.5))

df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

df.fig1.summary <- df.fig1.summary[,-c(5)]
df.fig1.summary$arrangement.id <- nrow(df.fig1.summary):1

df.fig1[, c("chromosome", "plasmid.1", "plasmid.2", "plasmid.3")] <- lapply(df.fig1[, c("chromosome", "plasmid.1", "plasmid.2", "plasmid.3")], as.character)
df.fig1 <- merge(df.fig1, df.fig1.summary, by=c("chromosome", "plasmid.1", "plasmid.2", "plasmid.3"), all.x=TRUE)

# Plot TEM-1 arrangement group phylogroups for Figure 1c

color_palette <- brewer.pal(n = 9, name = "Set3")

p.fig1.c <- ggplot(df.fig1, aes(x = reorder(interaction(chromosome, plasmid.1, plasmid.2, plasmid.3), arrangement.id), fill = ezclermont.phylogroup)) +
  geom_bar(position = "fill", show.legend = FALSE) +
  geom_bar(position = "fill", color = "black", size = 0.5) +
  scale_fill_manual(values = color_palette) +
  labs(y = "Combination", x = "", fill = "Phylogroup") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_blank(),  # Remove background color
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1.5, color = "black")
  ) +
  coord_flip() +
  theme(legend.position = "none")

# Plot TEM-1 arrangement group co-amoxicillin MIC for Figure 1c

df.fig1$coamox.mic <- factor(df.fig1$coamox.mic, levels = c("<=2.2", "4.2", "8.2", ">8.2",
                                                            "16.2", "32.2", ">32.2"))

color_palette <- brewer.pal(n = 7, name = "PuBu") 

p.fig1.d <- ggplot(df.fig1, aes(x = reorder(interaction(chromosome, plasmid.1, plasmid.2, plasmid.3), arrangement.id), fill = coamox.mic)) +
  geom_bar(position = "fill", show.legend = FALSE) +
  geom_bar(position = "fill", color = "black", size = 0.5) +
  scale_fill_manual(values = color_palette) +
  labs(y = "Combination", x = "", fill = "Co-amoxiclav MIC") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_blank(),  # Remove background color
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1.5, color = "black")
  ) +
  coord_flip() +
  theme(legend.position = "none")

###
### blaTEM-1 associated with conjugative plasmids.
###

df.plasmid <- df[df$contig.type=="plasmid" & df$contig.topology=="circular",]

df.plasmid$predicted.mobility <- factor(df.plasmid$predicted.mobility, levels = c("non-mobilizable",
                                                                                  "mobilizable", "conjugative"))
df.plasmid$carries.tem1 <- ifelse(df.plasmid$tem1.replicon>0, 1, 0)
summary(as.factor(df.plasmid$carries.tem1))

tem1.counts <- c(200, 286)
totals <- c(117 + 386 + 200, 20 + 27 + 286)
prop.test.result <- prop.test(tem1.counts, totals)
print(prop.test.result)

df.plasmid.2 <- df[df$contig.type=="plasmid",]
summary(df.plasmid.2$contig.copy.number)
nrow(df.plasmid.2[df.plasmid.2$contig.copy.number<1,])

df.plasmid.3 <- df.plasmid[df.plasmid$contig.copy.number>=1,]
summary(df.plasmid.3[df.plasmid.3$contig.length>10000,]$contig.copy.number)
summary(df.plasmid.3[df.plasmid.3$contig.length<=10000,]$contig.copy.number)


###
### E. coli phylogeny structures ampC diversity.
###

df.ampc <- df %>%
  filter(contig.type=="chromosome")

homogeneity(df.ampc$ampc.snv, df.ampc$ezclermont.phylogroup)
completeness(df.ampc$ampc.snv, df.ampc$ezclermont.phylogroup)

sum(df.ampc$mlst.st=="-")

df.ampc <- df.ampc %>%
  filter(mlst.st != "-")

homogeneity(df.ampc$ampc.snv, df.ampc$mlst.st)
completeness(df.ampc$ampc.snv, df.ampc$mlst.st)

# promoter info
summary(as.factor(df[df$contig.type=="chromosome",]$ampc.promoter.snv))

promoter.mutations <- df %>%
  filter(contig.type == "chromosome") %>%
  select(isolate.id, ampc.promoter.snv) %>%
  mutate(position = row_number()) %>%
  separate(ampc.promoter.snv, into = paste0("position_", 1:12), sep = "") %>%
  select(isolate.id, paste0("position_", 2:12)) %>%
  mutate_all(as.factor)
summary(promoter.mutations)

homogeneity(df.ampc$ampc.snv, df.ampc$ampc.promoter.snv)
completeness(df.ampc$ampc.snv, df.ampc$ampc.promoter.snv)

###
### Figure S5
###

library(ape)
library(phytools)
library(readxl)
library(ggtree)
library(ggheatmap)

phylo <- read.tree("data/ml_tree.txt")

exp_df <- read_excel("data/supplementary_data_3.xlsx")

df.plot <- df %>%
  mutate(tem1.replicon.copy.number = tem1.replicon * contig.copy.number,
         ) %>%
  group_by(isolate.assembly) %>%
  mutate(tem1.cell.copy.number = sum(tem1.replicon.copy.number)) %>%
  filter(contig.type == "chromosome") %>%
  select(isolate.assembly, isolate.id, ezclermont.phylogroup, coamox.mic, tem1.isolate, tem1.cell.copy.number, accession)

exp.plot <- exp_df
exp.plot$Isolate <- sapply(strsplit(as.character(exp.plot$Isolate), " "), function(x) x[1])
exp.plot <- exp.plot[c("Isolate", "delta Ct (control)")]
colnames(exp.plot) <- c("accession", "exp")

df.plot <- merge(df.plot, exp.plot, by="accession", all.x=TRUE) %>%
  select(isolate.assembly, isolate.id, ezclermont.phylogroup, tem1.isolate, tem1.cell.copy.number, exp, coamox.mic) %>%
  as.data.frame() 

phylo <- midpoint.root(phylo)
phylo$tip.label <- as.character(phylo$tip.label)

p <- ggtree(phylo)

p1 <- p %<+% df.plot + 
  geom_tippoint(aes(fill = ezclermont.phylogroup), shape = 21, size = 2, stroke = 0.5) +
  scale_fill_brewer(palette = "Set3", name="Phylogroup") + theme_tree2()

d2 <- df.plot %>%
  select("label" = isolate.assembly, "value" = tem1.isolate) %>%
  unique()

p2 <- facet_plot(p1, panel="Genome copies", data=d2, geom=geom_point, mapping = aes(x=value), color="firebrick") + theme_tree2()

d3 <- df.plot %>%
  select("label" = isolate.assembly, "value" = tem1.cell.copy.number)  %>%
  unique()

p3 <- facet_plot(p2, panel="Cell copies", data=d3, geom=geom_point, mapping = aes(x=value), color="firebrick") + theme_tree2()

d4 <- df.plot %>%
  select("label" = isolate.assembly, "value" = exp) %>%
  group_by(label) %>%
  summarise(value = mean(value))

p4 <- facet_plot(p3, panel="Mean delta Ct", data=d4, geom=geom_point, mapping = aes(x=value), color="firebrick") + 
  theme_tree2() + theme(legend.position='bottom') 

p4

d5 <- df.plot %>%
  select(label = isolate.assembly, value = coamox.mic) %>%
  mutate(value = factor(value, levels = c("<=2.2", "4.2", "8.2", "16.2", "32.2", ">32.2"))) %>%
  unique()

tip.order <- phylo$tip.label

p5 <- ggplot(d5, aes(x=value, y=factor(label, levels=tip.order))) +
  geom_point(color="firebrick") +
  scale_x_discrete(limits=levels(d5$value), expand=expansion(add=0.5)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(
    axis.title.y   = element_blank(),
    axis.text.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +   theme_tree2()


library(patchwork)

p_final <- p4 + p5 +
  plot_layout(ncol=2, widths=c(3, 0.7))
p_final

###
### Plot Figure S6
###

plasmid_info <- df %>%
  filter(contig.topology == "circular", contig.type == "plasmid") %>%
  select(isolate.id, contig.id, plasmidfinder)

chromosome_info <- df %>%
  filter(contig.type == "chromosome") %>%
  select(isolate.id, mlst.st, ezclermont.phylogroup)

df.figS6 <- merge(plasmid_info, chromosome_info, by="isolate.id", all.x=TRUE)

heatmap_phylo <- df.figS6 %>%
  count(ezclermont.phylogroup, plasmidfinder) %>%
  ggplot(aes(x = plasmidfinder, y = ezclermont.phylogroup, fill = n)) +
  geom_tile(color = "black") + 
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(x = "Plasmid replicon", y = "Phylogroup", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))

heatmap_st <- df.figS6 %>%
  count(mlst.st, plasmidfinder) %>%
  ggplot(aes(x = plasmidfinder, y = factor(mlst.st), fill = n)) +
  geom_tile(color = "black") +  
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = NULL, y = "ST", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
phylogroup.heatmap <- df.figS6 %>%
  count(ezclermont.phylogroup, plasmidfinder)

(heatmap_st / heatmap_phylo) + plot_layout(heights = c(7.75, 1))

###
### Plot Figure S7
###


p.fig.s7 <- ggplot(df.plasmid.3, aes(x = contig.copy.number, y = contig.length)) +
  geom_density_2d(h=NULL, colour='#ef8a62', size=1) +
  geom_point(alpha=0.5, size=2, stroke=0) +
  labs(x = "Plasmid copy number (log<sub>10</sub>)",
       y = "Plasmid length (bp, log<sub>10</sub>)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=20, color='black'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(color="white"), 
        axis.ticks = element_line(size = 1, color='black'),
        axis.text.x = element_text(color='black'),
        axis.text.y = element_text(color='black'),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown()) +
  scale_y_continuous(
    trans = "log10",
    breaks = 10^seq(4,5),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.1, 1, 10),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))
