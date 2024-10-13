# Breast Cancer Aggregated Network Analysis Script
# This script analyzes the proportion of CIS/TRANS chromosomal interactions and 
# intra-/inter-cytoband interactions in breast cancer subtypes (Normal, Luminal A, Luminal B, Her2, Basal).
# Input: 
# - Expression matrix for each subtype (5k+ variables).
# Output:
# - Bar plots in ggplot comparing CIS/TRANS and intra-/inter-cytoband interactions across subtypes.

library(tidyr)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(igraph)

# Load gene annotation from Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
features <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "hgnc_symbol", 
              "percentage_gene_gc_content", "gene_biotype", "band")
chrs <- c(1:22, "X", "Y")
annot <- getBM(attributes = features, filters = "chromosome_name", values = chrs, mart = ensembl)
colnames(annot) <- c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Cytoband")
annot$Chr_Cytoband <- paste(annot$Chr, annot$Cytoband, sep = "")

# Function to merge and process data for chromosomal and cytoband interactions
Edges_Chromosomes_Cytobands <- function(x, y) {
  Edges <- merge(x, y, by.x = "source", by.y = "ensembl_gene_id")
  Edges <- dplyr::select(Edges, source, mi, target, Chr, Chr_Cytoband)
  colnames(Edges) <- c("source", "mi", "target", "Chromosome_Source", "Citoband_source")
  Edges <- merge(Edges, y, by.x = "target", by.y = "ensembl_gene_id")
  Edges <- dplyr::select(Edges, source, mi, target, Chromosome_Source, Chr, Citoband_source, Chr_Cytoband)
  colnames(Edges) <- c("source", "mi", "target", "Chromosome_Source", "Chromosome_Target", "Citoband_source", "Citoband_target")
  return(distinct(Edges)[order(abs(Edges$mi), decreasing = TRUE), ])
}

# Dataframes to store results
resultados_cromosomas <- data.frame()
resultados_citobandas <- data.frame()

# Subtypes to process
subtipos_procesados <- c("Normal", "Luminal A", "Luminal B", "Her2", "Basal")

# Process each subtype
for (subtipo in subtipos_procesados) {
  cat("Processing subtype:", subtipo, "\n")
  archivo_input <- file.choose()  # Choose file for each subtype
  subtipo_data <- read.table(archivo_input, header = TRUE, sep = "\t")
  subtipo_data <- dplyr::select(subtipo_data, source, mi, target)
  
  # Generate edges and calculate CIS/TRANS and cytoband interaction proportions
  subtipo_edges <- Edges_Chromosomes_Cytobands(subtipo_data, annot)
  
  cis_trans_cromosomas <- subtipo_edges %>%
    mutate(Link_Type = ifelse(Chromosome_Source == Chromosome_Target, "CIS", "TRANS")) %>%
    group_by(Link_Type) %>%
    summarise(Count = n())
  
  cis_trans_citobandas <- subtipo_edges %>%
    mutate(
      Link_Type = case_when(
        Chromosome_Source == Chromosome_Target & Citoband_source == Citoband_target ~ "Intra-cytoband",
        Chromosome_Source == Chromosome_Target & Citoband_source != Citoband_target ~ "Inter-cytoband",
        TRUE ~ "Inter-chromosome"
      )
    ) %>%
    group_by(Link_Type) %>%
    summarise(Count = n())
  
  resultados_cromosomas <- rbind(resultados_cromosomas, cis_trans_cromosomas %>% mutate(Subtipo = subtipo))
  resultados_citobandas <- rbind(resultados_citobandas, cis_trans_citobandas %>% mutate(Subtipo = subtipo))
}

# Function to create plots
crear_grafico <- function(data, title, filename) {
  data$Subtipo <- factor(data$Subtipo, levels = c("Normal", "Luminal A", "Luminal B", "Her2", "Basal"))
  
  if (all(data$Link_Type %in% c("CIS", "TRANS"))) {
    data$Link_Type <- factor(data$Link_Type, levels = c("CIS", "TRANS"))
  } else {
    data$Link_Type <- factor(data$Link_Type, levels = c("Intra-cytoband", "Inter-cytoband", "Inter-chromosome"))
  }
  
  p <- ggplot(data, aes(x = Subtipo, y = Count, fill = Link_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Breast cancer subtypes",
      y = "Number of interactions"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 2.5),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks.y = element_line(color = "black", size = 0.5),  
      axis.ticks.length = unit(0.25, "cm")
    ) +
    scale_y_continuous(
      limits = c(0, 10000),
      breaks = seq(0, 10000, by = 2500),
      expand = c(0, 0)
    ) +
    guides(fill = guide_legend(title = "Type of interaction"))
  
  ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)
}

# Create and save plots
crear_grafico(resultados_cromosomas, "CIS/TRANS Comparison Across Subtypes", "Comparacion_CIS_TRANS_Cromosomas.jpeg")
crear_grafico(resultados_citobandas, "Cytoband Interaction Comparison Across Subtypes", "Comparacion_CIS_TRANS_Citobandas.jpeg")
