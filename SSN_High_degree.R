# This script processes unique networks of breast cancer subtypes to identify
# the top 10 highest-degree genes from the largest component of the network.
# The results are joined with cytoband information and HGNC symbols and 
# visualized as bar charts.
#
# Input: 
# - Network files (.tsv) for each subtype.
# - Gene annotation from Ensembl (biomaRt).
#
# Output: 
# - TSV files with top 10 highest-degree genes for each network.
# - Bar plots showing the frequency of top genes and cytobands.
#
################################################################################

# Load required libraries
library(igraph)
library(dplyr)
library(readr)
library(biomaRt)
library(ggplot2)

# Load gene annotation from Ensembl (biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
annot <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", 
                 "end_position", "hgnc_symbol", "gene_biotype", "band"), 
  filters = "chromosome_name", values = c(1:22, "X", "Y"), 
  mart = ensembl
) %>% 
  filter(gene_biotype == "protein_coding") %>%
  mutate(Chr_Cytoband = paste(chromosome_name, band, sep=""))

# Function to process each network file and extract the top 10 highest-degree genes
procesar_red <- function(archivo, ruta_directorio) {
  g <- graph_from_data_frame(read_tsv(archivo, col_types = cols(.default = "c")), directed = FALSE)
  
  # Select the largest component of the network
  g_largest <- induced_subgraph(g, which(components(g)$membership == which.max(components(g)$csize)))
  
  # Extract the top 10 highest-degree genes
  top_genes <- sort(degree(g_largest), decreasing = TRUE)[1:10]
  
  # Join with cytoband and HGNC symbol information
  resultado <- data.frame(Gene = names(top_genes), Degree = top_genes) %>%
    left_join(annot, by = c("Gene" = "ensembl_gene_id")) %>%
    dplyr::select(Gene, Degree, hgnc_symbol, Chr_Cytoband)
  
  # Create subdirectory and save the results
  nombre_directorio <- file.path(ruta_directorio, tools::file_path_sans_ext(basename(archivo)))
  dir.create(nombre_directorio, showWarnings = FALSE)
  write_tsv(resultado, file.path(nombre_directorio, paste0("top_genes_", basename(archivo))))
}

# Prompt the user to select the main directory
ruta_directorio <- dirname(file.choose())

# List the network files in the selected directory
archivos_red <- list.files(path = ruta_directorio, pattern = "TCGA-.*\\.tsv$", full.names = TRUE)

# Apply the function to all network files
lapply(archivos_red, procesar_red, ruta_directorio = ruta_directorio)

############################# PLOTTING FUNCTIONS ###############################
# Function to create bar plots and save them as JPEG
graficar_barras <- function(data, x_var, y_var, x_lab, y_lab, title, subtipo) {
  plot <- ggplot(data, aes(x = reorder(!!sym(x_var), -!!sym(y_var)), y = !!sym(y_var))) +
    geom_bar(stat = "identity", fill = "#0073C2FF", color = "black", width = 0.7) +  # Improve colors and borders
    geom_text(aes(label = !!sym(y_var)), vjust = -0.3, size = 4) +  # Add value labels on bars
    labs(title = paste(title, subtipo), x = x_lab, y = y_lab) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Prominent title
      axis.title = element_text(size = 14, face = "bold"),  # Strong axis titles
      axis.text = element_text(size = 12),  # Larger axis text for readability
      panel.grid.major = element_line(color = "grey90"),  # Lighter grid lines
      panel.grid.minor = element_blank(),  # No minor grid lines for clarity
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotated x-axis text
      axis.ticks = element_line(size = 0.5),  # Clear axis ticks
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Border around the plot
    )
  
  # Save the plot as a JPEG file
  ggsave(filename = paste0("grafico_", subtipo, "_", x_lab, ".jpeg"), 
         plot = plot, 
         dpi = 300, 
         width = 10, 
         height = 7, 
         device = "jpeg")
}

######################### FREQUENCY OF HIGH-DEGREE GENES ########################
# Function to count gene frequencies across multiple network files
contar_frecuencias <- function(directorio) {
  archivos <- list.files(path = directorio, pattern = "top_genes_.*\\.tsv$", full.names = TRUE)
  
  # Read and combine the genes from all files
  genes <- lapply(archivos, function(archivo) {
    read_tsv(archivo, col_types = cols()) %>% pull(Gene)
  }) %>% unlist()
  
  # Count the frequency of each gene
  tabla_frecuencias <- as.data.frame(table(genes))
  colnames(tabla_frecuencias) <- c("Gen", "Frecuencia")
  return(tabla_frecuencias)
}

# Prompt the user to select the main directory containing the TCGA files
directorio_principal <- dirname(file.choose())
subdirectorios <- list.dirs(path = directorio_principal, full.names = TRUE, recursive = FALSE)
subdirectorios_tcga <- subdirectorios[grepl("^TCGA-", basename(subdirectorios))]

# Process all subdirectories and count gene frequencies
frecuencias_genes <- lapply(subdirectorios_tcga, contar_frecuencias)
frecuencias_genes_unidas <- bind_rows(frecuencias_genes)

# Sum the frequencies of each gene across all directories
frecuencias_totales <- frecuencias_genes_unidas %>%
  group_by(Gen) %>%
  summarise(Frecuencia = sum(Frecuencia)) %>%
  arrange(desc(Frecuencia))

########################### TOP 10 HIGH-DEGREE GENES ############################
# Top 10 highest-frequency genes joined with HGNC symbols
frecuencias_top10 <- frecuencias_totales %>%
  left_join(annot, by = c("Gen" = "ensembl_gene_id")) %>%
  filter(!is.na(hgnc_symbol)) %>%
  group_by(hgnc_symbol) %>%
  summarise(Frecuencia = sum(Frecuencia)) %>%
  arrange(desc(Frecuencia)) %>%
  top_n(10, Frecuencia)

# Prompt the user to input the subtype for the plot
subtipo <- readline(prompt = "Ingresa el subtipo del fenotipo para el gráfico: ")

# Plot the top 10 high-degree genes
graficar_barras(frecuencias_top10, "hgnc_symbol", "Frecuencia", "Genes", "Frequency", "Top 10 Genes by Degree", subtipo)

########################## TOP 10 CYTOBAND FREQUENCIES ##########################
# Top 10 cytobands where the highest-degree genes are located
frecuencias_top10_citobanda <- frecuencias_totales %>%
  left_join(annot, by = c("Gen" = "ensembl_gene_id")) %>%
  filter(!is.na(hgnc_symbol)) %>%
  group_by(Chr_Cytoband) %>%
  summarise(Frecuencia = sum(Frecuencia)) %>%
  arrange(desc(Frecuencia)) %>%
  top_n(10, Frecuencia)

# Plot the top 10 cytobands with the highest gene frequencies
graficar_barras(frecuencias_top10_citobanda, "Chr_Cytoband", "Frecuencia", "Cytoband", "Frequency", "Top 10 Cytobands by Gene Frequency", subtipo)

################################################################################
