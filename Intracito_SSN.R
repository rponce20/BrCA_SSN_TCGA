# Single-Sample Network Analysis for CIS/TRANS and Intra-/Inter-Cytoband Interactions
# This script analyzes single-sample coexpression networks across breast cancer subtypes (Normal, Luminal A, Luminal B, Her2, Basal).
# The analysis focuses on the proportion of CIS/TRANS chromosomal interactions and intra-/inter-cytoband interactions.
# Input:
# - Expression matrices for each subtype.
# Output:
# - Bar plots and statistical significance tests for the proportion of CIS/TRANS and intra-/inter-cytoband interactions.

library(tidyr)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(readr)
library(ggsignif)

# Load gene annotation from Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
annot <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", 
                              "end_position", "hgnc_symbol", "gene_biotype", "band"), 
               filters = "chromosome_name", values = c(1:22, "X", "Y"), mart = ensembl) %>% 
  filter(gene_biotype == "protein_coding") %>%
  mutate(Chr_Cytoband = paste(chromosome_name, band, sep=""))

# Function to merge data and obtain chromosome and cytoband interactions
Edges_Chromosomes_Cytobands <- function(x, y) {
  x %>% inner_join(y, by = c("source" = "ensembl_gene_id")) %>% 
    rename(Chromosome_Source = chromosome_name, Citoband_source = Chr_Cytoband) %>%
    inner_join(y, by = c("target" = "ensembl_gene_id")) %>% 
    rename(Chromosome_Target = chromosome_name, Citoband_target = Chr_Cytoband) %>%
    dplyr::select(source, target, Chromosome_Source, Chromosome_Target, Citoband_source, Citoband_target) %>% 
    distinct()
}

# Select directories for each subtype
subtipos_directorios <- setNames(lapply(c("Normal", "Luminal A", "Luminal B", "Her2", "Basal"), function(subtipo) {
  cat("Select the directory for breast cancer subtype:", subtipo, "\n")
  dirname(file.choose())
}), c("Normal", "Luminal A", "Luminal B", "Her2", "Basal"))

# Process files and retrieve results for each subtype
procesar_archivos <- function(directorio_redes, subtipo) {
  archivos_redes <- list.files(path = directorio_redes, pattern = "TCGA-.*\\.tsv$", full.names = TRUE)
  resultados <- lapply(archivos_redes, function(archivo) {
    cat("Processing file:", archivo, "for subtype:", subtipo, "\n")
    red_data <- read_tsv(archivo, col_types = cols(.default = "c"))
    red_edges <- Edges_Chromosomes_Cytobands(red_data, annot)
    
    list(
      cromosoma = red_edges %>% 
        mutate(Chromosome_Link_Type = ifelse(Chromosome_Source == Chromosome_Target, "CIS", "TRANS")) %>% 
        count(Chromosome_Link_Type) %>% 
        mutate(Subtipo = subtipo, Archivo = basename(archivo)),
      citobanda = red_edges %>%
        mutate(Cytoband_Link_Type = case_when(
          Chromosome_Source == Chromosome_Target & Citoband_source == Citoband_target ~ "Intracitoband",
          Chromosome_Source == Chromosome_Target & Citoband_source != Citoband_target ~ "Intercitoband",
          TRUE ~ "Interchromosomal")) %>%
        count(Cytoband_Link_Type) %>% 
        mutate(Subtipo = subtipo, Archivo = basename(archivo))
    )
  })
  list(cromosoma = bind_rows(lapply(resultados, `[[`, "cromosoma")),
       citobanda = bind_rows(lapply(resultados, `[[`, "citobanda")))
}

resultados <- lapply(names(subtipos_directorios), function(subtipo) {
  procesar_archivos(subtipos_directorios[[subtipo]], subtipo)
})

# Combine results into final dataframes
resultados_cromosoma_final <- bind_rows(lapply(resultados, `[[`, "cromosoma"))
resultados_citobanda_final <- bind_rows(lapply(resultados, `[[`, "citobanda"))

# Order results by subtype
niveles_subtipos <- c("Normal", "Luminal A", "Luminal B", "Her2", "Basal")
resultados_cromosoma_final$Subtipo <- factor(resultados_cromosoma_final$Subtipo, levels = niveles_subtipos)
resultados_citobanda_final$Subtipo <- factor(resultados_citobanda_final$Subtipo, levels = niveles_subtipos)

############################# STATISTICAL ANALYSIS ######################################
# Filter healthy and cancerous samples
datos_sanos <- resultados_cromosoma_final %>% filter(Subtipo == "Normal")
datos_enfermos <- resultados_cromosoma_final %>% filter(Subtipo != "Normal")

# Function to perform statistical tests and adjust p-values
realizar_pruebas <- function(test_func, data_sanos, data_enfermos) {
  subtipos <- unique(data_enfermos$Subtipo)
  tipos_enlace <- c("CIS", "TRANS")
  
  resultados_pruebas <- lapply(subtipos, function(subtipo) {
    lapply(tipos_enlace, function(tipo) {
      datos_subtipo <- filter(data_enfermos, Subtipo == subtipo, Chromosome_Link_Type == tipo)
      datos_normal <- filter(data_sanos, Chromosome_Link_Type == tipo)
      
      prueba <- test_func(datos_subtipo$n, datos_normal$n)
      return(prueba$p.value)
    })
  })
  
  p_valores <- unlist(resultados_pruebas)
  p_ajustados <- p.adjust(p_valores, method = "BH")
  
  return(list(resultados_pruebas = resultados_pruebas, p_ajustados = p_ajustados))
}

# Perform t-test and Mann-Whitney U test
resultados_t_test <- realizar_pruebas(t.test, datos_sanos, datos_enfermos)
resultados_mann_whitney <- realizar_pruebas(wilcox.test, datos_sanos, datos_enfermos)

# Generate significance tables
generar_tabla_significancia <- function(resultados_pruebas, subtipos, tipos_enlace) {
  Comparación <- expand.grid(subtipos, tipos_enlace) %>%
    apply(1, function(x) paste("Normal vs", x[1], x[2])) %>%
    unlist()
  
  P_valor <- unlist(resultados_pruebas$resultados_pruebas)
  
  tabla_significancia <- data.frame(
    Comparación = Comparación,
    P_valor = P_valor,
    P_ajustado = resultados_pruebas$p_ajustados,
    Significancia = ifelse(resultados_pruebas$p_ajustados < 0.001, "***",
                           ifelse(resultados_pruebas$p_ajustados < 0.01, "**",
                                  ifelse(resultados_pruebas$p_ajustados < 0.05, "*", "ns")))
  )
  return(tabla_significancia)
}

subtipos <- unique(datos_enfermos$Subtipo)
tipos_enlace <- c("CIS", "TRANS")

tabla_significancia_t_test <- generar_tabla_significancia(resultados_t_test, subtipos, tipos_enlace)
tabla_significancia_mann_whitney <- generar_tabla_significancia(resultados_mann_whitney, subtipos, tipos_enlace)

# Display significance tables
print(tabla_significancia_t_test)
print(tabla_significancia_mann_whitney)

###################### SIGNIFICANCE PLOTS ################################ 
# Ensure no NA values in 'resultados_cromosoma_final'
resultados_cromosoma_final <- resultados_cromosoma_final %>%
  filter(!is.na(Subtipo) & !is.na(Chromosome_Link_Type) & !is.na(n))

# Adjust factor levels for consistency
resultados_cromosoma_final$Subtipo <- factor(resultados_cromosoma_final$Subtipo, 
                                             levels = c("Normal", "Luminal A", "Luminal B", "Her2", "Basal"))
resultados_cromosoma_final$Chromosome_Link_Type <- factor(resultados_cromosoma_final$Chromosome_Link_Type, 
                                                          levels = c("CIS", "TRANS"))

# Ensure no NA values in 'resultados_citobanda_final'
resultados_citobanda_final <- resultados_citobanda_final %>%
  filter(!is.na(Subtipo) & !is.na(Cytoband_Link_Type) & !is.na(n))

# Adjust factor levels for consistency
resultados_citobanda_final$Subtipo <- factor(resultados_citobanda_final$Subtipo, 
                                             levels = c("Normal", "Luminal A", "Luminal B", "Her2", "Basal"))
resultados_citobanda_final$Cytoband_Link_Type <- factor(resultados_citobanda_final$Cytoband_Link_Type, 
                                                        levels = c("Intracitoband", "Intercitoband", "Interchromosomal"))

# Boxplot with jitter for chromosomal interactions with significance
ggplot(resultados_cromosoma_final, aes(x = Subtipo, y = n, fill = Subtipo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), size = 0.5, alpha = 0.5) +
  facet_wrap(~ Chromosome_Link_Type, scales = "free") +
  labs(x = "Subtypes", y = "Count of edges") +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank()) +  
  geom_signif(comparisons = list(c("Normal", "Luminal A"), c("Normal", "Luminal B"), 
                                 c("Normal", "Her2"), c("Normal", "Basal")),
              map_signif_level = TRUE, test = "wilcox.test", 
              y_position = c(9500, 10000, 10500, 11000),
              step_increase = 0.1)

ggsave("chromosome_links_single_samples_jitter_signif.jpeg", width = 10, height = 8, dpi = 300)

# Boxplot with jitter for cytoband interactions with significance
ggplot(resultados_citobanda_final, aes(x = Subtipo, y = n, fill = Subtipo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), size = 0.5, alpha = 0.5) +
  facet_wrap(~ Cytoband_Link_Type, scales = "free") +
  labs(x = "Subtypes", y = "Count of edges") +
  scale_x_discrete(labels = c("Intracitoband" = "Intra-cytoband", 
                              "Intercitoband" = "Inter-cytoband", 
                              "Interchromosomal" = "Inter-chromosome")) +
  scale_fill_discrete(name = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank()) +  
  geom_signif(comparisons = list(c("Normal", "Luminal A"), c("Normal", "Luminal B"), 
                                 c("Normal", "Her2"), c("Normal", "Basal")),
              map_signif_level = TRUE, test = "wilcox.test", 
              y_position = c(9000, 9500, 10000, 10500),
              step_increase = 0.1)

ggsave("cytoband_links_single_samples_jitter_signif.jpeg", width = 10, height = 8, dpi = 300)
