################################################################################
############### NETWORK METRICS ANALYSIS ###################################

# This section computes six key network metrics: Clustering coefficient, Modularity, 
# Closeness, Degree, Global Efficiency, and Density. The analysis is performed on 
# the largest component of each network. Statistical comparisons between normal 
# and breast cancer subtypes are included, along with significance visualization.

library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(ggsignif)

# Function to compute network metrics, focusing on the largest component
calcular_metricas <- function(archivo, subtipo) {
  datos_red <- read_tsv(archivo, col_types = cols(.default = "c"))
  
  if (ncol(datos_red) < 2) {
    stop("El archivo no tiene el formato esperado.")
  }
  
  g <- graph_from_data_frame(datos_red, directed = FALSE)
  
  # Select the largest component
  comp <- components(g)
  g_largest <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  
  # Calculate only the requested metrics
  distances <- distances(g_largest)
  efficiency_global <- mean(1/distances[distances != 0])
  
  metricas <- data.frame(
    Subtipo = factor(subtipo, levels = c("Normal", "Luminal A", "Luminal B", "Her2", "Basal")),
    TCGA_ID = tools::file_path_sans_ext(basename(archivo)),
    Clustering_coeficient = mean(transitivity(g_largest, type = "average")),
    Modularity = modularity(cluster_fast_greedy(g_largest)),
    Closeness = mean(closeness(g_largest)),
    Degree = mean(degree(g_largest)),
    Global_Efficiency = efficiency_global,
    Density = edge_density(g_largest)
  )
  
  return(metricas)
}

datasets <- list()
repeat {
  subtipo <- readline(prompt = "Introduzca el subtipo (Normal, Luminal A, Luminal B, Her2, Basal) o 'salir' para terminar: ")
  if (tolower(subtipo) == "salir") { break }
  
  archivo <- file.choose()
  datasets[[length(datasets) + 1]] <- list(subtipo = subtipo, archivo = archivo)
}

# Process the selected datasets
resultados <- data.frame()
for (dataset in datasets) {
  carpeta <- dirname(dataset$archivo)
  archivos <- list.files(path = carpeta, pattern = "TCGA-.*\\.tsv$", full.names = TRUE)
  
  resultados_subtipo <- lapply(archivos, function(x) calcular_metricas(x, dataset$subtipo))
  resultados <- rbind(resultados, do.call(rbind, resultados_subtipo))
}

############################# STATISTICAL ANALYSIS ######################################
# Statistical comparisons - t-test and Wilcoxon test - Normal vs. Subtypes
metricas_a_comparar <- c("Clustering_coeficient", "Modularity", "Closeness", "Degree", 
                         "Global_Efficiency", "Density")

comparar_metricas <- function(data_normal, data_subtipo) {
  if (sd(data_normal) < 1e-10 || sd(data_subtipo) < 1e-10) {
    return(list(
      p_valor_t = NA,  # Return NA if the data is constant
      p_valor_w = NA
    ))
  }
  list(
    p_valor_t = t.test(data_normal, data_subtipo)$p.value,
    p_valor_w = wilcox.test(data_normal, data_subtipo)$p.value
  )
}

resultados_significancia <- do.call(rbind, lapply(metricas_a_comparar, function(metrica) {
  do.call(rbind, lapply(setdiff(unique(resultados$Subtipo), "Normal"), function(subtipo) {
    data_normal <- filter(resultados, Subtipo == "Normal", !is.na(!!sym(metrica)))[[metrica]]
    data_subtipo <- filter(resultados, Subtipo == subtipo, !is.na(!!sym(metrica)))[[metrica]]
    pruebas <- comparar_metricas(data_normal, data_subtipo)
    data.frame(Subtipo = subtipo, Metrica = metrica, pruebas)
  }))
}))

# Adjustment for multiple comparisons
resultados_significancia <- resultados_significancia %>%
  group_by(Metrica, Subtipo) %>%
  mutate(
    Ajuste_p_valor_t = p.adjust(p_valor_t, method = "BH"),
    Ajuste_p_valor_w = p.adjust(p_valor_w, method = "fdr"),
    Significancia_t = ifelse(Ajuste_p_valor_t < 0.001, "***", 
                             ifelse(Ajuste_p_valor_t < 0.01, "**", 
                                    ifelse(Ajuste_p_valor_t < 0.05, "*", "ns"))),
    Significancia_w = ifelse(Ajuste_p_valor_w < 0.001, "***", 
                             ifelse(Ajuste_p_valor_w < 0.01, "**", 
                                    ifelse(Ajuste_p_valor_w < 0.05, "*", "ns")))
  ) %>%
  ungroup()

# Reformat the results table for presentation
resultados_significancia_presentacion <- resultados_significancia %>%
  dplyr::select(Subtipo, Metrica, p_valor_t, p_valor_w, Ajuste_p_valor_t, Ajuste_p_valor_w, Significancia_t, Significancia_w) %>%
  rename(`p-valor t` = p_valor_t, `p-valor w` = p_valor_w, 
         `p-ajustado t` = Ajuste_p_valor_t, `p-ajustado w` = Ajuste_p_valor_w, 
         `Sig t` = Significancia_t, `Sig w` = Significancia_w)

print(resultados_significancia_presentacion)

############################# PLOTTING ######################################
# Preprocess data for plotting
df_grafico_boxplot <- resultados %>%
  pivot_longer(cols = -c(TCGA_ID, Subtipo), names_to = "Metrica", values_to = "Valor") %>%
  mutate(
    Metrica = factor(Metrica, levels = c("Clustering_coeficient", "Modularity", "Closeness", 
                                         "Degree", "Global_Efficiency", "Density")),
    Subtipo = factor(Subtipo, levels = c("Normal", "Luminal A", "Luminal B", "Her2", "Basal"))
  )

# Define colors for subtypes
subtipo_colors <- c("Normal" = "#FF0000", "Luminal A" = "#00FF00", "Luminal B" = "#0000FF", "Her2" = "#FFFF00", "Basal" = "#FF00FF")

# Significance comparisons for ggsignif
comparisons <- list(
  c("Normal", "Luminal A"),
  c("Normal", "Luminal B"),
  c("Normal", "Her2"),
  c("Normal", "Basal")
)

# Plot boxplots with jitter and significance
p <- ggplot(df_grafico_boxplot, aes(x = Subtipo, y = Valor, fill = Subtipo)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.5) +
  scale_fill_manual(values = subtipo_colors) +
  facet_wrap(~ Metrica, scales = "free_y") +
  theme_bw() +
  labs(title = "",
       x = NULL, y = "Value") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.text.x = element_blank(),    # Remove x-axis labels
    axis.ticks.x = element_blank(),   # Remove x-axis ticks
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel borders
    axis.line.x = element_line(col
