# This script calculates the proportion of CIS links in breast cancer subtypes
# and performs survival analysis to determine the association between high/low 
# CIS proportions and patient survival over 5 years.
#
# Input: 
# - Aggregated network files (.tsv) for each subtype (Normal, Luminal A, 
#   Luminal B, Her2, Basal).
# - Clinical data for survival analysis.
#
# Output: 
# - CSV files of high/low CIS proportion groups for each subtype.
# - Kaplan-Meier survival plots.
#
################################################################################

# Load required libraries
library(dplyr)
library(biomaRt)
library(readr)
library(data.table)

# Load gene annotation data from Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
annot <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "gene_biotype"),
  filters = "chromosome_name", 
  values = c(1:22, "X", "Y"), 
  mart = ensembl
) %>%
  dplyr::select(ensembl_gene_id, Chr = chromosome_name)

# Function to calculate CIS proportion for each network file
calc_cis_proportion <- function(file, annot) {
  red_data <- fread(file, col.names = c("source", "target", "mi"))
  red_edges <- merge(red_data, annot, by.x = "source", by.y = "ensembl_gene_id") %>%
    merge(annot, by.x = "target", by.y = "ensembl_gene_id", suffixes = c("_source", "_target"))
  
  total_enlaces <- nrow(red_edges)
  total_cis <- sum(red_edges$Chr_source == red_edges$Chr_target)
  
  data.frame(
    TCGA_ID = sub("\\.tsv$", "", basename(file)),
    Total_Enlaces = total_enlaces,
    Total_CIS = total_cis,
    Proporcion_CIS = total_cis / total_enlaces
  )
}

# Define breast cancer subtypes
subtipos <- c("Normal", "Luminal A", "Luminal B", "Her2", "Basal")

# Select directories for each subtype
subtipos_directorios <- lapply(subtipos, function(subtipo) {
  cat("Seleccione el directorio para el subtipo:", subtipo, "\n")
  return(dirname(file.choose()))
})

# Process all files and calculate CIS proportions
resultados_proporciones_cis <- lapply(seq_along(subtipos), function(i) {
  directorio_redes <- subtipos_directorios[[i]]
  archivos_redes <- list.files(path = directorio_redes, pattern = "TCGA-.*\\.tsv$", full.names = TRUE)
  
  resultados_subtipo <- bind_rows(lapply(archivos_redes, calc_cis_proportion, annot = annot))
  
  return(resultados_subtipo %>% mutate(Subtipo = subtipos[i]))
})

# Combine results into a single dataframe
resultados_proporciones_cis <- bind_rows(resultados_proporciones_cis)

# Function to split data into high and low CIS groups by median
split_cis_by_median <- function(subtipo_df) {
  mediana <- median(subtipo_df$Proporcion_CIS, na.rm = TRUE)
  
  alto_cis <- subtipo_df %>% filter(Proporcion_CIS >= mediana)
  bajo_cis <- subtipo_df %>% filter(Proporcion_CIS < mediana)
  
  list(alto_cis = alto_cis, bajo_cis = bajo_cis)
}

# Apply function to each subtype and save the results
lapply(subtipos, function(subtipo) {
  subtipo_df <- resultados_proporciones_cis %>% filter(Subtipo == subtipo)
  split_result <- split_cis_by_median(subtipo_df)
  
  # Save high and low CIS data to CSV
  write.csv(split_result$alto_cis, file = paste0(subtipo, "_alto_CIS.csv"), row.names = FALSE)
  write.csv(split_result$bajo_cis, file = paste0(subtipo, "_bajo_CIS.csv"), row.names = FALSE)
})

# Save overall CIS proportions for all subtypes
write.csv(resultados_proporciones_cis, "Proporciones_CIS_por_Subtipo.csv", row.names = FALSE)

cat("CIS proportion calculation completed and files have been saved.\n")

################################################################################
########################### SURVIVAL ANALYSIS ##################################
# This script performs survival analysis for high and low CIS proportion groups 
# using clinical data and generates Kaplan-Meier survival plots.
#
# Input: 
# - High and low CIS proportion files for each subtype.
# - Clinical data for TCGA BRCA patients.
#
# Output: 
# - Kaplan-Meier survival plots for each subtype.
# - Cox proportional hazards model for survival analysis.
#
################################################################################

library(dplyr)
library(survival)
library(survminer)

# Load clinical data
clinical_data <- read.csv("clinical_data_TCGA_BRCA_completo.csv", stringsAsFactors = FALSE)

# Prompt user to select CIS high and low files
cat("Seleccione el archivo de CIS alto para el subtipo de cáncer a analizar.\n")
alto_cis_file <- file.choose()

cat("Seleccione el archivo de CIS bajo para el mismo subtipo de cáncer.\n")
bajo_cis_file <- file.choose()

# Request subtype name from the user
subtipo <- readline(prompt = "Ingrese el nombre del subtipo que está analizando: ")

# Create a directory to save the survival plots
dir.create(subtipo, showWarnings = FALSE)

# Load CIS high and low data
alto_cis <- read.csv(alto_cis_file, stringsAsFactors = FALSE)
bajo_cis <- read.csv(bajo_cis_file, stringsAsFactors = FALSE)

# Merge clinical data with high and low CIS data
merged_alto <- merge(alto_cis, clinical_data, by.x = "TCGA_ID", by.y = "submitter_id")
merged_bajo <- merge(bajo_cis, clinical_data, by.x = "TCGA_ID", by.y = "submitter_id")

# Add a column to indicate CIS group
merged_alto <- merged_alto %>% mutate(CIS_Group = "Alto_CIS")
merged_bajo <- merged_bajo %>% mutate(CIS_Group = "Bajo_CIS")

# Combine high and low CIS datasets
combined_data <- bind_rows(merged_alto, merged_bajo)

# Ensure survival data is properly formatted
combined_data <- combined_data %>%
  mutate(
    overall_survival = as.numeric(overall_survival),  # Survival in days
    vital_status = ifelse(vital_status == "Alive", 0, 1)  # 0 = Alive, 1 = Deceased
  )

# Filter data for a maximum of 1825 days (5 years)
combined_data_filtered <- combined_data %>%
  filter(overall_survival <= 1825)

# Create a survival object
surv_object <- Surv(time = combined_data_filtered$overall_survival, event = combined_data_filtered$vital_status)

# Fit the Kaplan-Meier survival model
fit <- survfit(surv_object ~ CIS_Group, data = combined_data_filtered)

# Plot Kaplan-Meier survival curves and save them
plot <- ggsurvplot(
  fit, 
  data = combined_data_filtered, 
  pval = TRUE, 
  conf.int = TRUE,
  risk.table = TRUE, 
  risk.table.col = "strata",
  title = paste("", subtipo),
  xlab = "Days", 
  ylab = "Survival Probability",
  legend.labs = c("High CIS", "Low CIS"),
  palette = c("#E7B800", "#2E9FDF"),
  ggtheme = theme_minimal(),  
  font.title = c(14, "bold"),  
  font.x = c(12),  
  font.y = c(12),  
  font.tickslab = c(10),  
  font.legend = c(10)
)

# Display the plot
print(plot)

# Save the Kaplan-Meier plot as a JPEG file
ggsave(filename = file.path(subtipo, paste0(subtipo, "_5_years_Survival_Plot.jpeg")), plot = plot$plot, dpi = 300, device = "jpeg")

# Perform log-rank test
survdiff_test <- survdiff(surv_object ~ CIS_Group, data = combined_data_filtered)
print(survdiff_test)

# Fit the Cox proportional hazards model
cox_model <- coxph(surv_object ~ CIS_Group, data = combined_data_filtered)
print(cox_model)

# Check Cox model assumptions
cox_test <- cox.zph(cox_model)
print(cox_test)

################################################################################
