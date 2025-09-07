# ===================================================================
#
# TFM: Estudio de la Hipermetilación de BRCA1
# 02_analyze_and_plot_methylation.R
# AUTORA: Jana Vázquez Navarro
# DESCRIPCIÓN:
# Este script carga los datos de metilación procesados de la cohorte TCGA-BRCA
# y realiza los siguientes análisis y visualizaciones:
#   1. Análisis de Componentes Principales (PCA) para explorar la varianza.
#   2. Cálculo y visualización de la Entropía de Metilación para evaluar
#      la heterogeneidad epigenética.
#   3. Creación de gráfico para visualizar la distribución
#      de la metilación en cada CpG del promotor de BRCA1.
#
# ===================================================================

#### PASO 1: Cargar Librerías ####

library(dplyr)
library(ggplot2)
library(tibble)
library(data.table)
library(tidyr)

#### PASO 2: Definir Rutas y Cargar Datos Procesados ####

# --- Definir rutas de entrada y salida ---
input_dir <- "data/processed_data/"
output_dir_figs <- "results/figures/"

# --- Cargar los datos de metilación limpios y unidos con metadatos ---
# Este archivo fue generado por el script anterior: 01_process_methylation.R
file_clean_data <- file.path(input_dir, "tcga_methylation_clean.csv")
methylation_data_merged <- fread(file_clean_data)


#### PASO 3: Análisis de Componentes Principales (PCA) ####

message("Iniciando Análisis de Componentes Principales (PCA)...")

# --- 3.1 Preparar la matriz de datos para el PCA ---
# Seleccionamos solo las columnas que corresponden a los CpGs
cpg_columns <- colnames(methylation_data_merged)[grepl("^cg", colnames(methylation_data_merged))]
pca_input_matrix <- methylation_data_merged %>%
  select(all_of(cpg_columns)) %>%
  as.matrix()

# --- 3.2 Limpieza de la matriz: varianza e imputación ---
# Eliminar CpGs con varianza cero para evitar errores en prcomp
variances <- apply(pca_input_matrix, 2, var, na.rm = TRUE)
zero_var_cols <- which(variances < 1e-10)
if (length(zero_var_cols) > 0) {
  pca_input_matrix <- pca_input_matrix[, -zero_var_cols]
}

# Imputar valores NA con la media de la columna (del CpG)
for(i in 1:ncol(pca_input_matrix)){
  pca_input_matrix[is.na(pca_input_matrix[,i]), i] <- mean(pca_input_matrix[,i], na.rm = TRUE)
}

# --- 3.3 Ejecutar PCA ---
pca_results <- prcomp(pca_input_matrix, center = TRUE, scale. = TRUE)

# --- 3.4 Preparar datos para la visualización ---
pca_scores <- as.data.frame(pca_results$x)
pca_scores$classification_receptor <- methylation_data_merged$classification_receptor

variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2)
pc1_var <- round(variance_explained[1] * 100, 1)
pc2_var <- round(variance_explained[2] * 100, 1)

# --- 3.5 Generar y guardar el gráfico PCA ---
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = classification_receptor)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA de Metilación en el Promotor de BRCA1 (TCGA-BRCA)",
    x = paste0("PC1 (", pc1_var, "% de varianza)"),
    y = paste0("PC2 (", pc2_var, "% de varianza)"),
    color = "Subtipo Molecular"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set1")

ggsave(
  filename = file.path(output_dir_figs, "PCA_TCGA_Methylation.png"),
  plot = pca_plot,
  width = 10, height = 7
)

message("PCA completado y guardado.")


#### PASO 4: Análisis de Entropía de Metilación ####

message("Iniciando análisis de Entropía de Metilación...")

# --- 4.1 Definir la función de entropía ---
calculate_entropy <- function(p) {
  p <- pmax(pmin(p, 1 - 1e-9), 1e-9)
  - (p * log2(p) + (1 - p) * log2(1 - p))
}

# --- 4.2 Calcular la entropía para cada CpG y el promedio por muestra ---
entropy_matrix <- apply(pca_input_matrix, 2, calculate_entropy)
sample_avg_entropy <- data.frame(
  sampleID = methylation_data_merged$sampleID,
  average_entropy = rowMeans(entropy_matrix, na.rm = TRUE),
  classification_receptor = methylation_data_merged$classification_receptor
)

# --- 4.3 Filtrar para el subtipo Triple Negativo ---
tn_entropy_data <- sample_avg_entropy %>%
  filter(classification_receptor == "Triple Negative")

# --- 4.4 Generar y guardar el histograma de entropía ---
entropy_histogram <- ggplot(tn_entropy_data, aes(x = average_entropy)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 20,
    fill = "cyan4",
    alpha = 0.8
  ) +
  geom_density(color = "black", linewidth = 1) +
  labs(
    title = "Distribución de la Entropía de Metilación en Tumores Triple Negativo",
    x = "Entropía de Metilación Promedio",
    y = "Densidad"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path(output_dir_figs, "Histogram_TCGA_Entropy_TNBC.png"),
  plot = entropy_histogram,
  width = 8, height = 6
)

message("Análisis de entropía completado y guardado.")


#### PASO 5: Gráficos de Violín y Jitter por CpG (Subtipo TNBC) ####

message("Generando gráficos de distribución por CpG...")

# --- 5.1 Preparar los datos en formato largo ---
tn_methylation_data <- methylation_data_merged