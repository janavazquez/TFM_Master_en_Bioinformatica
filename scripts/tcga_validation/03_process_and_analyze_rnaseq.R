# ===================================================================
#
# TFM: Estudio de la Hipermetilación de BRCA1
#
# 03: PROCESADO Y ANÁLISIS DE DATOS DE RNA-seq DE TCGA
#
# AUTORA: Jana Vázquez Navarro
#
# DESCRIPCIÓN:
# Este script carga los datos de expresión (RNA-seq) y los metadatos clínicos
# de la cohorte TCGA-BRCA. Realiza los siguientes pasos:
#   1. Limpia y formatea los metadatos clínicos.
#   2. Crea una clasificación de subtipos moleculares (p. ej., Triple Negativo).
#   3. Une los datos de expresión con los metadatos.
#   4. Realiza un Análisis de Componentes Principales (PCA) sobre un panel
#      de genes de reparación del ADN.
#   5. Genera y guarda las figuras del PCA.
#
# ===================================================================


#### PASO 1: Cargar Librerías ####

# install.packages(c("data.table", "dplyr", "ggplot2", "stats"))
library(data.table)
library(dplyr)
library(ggplot2)
library(stats)

#### PASO 2: Definir Rutas y Cargar Datos ####

# --- Definir rutas de entrada y salida ---
# Esto hace que el script sea más fácil de adaptar en el futuro.
input_dir <- "data/raw_data/tcga_rnaseq/"
output_dir_figs <- "results/figures/"
output_dir_data <- "data/processed_data/"

# --- Cargar datos de expresión y metadatos ---
# Usamos fread por su eficiencia con archivos grandes.
file_expr <- file.path(input_dir, "HiSeqV2")
file_meta <- file.path(input_dir, "TCGA.BRCA.sampleMap-BRCA_clinicalMatrix")

datos_expr <- fread(file_expr)
metadatos <- fread(file_meta)


#### PASO 3: Limpieza y Procesamiento de Metadatos ####

# --- Seleccionar y renombrar columnas de interés ---
metadata_clean <- metadatos %>%
  select(sampleID, sample_type, pathologic_stage, 
         breast_carcinoma_estrogen_receptor_status,
         breast_carcinoma_progesterone_receptor_status,
         lab_proc_her2_neu_immunohistochemistry_receptor_status) %>%
  rename(
    ajcc_stage = pathologic_stage,
    er_status = breast_carcinoma_estrogen_receptor_status,
    pr_status = breast_carcinoma_progesterone_receptor_status,
    her2_status = lab_proc_her2_neu_immunohistochemistry_receptor_status
  )

# --- Crear la columna de clasificación molecular ---
metadata_clean <- metadata_clean %>%
  mutate(
    classification_receptor = case_when(
      er_status == "Negative" & pr_status == "Negative" & her2_status == "Negative" ~ "Triple Negative",
      er_status == "Negative" & pr_status == "Negative" & her2_status == "Positive" ~ "HER2-enriched",
      (er_status == "Positive" | pr_status == "Positive") & her2_status == "Positive" ~ "Luminal B (HER2+)",
      (er_status == "Positive" | pr_status == "Positive") & her2_status == "Negative" ~ "Luminal A/B (HER2-)",
      TRUE ~ "Indefinido/NA"
    )
  )

# --- Filtrar para quedarse solo con tumores primarios ---
metadata_clean <- metadata_clean %>%
  filter(sample_type == "Primary Tumor")


#### PASO 4: Unir Datos de Expresión y Metadatos ####

# --- Transponer la matriz de expresión ---
datos_expr_t <- as.data.frame(t(datos_expr))
colnames(datos_expr_t) <- datos_expr_t[1, ]
datos_expr_t <- datos_expr_t[-1, ]
# Convertir rownames a una columna para poder unir los datos
datos_expr_t <- tibble::rownames_to_column(datos_expr_t, var = "sampleID")


# --- Encontrar y filtrar por las muestras comunes ---
common_samples <- intersect(datos_expr_t$sampleID, metadata_clean$sampleID)

metadata_final <- metadata_clean %>% filter(sampleID %in% common_samples)
datos_expr_final <- datos_expr_t %>% filter(sampleID %in% common_samples)

# --- Alinear y combinar ---
# Se ordena la metadata para que coincida con el orden de la matriz de expresión
metadata_final <- metadata_final[match(datos_expr_final$sampleID, metadata_final$sampleID), ]

# Verificación
stopifnot(all(datos_expr_final$sampleID == metadata_final$sampleID))

# Combinación final
datos_merged <- cbind(metadata_final, datos_expr_final[, -1]) # Se excluye la columna sampleID duplicada

# Guardar el archivo procesado para uso futuro (opcional pero recomendado)
fwrite(datos_merged, file.path(output_dir_data, "tcga_rnaseq_clinical_merged.csv"))


#### PASO 5: Análisis de Componentes Principales (PCA) ####

# --- Definir el panel de genes de interés ---
repair_genes <- c("BRCA1", "BARD1", "PALB2", "BRCA2", "RAD51", "RAD51B", "RAD51C", 
                  "RAD51D", "RAD50", "MRE11A", "NBN", "ATM", "ATR", "CHEK1", "CHEK2", 
                  "MDC1", "FANCD2", "FANCI", "XRCC6", "XRCC5", "PRKDC", "LIG4", 
                  "XRCC4", "NHEJ1", "TP53BP1", "RAD52", "XRCC2", "XRCC3", "ERCC1", 
                  "ERCC4", "PNKP")

# --- Preparar la matriz de datos para el PCA ---
genes_in_data <- intersect(repair_genes, colnames(datos_merged))
expr_matrix_for_pca <- datos_merged %>%
  select(all_of(genes_in_data)) %>%
  mutate(across(everything(), as.numeric))

# --- Ejecutar PCA ---
pca_results <- prcomp(expr_matrix_for_pca, center = TRUE, scale. = TRUE)

# --- Preparar datos para la visualización ---
pca_df <- as.data.frame(pca_results$x)
pca_df$classification_receptor <- datos_merged$classification_receptor

# Calcular la varianza explicada para las etiquetas de los ejes
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2)
pc1_label <- paste0("PC1 (", round(variance_explained[1] * 100, 1), "%)")
pc2_label <- paste0("PC2 (", round(variance_explained[2] * 100, 1), "%)")


#### PASO 6: Generar y Guardar Gráficos PCA ####

# --- Gráfico 1: PCA con todos los subtipos ---
pca_plot_all <- ggplot(pca_df, aes(x = PC1, y = PC2, color = classification_receptor)) + 
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA de Genes de Reparación del ADN en TCGA-BRCA",
    subtitle = "Todos los subtipos moleculares",
    x = pc1_label, 
    y = pc2_label, 
    color = 'Clasificación Molecular'
  ) +
  theme_light(base_size = 14) +
  scale_color_brewer(palette = "Set1")

# --- Gráfico 2: PCA solo de Tumores Triple Negativo ---
pca_df_tnbc <- pca_df %>% filter(classification_receptor == "Triple Negative")

pca_plot_tnbc <- ggplot(pca_df_tnbc, aes(x = PC1, y = PC2)) + 
  geom_point(size = 3, color = "black") +
  labs(
    title = "PCA de Genes de Reparación del ADN en TCGA-BRCA",
    subtitle = "Subtipo Triple Negativo",
    x = pc1_label, 
    y = pc2_label
  ) +
  theme_light(base_size = 14)

# --- Guardar las figuras ---
ggsave(
  filename = file.path(output_dir_figs, "PCA_TCGA_RNAseq_All_Subtypes.png"),
  plot = pca_plot_all,
  width = 10, height = 7
)

ggsave(
  filename = file.path(output_dir_figs, "PCA_TCGA_RNAseq_TNBC_Only.png"),
  plot = pca_plot_tnbc,
  width = 8, height = 6
)

message("Análisis de RNA-seq de TCGA completado. Las figuras se han guardado en 'results/figures/'.")