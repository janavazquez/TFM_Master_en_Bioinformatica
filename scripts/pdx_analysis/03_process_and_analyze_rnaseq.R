# ===================================================================
#
# TFM: Estudio de la Hipermetilación de BRCA1
#
# SCRIPT 03: PROCESADO Y ANÁLISIS DE DATOS DE RNA-seq DE LOS PDX
#
# AUTORA: Jana Vázquez Navarro
#
# DESCRIPCIÓN:
# Este script realiza un flujo de trabajo completo para los datos de RNA-seq
# de los modelos PDX. Incluye:
#   1. Carga de datos de conteos crudos y normalización a TPM (Transcritos Por Millón).
#   2. Carga de un segundo conjunto de datos de TPMs ya normalizados.
#   3. Integración de ambos conjuntos de datos.
#   4. Filtrado por un panel de genes de reparación del ADN.
#   5. Generación y guardado de un heatmap de expresión génica.
#
# ===================================================================


#### PASO 1: Cargar Librerías ####

# install.packages(c("readxl", "data.table", "dplyr", "tibble", "pheatmap"))
library(readxl)
library(data.table)
library(dplyr)
library(tibble)
library(pheatmap)

# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


#### PASO 2: Definir Rutas ####

input_dir <- "data/raw_data/pdx_rnaseq/"
output_dir_figs <- "results/figures/"
output_dir_data <- "data/processed_data/"


#### PASO 3: Procesar Datos de Conteos Crudos a TPM ####

# --- 3.1 Cargar y preparar base de datos de longitud de genes ---
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths_list <- sum(width(reduce(exons_by_gene)))
gene_lengths_df <- data.frame(
  gene_id = names(gene_lengths_list),
  length = as.vector(gene_lengths_list)
)
gene_lengths_df$ensembl_id <- mapIds(org.Hs.eg.db,
                                     keys = gene_lengths_df$gene_id,
                                     column = "ENSEMBL",
                                     keytype = "ENTREZID",
                                     multiVals = "first")
gene_lengths_final_df <- gene_lengths_df[!is.na(gene_lengths_df$ensembl_id), 
                                         c("ensembl_id", "length")]

# --- 3.2 Definir la función de normalización ---
counts_to_tpm <- function(counts, lengths) {
  lengths_kb <- lengths / 1000
  rpk <- sweep(counts, 1, lengths_kb, FUN = "/")
  scaling_factor <- colSums(rpk, na.rm = TRUE) / 1e6
  tpm <- sweep(rpk, 2, scaling_factor, FUN = "/")
  return(tpm)
}

# --- 3.3 Cargar y procesar matriz de conteos ---
file_counts <- file.path(input_dir, "COUNTS_SERRAVIO_03_04_09_10_11.xlsx")
counts_raw <- read_excel(file_counts)

counts <- counts_raw[-1, ]
colnames(counts) <- c("id_gene", "gene_name", "gene_type", colnames(counts_raw)[-c(1:3)])
ensembl_ids_clean <- gsub("\\..*$", "", counts$id_gene)
counts_matrix <- counts %>%
  select(starts_with("PDX")) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
rownames(counts_matrix) <- ensembl_ids_clean

# --- 3.4 Alinear conteos y longitudes ---
genes_in_counts <- rownames(counts_matrix)
lengths_for_matrix <- gene_lengths_final_df[gene_lengths_final_df$ensembl_id %in% genes_in_counts, ]
ordered_lengths_df <- lengths_for_matrix[match(genes_in_counts, lengths_for_matrix$ensembl_id), ]
genes_con_longitud <- !is.na(ordered_lengths_df$length)
ordered_lengths_vector <- ordered_lengths_df$length[genes_con_longitud]
counts_matrix_filtrada <- counts_matrix[genes_con_longitud, ]

# --- 3.5 Normalizar, transformar a log2 y mapear a símbolos de genes ---
tpm_matrix <- counts_to_tpm(counts_matrix_filtrada, ordered_lengths_vector)
log2_tpm_df <- log2(tpm_matrix + 1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensembl_id")

log2_tpm_df$gene_symbol <- mapIds(org.Hs.eg.db,
                                  keys = log2_tpm_df$ensembl_id,
                                  column = "SYMBOL",
                                  keytype = "ENSEMBL",
                                  multiVals = "first")

# Dataframe final del primer set de datos, con símbolos de gen como rownames
log2_tpm_from_counts <- log2_tpm_df %>%
  filter(!is.na(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  select(-ensembl_id) %>%
  column_to_rownames("gene_symbol")


#### PASO 4: Procesar Segundo Conjunto de Datos (TPMs pre-calculados) ####

file_tpms <- file.path(input_dir, "RNA_seq _021_unlogged_tpm.xlsx")
tpms_raw <- read_excel(file_tpms)

tpms_clean <- tpms_raw[-1, ]
colnames(tpms_clean) <- tpms_raw[1, ]

# Dataframe final del segundo set de datos, ya procesado y en log2
log2_tpm_from_excel <- tpms_clean %>%
  mutate(across(-hgnc_symbol, as.numeric)) %>%
  filter(!is.na(hgnc_symbol)) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  column_to_rownames(var = "hgnc_symbol") %>%
  mutate(across(everything(), ~ log2(.x + 1)))


#### PASO 5: Integrar Datos y Generar Heatmap ####

# --- 5.1 Definir genes de interés ---
repair_genes <- c("BRCA1", "BARD1", "PALB2", "BRCA2", "RAD51", "RAD51B", "RAD51C", 
                  "RAD51D", "RAD50", "MRE11A", "NBN", "ATM", "ATR", "CHEK1", "CHEK2", 
                  "MDC1", "FANCD2", "FANCI", "XRCC6", "XRCC5", "PRKDC", "LIG4", 
                  "XRCC4", "NHEJ1", "TP53BP1", "RAD52", "XRCC2", "XRCC3", "ERCC1", 
                  "ERCC4", "PNKP")

# --- 5.2 Alinear y combinar ambos dataframes ---
common_genes <- intersect(rownames(log2_tpm_from_counts), rownames(log2_tpm_from_excel))
common_genes <- intersect(common_genes, repair_genes) # Filtrar solo por los genes de reparación

# Alinear ambos dataframes por los genes comunes
df1_aligned <- log2_tpm_from_counts[common_genes, , drop = FALSE]
df2_aligned <- log2_tpm_from_excel[common_genes, , drop = FALSE]

# Combinar
combined_df <- cbind(
  df1_aligned[, c("PDX302_333-16_1R", "PDX270_225-15_2R")],
  df2_aligned[, c("PDX302", "PDX302OR2", "STG201")]
)

# Renombrar columnas para mayor claridad
colnames(combined_df) <- c("PDX302_count", "PDX270_count", "PDX302", "PDX302OR2", "PDX201")
# (Opcional) Si una muestra es redundante, se puede seleccionar una
combined_final <- combined_df %>% select(PDX302, PDX302OR2, PDX201, PDX270_count) %>%
  rename(PDX270 = PDX270_count)

# Guardar la matriz de datos final
fwrite(as.data.frame(combined_final) %>% rownames_to_column("gene"), 
       file.path(output_dir_data, "pdx_rnaseq_integrated_log2tpm.csv"))

# --- 5.3 Generar el Heatmap ---
heatmap_matrix <- as.matrix(combined_final)

# Crear anotación para las columnas
col_annotation <- data.frame(Estado = c("sensible", "resistente", "sensible", "resistente"))
rownames(col_annotation) <- colnames(heatmap_matrix)

# Guardar el heatmap en un archivo PNG
png(file.path(output_dir_figs, "Heatmap_PDX_RNAseq_RepairGenes.png"), width = 6, height = 8, units = "in", res = 300)

pheatmap(heatmap_matrix, 
         scale = "row", # Escalar por fila (gen) es más común para ver patrones de expresión
         annotation_col = col_annotation,
         main = "Expresión de Genes de Reparación en Modelos PDX",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 8,
         fontsize_col = 10)

dev.off()

message("Análisis de RNA-seq de PDX completado. El heatmap se ha guardado en 'results/figures/'.")