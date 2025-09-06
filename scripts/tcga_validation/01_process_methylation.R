# ===================================================================
#
# TFM: Estudio de la Hipermetilación de BRCA1
# 01_process_methylation.R
# AUTORA: Jana Vázquez Navarro
# DESCRIPCIÓN: Cargar y procesar los datos de metilación y los metadatos de TCGA.
# 
# ===================================================================

#### ---- cargar librerias ---- ####

library(readxl)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stats)

#### ---- read data and metadata ---- ####

archivo_metyl <- "HumanMethylation450" 
archivo_metadata <- "TCGA.BRCA.sampleMap-BRCA_clinicalMatrix"
archivo_coordenadas <- "HM450.hg38.manifest.tsv"
datos_methyl <- fread(archivo_metyl)
metadatos <- fread(archivo_metadata)
coordenadas <- fread(archivo_coordenadas)

## METADATA

#select interesting columns

colnames(metadatos)

metadata_clean <- metadatos %>%
  dplyr::select(1, 9, 14, 20, 21, 40, 56, 61, 92, 94, 105, 147, 161) %>%
  dplyr::rename(ajcc_stage = pathologic_stage) %>%
  dplyr::rename(er_status = breast_carcinoma_estrogen_receptor_status) %>%
  dplyr::rename(gender = Gender_nature2012) %>%
  dplyr::rename(her2_status = lab_proc_her2_neu_immunohistochemistry_receptor_status) %>%
  dplyr::rename(metastasis = Metastasis_Coded_nature2012) %>%
  dplyr::rename(pam50_RNAseq = PAM50Call_RNAseq) %>%
  dplyr::rename(pam50_mRNA = PAM50_mRNA_nature2012) %>%
  dplyr::rename(pr_status = breast_carcinoma_progesterone_receptor_status)

# crear la nueva columna 'ihc_classification' en metadata_clean
metadata_clean <- metadata_clean %>%
  mutate(
    classification_receptor = case_when(
      # Triple Negativo (TNBC): ER-, PR-, HER2-
      er_status == "Negative" & pr_status == "Negative" & her2_status == "Negative" ~ "Triple Negative",
      # HER2-enriched: ER-, PR-, HER2+
      er_status == "Negative" & pr_status == "Negative" & her2_status == "Positive" ~ "HER2-enriched",
      # Luminal B: ER+ y/o PR+, HER2+ (independientemente de Ki67)
      (er_status == "Positive" | pr_status == "Positive") & her2_status == "Positive" ~ "Luminal B (HER2+)",
      # Luminal A/B (HER2-): ER+ y/o PR+, HER2-
      # Aquí agrupo los Luminal A y B HER2-negativos. Para diferenciarlos,
      # necesitarías Ki67 y/o grado, que no tenemos en las columnas seleccionadas.
      (er_status == "Positive" | pr_status == "Positive") & her2_status == "Negative" ~ "Luminal A/B (HER2-)",
      # Otras combinaciones o si hay NAs en el estado de los receptores
      TRUE ~ "Indefinido/NA"
    )
  )

head(metadata_clean %>% dplyr::select(er_status, pr_status, her2_status, classification_receptor))
table(metadata_clean$classification_receptor)

# filtrar metadata solo con muestras de tumor primario
metadata_clean <- metadata_clean %>%
  filter(sample_type %in% c("Primary Tumor"))

## MANIFEST FILE

coord_clean <- coordenadas %>%
  dplyr::select(CpG_chrm, CpG_beg, CpG_end, Probe_ID, mapYD_A) %>%
  dplyr::filter(CpG_chrm == "chr17", CpG_beg > 43124900, CpG_beg < 43125600) %>%
  dplyr::rename(chr = CpG_chrm) %>%
  dplyr::rename(strand = mapYD_A)

ids_chr17 <- coord_clean$Probe_ID
coord_clean$position <- coord_clean$CpG_beg + 1

## METHYLATION DATA

# filter as much as possible
datos_methyl_clean <- datos_methyl %>%
  dplyr::select(sample, ends_with("01")) %>%
  dplyr::filter(sample %in% ids_chr17)

# transponer datos y editar matriz
met_t <- t(datos_methyl_clean)
met_t <- as.data.frame(met_t)

colnames(met_t) <- met_t[1, ]
met_t <- met_t[-1,]

# Filtrar metadata_clean para incluir solo las muestras comunes
metadata_filtered <- metadata_clean %>%
  filter(sampleID %in% commons)

# Filtrar met_t (datos de metilación transpuesta) para incluir solo las comunes
met_t_filtered <- met_t[rownames(met_t) %in% commons, ]

# asegurarse de que los valores de metilación sean numéricos
met_t_filtered <- met_t_filtered %>%
  mutate(across(everything(), as.numeric))

# alinear los metadatos filtrados al orden exacto de las muestras en met_t_filtered
metadata_aligned <- metadata_filtered[match(rownames(met_t_filtered), metadata_filtered$sampleID), ]

# Verificación final de que la alineación es perfecta
stopifnot(all(rownames(met_t_filtered) == metadata_aligned$sampleID))
message("¡Datos de metilación y metadatos alineados y filtrados correctamente!")

# Combinar los metadatos alineados con los datos de metilación transpuesta
methylation_data_merged <- cbind(metadata_aligned, met_t_filtered)

# comprobaciones
dim(methylation_data_merged)
head(methylation_data_merged[, 1:10])

write.csv(methylation_data_merged, "data/processed_data/tcga_methylation_clean.csv")