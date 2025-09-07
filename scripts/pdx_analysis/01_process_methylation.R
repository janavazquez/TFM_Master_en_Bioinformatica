# 01_process_methylation.R

# Objetivo: Leer todos los datos crudos de metilación de los PDX, 
# filtrarlos para la región de BRCA1 y generar un único archivo de datos limpios.

#### ---- Cargar librerias ---- ####

library(dplyr)
library(ggplot2)
library(ggh4x)
library(tidyverse)
library(stats)

#### ---- Lectura y filtrado de datos ---- ####

# from GC reports

leer_filtrar_anotar <- function(archivo, modelo) {
  read.table(file = archivo, sep = "\t", header = TRUE) %>%
    filter(rname == "chr17", pos >= 43124900, pos <= 43125600) %>%
    mutate(modelo = modelo)
}

archivos <- list(
  PDX201 = "PDX-VHIO-PDX201-TU01-DNA01.CG-report.tsv",
  PDX270 = "PDX-VHIO-PDX270-TU01-DNA01.CG-report.tsv",
  PDX302 = "PDX-VHIO-PDX302-TU01-DNA01.CG-report.tsv",
  PDX302OR2 = "PDX-VHIO-PDX302OR2-TU01-DNA01.CG-report.tsv",
  PDX621CNT = "PDX-VHIO-PDX621CNT-TU01-DNA01.CG-report.tsv",
  PDX621SR1 = "PDX-VHIO-PDX621SR1-TU01-DNA01.CG-report.tsv"
)

brca1_combinado <- bind_rows(
  lapply(names(archivos), function(nombre) {
    leer_filtrar_anotar(archivos[[nombre]], nombre)
  })
)

# from VEF data

data_pdx <- read.table("PDX-summary.tsv", sep = "\t", header = TRUE)
data_pdx[, `:=` (seqnames=factor(seqnames, levels=paste0("chr", c(1:22, "M", "X", "Y"))))]  # no funciona

pdx_vef <- data_pdx %>%
  dplyr::select(target, starts_with("VEF"))

pdx_vef <- pdx_vef[-1, ]
pdx_vef[, -1] <- lapply(pdx_vef[, -1], function(x) as.numeric(as.character(x)))
pdx_vef_clean <- pdx_vef %>%
  rowwise() %>%
  dplyr::filter(!all(is.na(c_across(-target))) && !all(c_across(-target) == 0, na.rm = TRUE)) %>%
  ungroup()

genes_hm <- pdx_vef_clean$target
models <- pdx_vef_clean %>%
  select(starts_with("VEF")) %>%
  colnames()

# Guardar los datos procesados para usarlos en los siguientes scripts
write.csv(brca1_interseccion, "data/processed_data/pdx_methylation_clean.csv")
