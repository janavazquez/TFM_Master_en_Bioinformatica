# 02_analyze_and_plot_methylation.R
# Objetivo: Cargar los datos de metilación limpios y realizar todos los análisis 
# y visualizaciones (PCA, entropía, heatmap).

#### ---- Cargar librerias ---- ####

library(stats) 
library(dplyr)
library(data.table)
library(tibble)
library(ggplot2)

#### ---- Cargar datos ---- #####

pdx_methylation_clean <- read.csv("data/processed_data/pdx_methylation_clean.csv")

#### ---- preparar datos para PCA ---- ####

vef_matrix <- pdx_methylation_clean[, -1]
pdx_t <- t(vef_matrix)

pdx_t_clean <- pdx_t[, colSums(is.na(pdx_t) | is.infinite(as.matrix(pdx_t))) == 0]

#### ---- PCA ---- ####

pca_result <- prcomp(pdx_t_clean, scale. = TRUE)

pca.df <- data.frame(pca_result$x)
pca.df$models <- models

pca.df$models <- factor(pca.df$models,
                        levels = c("VEF_PDX201","VEF_PDX270","VEF_PDX302","VEF_PDX302OR2","VEF_PDX621CNT","VEF_PDX621SR1"),  # los valores originales
                        labels = c("PDX201","PDX270","PDX302","PDX302OR2","PDX621CNT","PDX621SR1"))  # los nombres que quieres que aparezcan

varianzas <- pca_result$sdev^2 
total.varianza <- sum(varianzas)
varianza.explicada <- varianzas/total.varianza 
varianza.acumulada <- cumsum(varianza.explicada) 
n.pc <- min(which(varianza.acumulada > 0.9)) 


x_label <- paste0(paste('PC1', round(varianza.explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(varianza.explicada[2] * 100, 2)), '%')


p <- ggplot(pca.df, aes(x=PC1, y=PC2, color=models)) +
  geom_point(size=3) +
  labs(title='PCA', x=x_label, y=y_label, color = models) +
  guides(color = guide_legend(title = "Model")) +
  theme_classic() +
  theme(panel.grid.major = element_line(color="gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title = element_text(hjust = 0.5))

p

#### ---- calculo de la entropia ---- ####

calculate_entropy <- function(p) {
  p <- pmax(pmin(p, 1 - 1e-9), 1e-9)
  - (p * log2(p) + (1 - p) * log2(1 - p))
}

data <- brca1_interseccion %>%
  select(pos, modelo, beta_value)

# Transformar los datos de formato largo a ancho

pdx_input_wide <- data %>%
  pivot_wider(
    names_from = modelo,
    values_from = beta_value
  )

# Preparar la Matriz de Datos 

pdx_input_final <- pdx_input_wide %>%
  column_to_rownames(var = "pos")

# Aplicar la función de entropía a la matriz limpia de PDX

pdx_entropy_matrix <- apply(pdx_input_final, 2, calculate_entropy)
rownames(pdx_entropy_matrix) <- rownames(pdx_input_clean)

# Calcular la entropía promedio para cada modelo PDX
pdx_avg_entropy <- data.frame(
  sampleID = rownames(pdx_entropy_matrix),
  average_entropy = rowMeans(pdx_entropy_matrix, na.rm = TRUE)
)

print(head(pdx_avg_entropy))

pdx_entropy_df <- as.data.frame(pdx_entropy_matrix)

# Pivotar de formato ancho a largo con `pivot_longer`

pdx_entropy_long <- pdx_entropy_df %>%
  pivot_longer(
    cols = everything(),
    names_to = "modelo",
    values_to = "entropy"
  )

# Calcular la mediana para cada modelo 

medians_per_model <- pdx_entropy_long %>%
  group_by(modelo) %>%
  summarise(median_entropy = median(entropy, na.rm = TRUE))

# determinar los colores
colores_manuales <- c(
  "PDX201" = "#D55E00",
  "PDX270" = "orchid",
  "PDX302" = "#009E73",
  "PDX302OR2" = "green3",
  "PDX621CNT" = "skyblue",
  "PDX621SR1" = "#56B4E9"
)

# Crear el gráfico con facetas 

p2 <- ggplot(pdx_entropy_long, aes(x = entropy)) +
  
  # 1. Histograma. Coloreamos cada uno de forma distinta.
  geom_histogram(
    aes(y = after_stat(density), fill = modelo),
    bins = 30,
    alpha = 0.8
  ) +
  
  # curva de densidad para ver la forma de la distribución.
  geom_density(color = "black") +
  
  # `facet_wrap` crea un panel para cada modelo.
  #    - `scales = "free_y"` permite que el eje Y de cada gráfico se ajuste a sus propios datos.
  facet_wrap(~ modelo, 
             scales = "free_y",
             ncol = 2) +
  scale_fill_manual(values = colores_manuales) +
  
  # 3. Etiquetas y tema.
  labs(
    x = "Entropía de Metilación por CpG",
    y = "Densidad"
  ) +
  theme_light(base_size = 14) +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "none") 

p2

#### ---- Heatmap ---- ####

# variables

pdxs <- c("VEF_PDX201", "VEF_PDX270","VEF_PDX302", "VEF_PDX302OR2", "VEF_PDX621CNT","VEF_PDX621SR1")

# Methylated genes for each model

for(model in pdxs) {
  positive_genes <- data_pdx$target[data_pdx[[model]] > 0]
  cat("\nModel:", model, "- Genes with >0 values:\n")
  #print(unique(positive_genes))
  print(length(positive_genes))
}

# Descriptive statistics

summary(pdx_vef_clean)

# 1. Lista de genes de reparación (sin duplicados)
repair_genes <- c("BRCA1", "BARD1", "PALB2", "BRCA2", "RAD51", "RAD51B", 
                  "RAD51C", "RAD51D", "RAD50", "MRE11A", "NBN", "ATM", "ATR", 
                  "CHEK1", "CHEK2", "MDC1", "FANCD2", "FANCI", "XRCC6", "XRCC5", 
                  "PRKDC", "LIG4", "XRCC4", "NHEJ1", "TP53BP1", "RAD52", "XRCC2", 
                  "XRCC3", "ERCC1", "ERCC4", "PNKP")

repair_genes <- unique(repair_genes)

# 2. Filtrar los genes que están presentes en pdx_vef$target
repair_present <- intersect(repair_genes, pdx_vef$target)

# 3. Filtrar la matriz de expresión para esos genes
# Suponiendo que los nombres de fila en gene_matrix_t son los nombres de los genes
filtered_repair_df <- pdx_vef[pdx_vef$target %in% repair_present, ]
rownames(filtered_repair_df) <- make.unique(as.character(filtered_repair_df$target))
filtered_repair_df$target <- NULL


# 4. Preparar matriz numérica para pheatmap
# Asegurarte de que todo es numérico y que no hay NAs
filtered_repair_df <- as.data.frame(filtered_repair_df)
filtered_repair_df[] <- lapply(filtered_repair_df, function(x) as.numeric(as.character(x)))
filtered_repair_df[is.na(filtered_repair_df)] <- 0
hm_matrix <- as.matrix(filtered_repair_df)

# 5. Etiquetas
genes_hm <- rownames(hm_matrix)
ordered_samples <- c("VEF_PDX621CNT", "VEF_PDX621SR1", "VEF_PDX302", 
                     "VEF_PDX302OR2", "VEF_PDX201", "VEF_PDX270")
# Reordenar columnas
hm_matrix <- hm_matrix[, ordered_samples]

# 6. Dibujar el heatmap
pheatmap(hm_matrix, 
         scale = "column",
         main = "Heatmap de genes de reparación del ADN",
         labels_row = genes_hm,
         cluster_rows = FALSE,
         color = colorRampPalette(c("white", "cyan4"))(100),
         fontsize_row = 6,
         fontsize_col = 8,
         legend_breaks = 0:3,
         legend_labels = c("0", "1e-4", "1e-3", "1e-2"),
         labels_col = c("PDX621CNT", "PDX621SR1", "PDX302", "PDX302OR2", "PDX201", "PDX270"))