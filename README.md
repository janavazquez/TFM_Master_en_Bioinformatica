# Estudio de la Hipermetilación de BRCA1 en Modelos de Cáncer de Mama Resistentes a la Inhibición de PARP

**Autora:** Jana Vázquez Navarro
**Trabajo de Fin de Máster | Máster Universitario en Bioinformática | UNIR**

---

### Descripción del Proyecto

Este repositorio contiene el código y los análisis realizados para el Trabajo de Fin de Máster (TFM) centrado en investigar los mecanismos de resistencia adquirida a los inhibidores de PARP (PARPi) en modelos de cáncer de mama.

El silenciamiento epigenético del gen supresor de tumores *BRCA1* mediante hipermetilación de su promotor confiere sensibilidad a los PARPi. Sin embargo, la adquisición de resistencia, a menudo asociada a la reversión de estos cambios epigenéticos, es un obstáculo clínico fundamental.

Este proyecto realiza un **análisis bioinformático multiómico** para:
1.  Evaluar los cambios en la metilación del promotor de *BRCA1* en modelos de xenoinjertos derivados de pacientes (PDX) sensibles y con resistencia adquirida a PARPi.
2.  Correlacionar el estado de metilación con los perfiles de expresión génica (RNA-seq).
3.  Utilizar métricas avanzadas como la **entropía de metilación** para caracterizar la heterogeneidad epigenética.
4.  Validar los hallazgos preclínicos utilizando la cohorte de pacientes de The Cancer Genome Atlas (TCGA-BRCA).

El principal hallazgo del estudio es que la **desmetilación del promotor de *BRCA1*** es un mecanismo molecular clave que impulsa la resistencia adquirida a los PARPi, lo que se manifiesta no solo como una reactivación transcripcional, sino también como un cambio drástico de un paisaje epigenético ordenado (baja entropía) a uno caótico y heterogéneo (alta entropía).

### Estructura del Repositorio

El proyecto está organizado en las siguientes carpetas:

-   `/data`: Contiene los datos de entrada.
    -   `/raw_data`: Archivos originales, separados por tipo y origen (PDX/TCGA, Metilación/RNA-seq).
    -   `/processed_data`: Archivos intermedios generados por los scripts de procesamiento.
-   `/scripts`: Contiene los scripts de R para realizar todo el análisis, organizados por cohorte.
    -   `/pdx_analysis`: Scripts para el análisis de los modelos PDX.
    -   `/tcga_validation`: Scripts para el análisis de la cohorte de validación de TCGA.
-   `/results`: Contiene los resultados finales generados por los scripts.
    -   `/figures`: Gráficos, heatmaps e histogramas.
    -   `/tables`: Tablas con resultados numéricos.
-   `README.md`: Este archivo.

### Requisitos e Instalación

Para ejecutar los análisis, se necesita un entorno de R (versión >= 4.0) con las siguientes librerías instaladas.

**Paquetes de CRAN:**
```R
install.packages(c("ggplot2", "data.table", "dplyr", "pheatmap", "ggh4x", "tidyverse", "readxl", "corrplot", "gridExtra", "tibble"))
```

**Paquetes de Bioconductor:**
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "epialleleR"))
```

### Flujo de Trabajo y uso

> **⚠️ NOTA IMPORTANTE SOBRE LOS DATOS**
>
> Debido a la naturaleza sensible de los datos derivados de pacientes, los archivos de datos brutos **no se incluyen en este repositorio**. Para poder ejecutar los scripts de análisis, primero deberás obtener los datos de sus fuentes originales y colocarlos en la carpeta `/data/raw_data/` siguiendo la estructura descrita. Consulta el archivo `DESCRIPCION_DE_DATOS.md` para obtener más detalles sobre el origen y la estructura de los datos necesarios.

Los análisis se pueden reproducir ejecutando los scripts de la carpeta `/scripts` en el siguiente orden:

#### 1. Análisis de los Modelos PDX (`/scripts/pdx_analysis/`)
1.  **`01_process_methylation.R`**: Procesa los archivos `.CG-report.tsv` para filtrar la región del promotor de *BRCA1* y calcular los valores beta de metilación.
2.  **`02_exploratory_analysis.R`**: Realiza análisis exploratorios como PCA y heatmaps sobre los datos de metilación (VEF).
3.  **`03_calculate_entropy.R`**: Utiliza los datos de metilación procesados para calcular la entropía de metilación por CpG y generar los gráficos de distribución.

#### 2. Análisis de Validación en TCGA (`/scripts/tcga_validation/`)
1.  **`01_process_methylation.R`**: Carga, limpia y une los datos de metilación y los metadatos clínicos de la cohorte TCGA-BRCA.
2.  **`02_analyze_methylation_PCA.R`**: Realiza el PCA sobre los datos de metilación de TCGA.
3.  **`03_calculate_entropy.R`**: Calcula la entropía de metilación para las muestras de TCGA y genera los gráficos de distribución.

Asegúrate de que los datos brutos se encuentran en las subcarpetas correspondientes de `/data/raw_data/` antes de ejecutar los scripts. Los resultados se guardarán automáticamente en la carpeta `/results/`.

### Cita

Si utilizas este código o los hallazgos de este trabajo, por favor, cita el siguiente documento:

> Vázquez Navarro, J. (2025). *Estudio de hipermetilación del gen BRCA1 en modelos de cáncer de mama resistentes a la inhibición de PARP* (Trabajo de Fin de Máster). Universidad Internacional de La Rioja, La Rioja, España.

---
