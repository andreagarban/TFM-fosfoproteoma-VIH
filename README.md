# TFM-fosfoproteoma-VIH

**Análisis proteómico y fosfoproteómico de PBMCs tratadas con inhibidores de tirosina quinasa como estrategia antiviral frente al VIH**

Trabajo de Fin de Máster en Bioinformática · Andrea Garban · Curso 2025–2026

## Descripción

Código en R y material suplementario del TFM. Se analiza el efecto de imatinib y dasatinib sobre el proteoma y fosfoproteoma de PBMCs de donantes sanos (TMT 11-plex), con foco en proteínas implicadas en el transporte intracelular y la entrada nuclear del VIH-1 (poro nuclear, complejo de preintegración, citoesqueleto, espliceosoma/procesamiento del ARN y señalización celular).

Se usa la hoja `Zpq` de `Results_Coiras_ProtAndPhosphoCorrect.xlsx` (fosfopéptidos corregidos por abundancia proteica).

| Condición | Muestras | Tratamiento | n |
|---|---|---|---|
| Control | 1–3 | PHA + IL-2 | 3 |
| Imatinib | 4–7 | Imatinib + PHA/IL-2 | 4 |
| Dasatinib | 8–11 | Dasatinib 75 nM + PHA/IL-2 | 4 |

Comparaciones: imatinib vs control y dasatinib vs control. Criterio de selección: p (t-test) < 0,05 y |AvA| > 1,3; el FDR (Benjamini-Hochberg) se calcula sobre todo el fosfoproteoma y se usa como capa de robustez interpretativa, no como filtro.

## Estructura

```
scripts/      Scripts R numerados en orden lógico de análisis
tablas/       Tablas de resultados (CSV/TXT) generadas por los scripts
figuras/      Figuras principales y heatmaps por categoría funcional
data_raw/     Bases de datos públicas NCBI (el Excel original no se incluye)
data_processed/  Objetos .rds intermedios (se regeneran al ejecutar los scripts)
```

## Pipeline

Ejecutar desde la raíz del proyecto (rutas relativas), en orden numérico. Cada script lee los `.rds` del anterior.

| Script | Contenido |
|---|---|
| `01_lectura_datos.R` | Lectura del Excel, renombrado de columnas, identificación de fosfopéptidos (`_Phospho`), eliminación de contaminantes |
| `02_PCA.R` | PCA de proteoma (Zq) y fosfoproteoma (Zpq) |
| `03_filtrado.R` | Filtrado p < 0,05 y \|AvA\| > 1,3; FDR (BH); dirección del cambio |
| `04_proteinas_VIH.R` | Selección por palabras clave y panel dirigido de proteínas VIH-relevantes |
| `05_cruce_bases_datos.R` | Cruce con la HIV-1 Human Interaction Database (niveles de evidencia 1–3) |
| `06_volcano.R` | Volcano plots |
| `07_tabla_uniprot.R` | Tabla integrada de resultados con niveles de evidencia |
| `08_tabla_fdr_dasatinib.R` | Fosfopéptidos de dasatinib que mantienen FDR < 0,10 |
| `09_heatmap_categorias_vih.R` | Heatmaps Zq y Zpq por categoría funcional |
| `10_go_kegg.R` | Enriquecimiento GO, KEGG y Reactome por fármaco y dirección |
| `11_terminos_interes_vih.R` | Términos enriquecidos relacionados con las categorías de interés |
| `12_cruce_go_ncbi.R` | Cruce de los genes de los términos con la evidencia NCBI |
| `13_comparacion_proteinas_virales.R` | Comparación de la distribución de proteínas virales |
| `14_string_red_global.R` | Redes STRING globales por fármaco |

`figuras/flujo_scripts.png` resume el flujo entre scripts.

## Requisitos

R 4.5.2. Paquetes principales: tidyverse 2.0.0, readxl, ggrepel, ComplexHeatmap, circlize, clusterProfiler, enrichplot, ReactomePA, GO.db, org.Hs.eg.db 3.21.0, AnnotationDbi 1.70.0, STRINGdb, igraph. Las versiones usadas en el enriquecimiento figuran en `tablas/10_versiones.txt`.

Nota: al cargar Bioconductor (`AnnotationDbi`, `S4Vectors`) se enmascaran verbos de dplyr; los scripts definen alias explícitos (`select <- dplyr::select`, etc.) al inicio.

## Datos

El Excel de partida pertenece al grupo que generó los datos y no se redistribuye aquí. Los objetos de `data_processed/` se generan al ejecutar los scripts.

## Licencia

Código bajo licencia MIT (ver `LICENSE`).

## Cita

Garban, A. (2026). *Análisis proteómico y fosfoproteómico de PBMCs tratadas con inhibidores de tirosina quinasa como estrategia antiviral frente al VIH* (Trabajo de Fin de Máster). Repositorio: TFM-fosfoproteoma-VIH.
