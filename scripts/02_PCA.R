# ==============================================================================
# 02_PCA.R
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con ITKs
#      como estrategia antiviral frente al VIH
# ==============================================================================

library(tidyverse)
library(ggrepel)

# Alias explicitos para evitar conflictos con otros paquetes
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
slice  <- dplyr::slice

zpq_clean <- readRDS("data_processed/zpq_clean.rds")


# 1. Correspondencia canal TMT -> muestra / condición

mapping <- tibble(
  canal     = c("126t", "127N", "127C", "128N", "128C",
                "129N", "129C", "130N", "130C", "131N", "131C"),
  muestra   = paste0("Muestra_", 1:11),
  condicion = c(rep("Control", 3), rep("Imatinib", 4), rep("Dasatinib", 4))
)


# 2. Función: preparar matriz muestras x péptidos y correr prcomp()
#    - Identificador de fila: RowProtein no es único (varios péptidos por
#      proteína), se combina con el número de fila.
#    - NAs: se filtran péptidos con algún NA en las 11 muestras
#      (complete.cases). Se informa del % de péptidos retenidos.
#    - center = TRUE, scale. = FALSE: los valores Zpq ya están en escala
#      comparable entre péptidos; escalar daría el mismo peso a péptidos
#      con señal plana (ruido) que a los que cambian biológicamente.

run_pca <- function(df, canal_prefix = "Zpq_") {

  cols_muestras <- paste0(canal_prefix, mapping$canal)

  df <- df %>% mutate(row_id = paste0(RowProtein, "_", row_number()))

  mat <- df %>% select(all_of(cols_muestras)) %>% as.matrix()
  rownames(mat) <- df$row_id

  n_total <- nrow(mat)
  mat_completa <- mat[complete.cases(mat), ]
  n_completa <- nrow(mat_completa)

  cat("Péptidos totales:", n_total, "\n")
  cat("Péptidos completos (sin NA), usados en el PCA:", n_completa,
      sprintf(" (%.1f%%)\n", 100 * n_completa / n_total))

  # Transponer: muestras en filas, péptidos en columnas
  mat_t <- t(mat_completa)

  pca <- prcomp(mat_t, center = TRUE, scale. = FALSE)

  list(pca = pca, n_peptidos = n_completa, canal_prefix = canal_prefix)
}


# 3. Ejecutar PCA: (a) proteoma (Zq, un valor por proteína),
#                  (b) fosfoproteoma (Zpq, solo fosfopéptidos)
# Zq es cuantificación a nivel de proteína: en zpq_clean aparece repetido en
# todas las filas de péptido de una misma proteína. Antes del
# PCA hay que quedarnos con una fila por proteína.

cols_zq <- paste0("Zq_", mapping$canal)
proteoma_dedup <- zpq_clean %>%
  mutate(.k  = do.call(paste, c(across(all_of(cols_zq)), sep = "|")),
         .na = if_any(all_of(cols_zq), is.na)) %>%
  group_by(protein, .k) %>% mutate(.n = n()) %>% ungroup() %>%
  arrange(protein, .na, desc(.n)) %>%     # sin NA primero, luego el valor mas frecuente
  distinct(protein, .keep_all = TRUE) %>%
  select(-.k, -.n, -.na)

pca_proteoma <- run_pca(proteoma_dedup, canal_prefix = "Zq_")
pca_phospho  <- run_pca(zpq_clean %>% filter(is_phospho == TRUE), canal_prefix = "Zpq_")


# 4. Función de plotting

plot_pca <- function(resultado, titulo) {

  scores <- as.data.frame(resultado$pca$x) %>%
    rownames_to_column("canal_col") %>%
    mutate(canal = str_remove(canal_col, paste0("^", resultado$canal_prefix))) %>%
    left_join(mapping, by = "canal")

  var_exp <- round(100 * summary(resultado$pca)$importance[2, 1:2], 1)

  ggplot(scores, aes(x = PC1, y = PC2, color = condicion, label = muestra)) +
    geom_point(size = 4) +
    geom_text_repel(show.legend = FALSE, size = 3.5) +
    scale_color_manual(values = c("Control"   = "grey40",
                                   "Imatinib"  = "#E69F00",
                                   "Dasatinib" = "#56B4E9")) +
    labs(
      title    = titulo,
      subtitle = paste0("n = ", resultado$n_peptidos, " péptidos"),
      x        = paste0("PC1 (", var_exp[1], "%)"),
      y        = paste0("PC2 (", var_exp[2], "%)"),
      color    = "Condición"
    ) +
    theme_minimal(base_size = 13)
}

plot_pca(pca_proteoma, "PCA - Proteoma (Zq)")
plot_pca(pca_phospho,  "PCA - Fosfoproteoma (Zpq)")


#Valores

#Varianza explicada (% PC1 y PC2)
round(100 * summary(pca_proteoma$pca)$importance[2, 1:2], 1)
round(100 * summary(pca_phospho$pca)$importance[2, 1:2], 1)

#Coordenadas por muestra
tabla_scores <- function(res) {
  as.data.frame(res$pca$x[, 1:2]) %>%
    rownames_to_column("canal_col") %>%
    mutate(canal = str_remove(canal_col, paste0("^", res$canal_prefix))) %>%
    left_join(mapping, by = "canal") %>%
    transmute(muestra, condicion, PC1 = round(PC1, 1), PC2 = round(PC2, 1))
}
tabla_scores(pca_proteoma)
tabla_scores(pca_phospho)

#Zq debe ser identico en todas las filas de una misma proteina
cols_zq <- paste0("Zq_", mapping$canal)
zpq_clean %>%
  mutate(k = do.call(paste, c(across(all_of(cols_zq)), sep = "|"))) %>%
  group_by(protein) %>%
  summarise(n_valores = n_distinct(k)) %>%
  count(n_valores)


# 5. Guardar 

saveRDS(pca_proteoma, "data_processed/pca_proteoma.rds")
saveRDS(pca_phospho,  "data_processed/pca_phospho.rds")
