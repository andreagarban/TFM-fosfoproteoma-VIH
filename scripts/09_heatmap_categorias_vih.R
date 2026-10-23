# ==============================================================================
# 09_heatmap_categorias_vih.R
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con ITKs
#      como estrategia antiviral frente al VIH
#
# Selección de proteínas (las tres condiciones a la vez):
#   (1) >= 1 fosfopéptido con p < 0.05 y |AvA| > 1.3 en imatinib o dasatinib
#   (2) relación documentada con VIH-1 en el cruce NCBI (Nivel 1 y 2)
#   (3) pertenencia a >= 1 de 5 categorías funcionales (GO principal, UniProt/
#       panel dirigido como complemento)
# ==============================================================================

#0. Paquetes y alias
suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(ComplexHeatmap)
  library(circlize)
})

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
slice  <- dplyr::slice


#1. Configuración 
umbral_p   <- 0.05
umbral_ava <- 1.3

RDS_FOSFO <- "data_processed/fosfoproteoma.rds"
RDS_TABLA <- "data_processed/tabla_uniprot_niveles.rds"

col_peptido <- "sequence"  
niveles_vih <- c("1. Relación directa con VIH-1",
                 "2. Relación indirecta o funcional")

canales   <- c("126t","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C")
muestras  <- paste0("Muestra_", 1:11)
condicion <- factor(c(rep("Control", 3), rep("Imatinib", 4), rep("Dasatinib", 4)),
                    levels = c("Control", "Imatinib", "Dasatinib"))

solo_sig <- TRUE

# Orden de prioridad para la categoría principal (de más específica a más general)
ORDEN_CAT <- c("Poro/transporte nuclear",
               "Complejo de preintegración",
               "Spliceosoma/procesamiento ARN",
               "Citoesqueleto/tráfico intracelular",
               "Señalización celular")


# 2. Carga y comprobaciones 
fosfo_all <- readRDS(RDS_FOSFO)
tabla     <- readRDS(RDS_TABLA)
prot      <- fosfo_all
fosfo     <- fosfo_all

cols_zpq <- paste0("Zpq_", canales)
cols_zq  <- paste0("Zq_", canales)
cols_stat_fosfo <- c("Zpq_AvA_imatinib", "Zpq_Ttest_imatinib",
                     "Zpq_AvA_dasatinib", "Zpq_Ttest_dasatinib")

comprobar <- function(df, cols, nombre) {
  falta <- setdiff(cols, names(df))
  if (length(falta) > 0) {
    cat("\nColumnas disponibles en", nombre, ":\n")
    print(names(df))
    stop("Faltan columnas en ", nombre, ": ", paste(falta, collapse = ", "),
         "\nAjusta la sección 1 de este script.")
  }
}
comprobar(fosfo, c(cols_zpq, cols_stat_fosfo, "protein", "gen", col_peptido), "fosfoproteoma")
comprobar(tabla, c("proteina", "nivel"), "tabla_completa")
comprobar(prot,  c(cols_zq, "protein", "gen"), "proteoma")

cat("\nValores de 'nivel' en tabla_completa:\n")
print(table(tabla$nivel, useNA = "ifany"))
if (!any(tabla$nivel %in% niveles_vih)) {
  stop("Ningún valor de 'nivel' coincide con NIVELES_VIH. Ajusta la sección 1.")
}


#3. Criterio (1): fosfopéptidos significativos
fosfo <- fosfo %>%
  mutate(
    sig_imatinib  = !is.na(Zpq_Ttest_imatinib)  & Zpq_Ttest_imatinib  < umbral_p &
                    !is.na(Zpq_AvA_imatinib)    & abs(Zpq_AvA_imatinib)  > umbral_ava,
    sig_dasatinib = !is.na(Zpq_Ttest_dasatinib) & Zpq_Ttest_dasatinib < umbral_p &
                    !is.na(Zpq_AvA_dasatinib)   & abs(Zpq_AvA_dasatinib) > umbral_ava
  )

prot_sig <- fosfo %>%
  group_by(protein, gen) %>%
  summarise(sig_imatinib  = any(sig_imatinib),
            sig_dasatinib = any(sig_dasatinib),
            n_fosfopeptidos = n(),
            n_sig = sum(sig_imatinib | sig_dasatinib),
            .groups = "drop") %>%
  filter(sig_imatinib | sig_dasatinib)

cat("\nProteínas con >=1 fosfopéptido significativo:", nrow(prot_sig), "\n")

#4. Criterio (2): relación documentada con VIH-1 (cruce NCBI)
nivel_prot <- tabla %>%
  filter(nivel %in% niveles_vih) %>%
  arrange(nivel) %>%
  distinct(proteina, .keep_all = TRUE) %>%
  select(protein = proteina, nivel_ncbi = nivel)

prot_vih <- prot_sig %>% inner_join(nivel_prot, by = "protein")
cat("Tras el filtro NCBI (", paste(niveles_vih, collapse = " + "), "):",
    nrow(prot_vih), "proteínas\n")

#5. Criterio (3): categorías funcionales (GO + complemento)
# 5.1 Anotaciones GO directas (org.Hs.eg.db) por símbolo génico
genes <- prot_vih %>% filter(!is.na(gen)) %>% pull(gen) %>% unique()
go_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
                        columns = c("GO", "ONTOLOGY"))
) %>% filter(!is.na(GO)) %>% distinct(SYMBOL, GO, ONTOLOGY)

go_nombres <- suppressMessages(
  AnnotationDbi::select(GO.db, keys = unique(go_map$GO), keytype = "GOID",
                        columns = "TERM")
) %>% rename(GO = GOID)

go_terms <- go_map %>% left_join(go_nombres, by = "GO")

# 5.2 Patrones por categoría sobre el nombre del término GO
#     (límites de palabra donde hay riesgo de falso positivo, p. ej. "actin")
#     "Señalización celular" se restringe a propósito: "signal transduction"
#     por sí sola capturaría cientos de proteínas.

patrones <- list(
  "Poro/transporte nuclear" = paste(
    "nuclear pore", "nuclear import", "nuclear export", "nucleocytoplasmic",
    "import into nucleus", "export from nucleus", "nuclear localization sequence",
    "nuclear transport", "karyopherin", sep = "|"),
  "Complejo de preintegración" = paste(
    "viral genome integration", "DNA integration", "provir", "preintegration",
    "pre-integration", "nonhomologous end joining", sep = "|"),
  "Spliceosoma/procesamiento ARN" = paste(
    "splic", "spliceosom", "mRNA processing", "mRNA 3'-end processing",
    "mRNA 5'-splice site", sep = "|"),
  "Citoesqueleto/tráfico intracelular" = paste(
    "cytoskelet", "\\bactin\\b", "microtubule", "intermediate filament",
    "kinesin", "dynein", "cell cortex", "lamellipodi", "stress fiber",
    "vesicle-mediated transport", "intracellular transport",
    "intracellular protein transport", "endosom", sep = "|"),
  "Señalización celular" = paste(
    "T cell receptor signaling", "T cell activation", "tyrosine kinase",
    "JAK-STAT", "cytokine-mediated signaling", "NF-kappaB", "MAPK cascade",
    "phosphatidylinositol 3-kinase", "ERK1 and ERK2", "interleukin.*signaling",
    "Fc receptor signaling", sep = "|")
)

go_hits <- purrr::imap(patrones, function(rx, cat) {
  d <- go_terms %>%
    filter(!is.na(TERM), str_detect(TERM, regex(rx, ignore_case = TRUE)))
  if (cat == "Poro/transporte nuclear") {
    d <- d %>% filter(!str_detect(TERM, regex("regulation of", ignore_case = TRUE)))
  }
  d %>% mutate(categoria = cat)
}) %>% bind_rows()

# 5.3 Complemento (UniProt / panel dirigido / hallazgo del cruce NCBI, script 04).

complemento <- tribble(
  ~gen,     ~categoria,                      ~fuente,
  "RANBP2", "Poro/transporte nuclear",       "Panel dirigido (script 03)",
  "NUP62",  "Poro/transporte nuclear",       "Panel dirigido (script 03)",
  "TNPO3",  "Poro/transporte nuclear",       "Panel dirigido (script 03)",
  "BANF1",  "Complejo de preintegración",    "Panel dirigido (script 03)",
  "HMGA1",  "Complejo de preintegración",    "Panel dirigido (script 03)",
  "PPIA",   "Complejo de preintegración",    "Panel dirigido (script 03)",
  "XRCC5",  "Complejo de preintegración",    "Panel dirigido (script 03)",
  "XRCC6",  "Complejo de preintegración",    "Panel dirigido (script 03)",
  "PSIP1",  "Complejo de preintegración",    "Panel dirigido (script 03)",
  "PRKDC",  "Complejo de preintegración",    "Anotación UniProt/GO (NHEJ)",
  "LCK",    "Señalización celular",          "Hallazgo cruce NCBI (script 04)",
  "NFKB1",  "Señalización celular",          "Hallazgo cruce NCBI (script 04)",
  "STAT1",  "Señalización celular",          "Hallazgo cruce NCBI (script 04)",
  "STAT3",  "Señalización celular",          "Hallazgo cruce NCBI (script 04)",
  "MAPK14", "Señalización celular",          "Hallazgo cruce NCBI (script 04)",
  "PIK3CD", "Señalización celular",          "Hallazgo cruce NCBI (script 04)"
)

# 5.4 Asignación completa gen -> categorías (con trazabilidad de la fuente)
asig_go <- go_hits %>%
  group_by(SYMBOL, categoria) %>%
  summarise(fuente = "GO", terminos = paste(head(unique(TERM), 5), collapse = "; "),
            .groups = "drop") %>%
  rename(gen = SYMBOL)

asig_comp <- complemento %>%
  filter(gen %in% genes) %>%
  mutate(terminos = NA_character_)

asignacion <- bind_rows(asig_go, asig_comp) %>%
  group_by(gen, categoria) %>%
  summarise(fuente = paste(unique(fuente), collapse = " + "),
            terminos = paste(na.omit(unique(terminos)), collapse = "; "),
            .groups = "drop")

write_csv(asignacion, "tablas/09_asignacion_categorias_heatmap.csv")
write_csv(go_hits %>% select(gen = SYMBOL, GO, ONTOLOGY, TERM, categoria),
          "tablas/09_go_terminos_por_categoria.csv")

# 5.5 Categoría principal (regla de prioridad) y categorías secundarias
cat_prin <- asignacion %>%
  mutate(rango = match(categoria, ORDEN_CAT)) %>%
  group_by(gen) %>%
  summarise(
    categoria_principal = categoria[which.min(rango)],
    categorias_todas = paste(categoria[order(rango)], collapse = " | "),
    n_categorias = n(),
    .groups = "drop"
  )

GENES_SENAL_T <- c("LCK", "NFKB1", "STAT1", "STAT3", "MAPK14", "PIK3CD")

cat_prin <- cat_prin %>%
  mutate(categoria_principal = ifelse(
    gen %in% GENES_SENAL_T & str_detect(categorias_todas, fixed("Señalización")),
    "Señalización celular", categoria_principal))

seleccion <- prot_vih %>%
  inner_join(cat_prin, by = "gen") %>%
  mutate(categoria_principal = factor(categoria_principal, levels = ORDEN_CAT)) %>%
  arrange(categoria_principal, gen)

# Una misma proteína (entry name) no puede aparecer dos veces en un heatmap
if (anyDuplicated(seleccion$protein) > 0) {
  cat("\nAviso: entradas duplicadas por gen compartido entre UniProt; se conserva la primera.\n")
  print(seleccion %>% filter(duplicated(protein) | duplicated(protein, fromLast = TRUE)) %>%
          select(protein, gen))
  seleccion <- seleccion %>% distinct(protein, .keep_all = TRUE)
}

write_csv(seleccion, "tablas/09_seleccion_proteinas_heatmap.csv")
write_csv(seleccion %>% filter(n_categorias > 1),
          "tablas/09_proteinas_multicategoria.csv")
saveRDS(seleccion, "data_processed/09_seleccion_heatmap_vih.rds")

cat("\nProteínas seleccionadas:", nrow(seleccion), "\n")
print(table(seleccion$categoria_principal))
cat("Proteínas en más de una categoría:", sum(seleccion$n_categorias > 1), "\n")

# Control: el bloque de linfocitos T debería aparecer si cumple (1) y (2)
ctrl <- c("LCK", "NFKB1", "STAT1", "STAT3", "MAPK14", "PIK3CD")
cat("Control señalización T (presentes en la selección):",
    paste(intersect(ctrl, seleccion$gen), collapse = ", "), "\n")


#6. Matriz Zpq: una fila por fosfopéptido
# Se incluyen solo los fosfopéptidos significativos en al menos un fármaco
# (solo_sig = TRUE); con FALSE se incluirían todos los de las proteínas seleccionadas.
fp <- fosfo %>%
  { if (solo_sig) filter(., sig_imatinib | sig_dasatinib) else . } %>%
  inner_join(seleccion %>% select(protein, categoria_principal, nivel_ncbi),
             by = "protein") %>%
  mutate(etiqueta = paste0(gen, " | ", str_remove(.data[[col_peptido]], "_Phospho")))

m_zpq <- as.matrix(fp[, cols_zpq]); colnames(m_zpq) <- muestras
rownames(m_zpq) <- make.unique(fp$etiqueta, sep = " #")

completas <- complete.cases(m_zpq)
cat("\nFosfopéptidos con Zpq:", nrow(fp), "| descartados por NA en alguna muestra:",
    sum(!completas), "\n")
  fp <- fp[completas, ]; m_zpq <- m_zpq[completas, , drop = FALSE]

zscore_filas <- function(m) {
  z <- t(scale(t(m))); z[is.nan(z)] <- NA
  z[rowSums(!is.na(z)) > 0, , drop = FALSE]
}
z_zpq <- zscore_filas(m_zpq)
fp <- fp[match(rownames(z_zpq), rownames(m_zpq)), ]


#7. Matriz Zq: una fila por proteína
zq_unico <- prot %>%
  mutate(.k = do.call(paste, c(across(all_of(cols_zq)), sep = "|"))) %>%
  group_by(protein, .k) %>% mutate(.n = n()) %>% ungroup() %>%
  arrange(protein, desc(.n)) %>%      # primero el Zq más frecuente
  distinct(protein, .keep_all = TRUE) %>%
  select(-.k, -.n)

cols_stat_zq <- c("Zq_AvA_imatinib", "Zq_Ttest_imatinib", "Zq_AvA_dasatinib", "Zq_Ttest_dasatinib")
tiene_stat_zq <- all(cols_stat_zq %in% names(zq_unico))
if (!tiene_stat_zq) {
  warning("No se encuentran columnas de estadísticos de Zq (", paste(cols_stat_zq, collapse = ", "),
          "); no se anotará la significación en Zq.")
}

zq_sel <- seleccion %>%
  select(protein, gen, categoria_principal, nivel_ncbi, sig_imatinib, sig_dasatinib) %>%
  inner_join(zq_unico %>% select(protein, any_of(c(cols_zq, cols_stat_zq))), by = "protein")

sin_zq <- seleccion %>% filter(!protein %in% zq_sel$protein)
write_csv(sin_zq %>% select(protein, gen, categoria_principal), "tablas/09_proteinas_sin_Zq.csv")
cat("\nProteínas seleccionadas sin Zq:", nrow(sin_zq), "\n")

if (tiene_stat_zq) {
  zq_sel <- zq_sel %>% mutate(
    zq_sig_imatinib  = !is.na(Zq_Ttest_imatinib)  & Zq_Ttest_imatinib  < umbral_p & abs(Zq_AvA_imatinib)  > umbral_ava,
    zq_sig_dasatinib = !is.na(Zq_Ttest_dasatinib) & Zq_Ttest_dasatinib < umbral_p & abs(Zq_AvA_dasatinib) > umbral_ava
  )
}

m_zq <- as.matrix(zq_sel[, cols_zq]); colnames(m_zq) <- muestras
etq <- ifelse(duplicated(zq_sel$gen) | duplicated(zq_sel$gen, fromLast = TRUE),
              paste0(zq_sel$gen, " (", zq_sel$protein, ")"), zq_sel$gen)
rownames(m_zq) <- etq
ok_zq <- complete.cases(m_zq)
cat("Proteínas con Zq descartadas por NA:", sum(!ok_zq), "\n")
zq_sel <- zq_sel[ok_zq, ]; m_zq <- m_zq[ok_zq, , drop = FALSE]
z_zq <- zscore_filas(m_zq)
zq_sel <- zq_sel[match(rownames(z_zq), rownames(m_zq)), ]


#8. Función de dibujo
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
col_cond <- c(Control = "#7F7F7F", Imatinib = "#E69F00", Dasatinib = "#009E73")
col_si   <- c("Sí" = "#CC0000", "No" = "grey90")
col_nivel <- setNames(c("#6A3D9A", "#B2ABD2")[seq_along(niveles_vih)], niveles_vih)

dibujar_heatmap <- function(z, meta, anot, titulo, archivo, tam_fuente = 7) {
  ha_col <- HeatmapAnnotation(Condición = condicion, col = list(Condición = col_cond),
                              annotation_name_side = "left")
  ha_row <- rowAnnotation(df = anot, col = col_list_anot(anot), border = TRUE,
                          annotation_name_gp = gpar(fontsize = 8))
  ht <- Heatmap(
    z, name = "Z-score\n(por fila)", col = col_fun,
    top_annotation = ha_col, right_annotation = ha_row,
    row_split = meta$categoria_principal, row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"), row_gap = unit(3, "mm"),
    column_split = condicion, cluster_columns = FALSE, column_gap = unit(2, "mm"),
    column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = tam_fuente),
    cluster_rows = TRUE, show_row_dend = FALSE, border = TRUE,
    column_title = titulo
  )
  alto <- max(6, 0.17 * nrow(z) + 3)
  for (ext in c("pdf", "png")) {
    f <- paste0("figuras/", archivo, ".", ext)
    if (ext == "pdf") pdf(f, width = 11, height = alto) else png(f, width = 11, height = alto, units = "in", res = 300)
    draw(ht, merge_legend = TRUE, heatmap_legend_side = "right",
         annotation_legend_side = "right")
    dev.off()
  }
  invisible(ht)
}

col_list_anot <- function(anot) {
  lapply(anot, function(x) if (all(unique(x) %in% names(col_nivel))) col_nivel else col_si)
}

por_categoria <- function(z, meta, anot, titulo, prefijo, tam_fuente) {
  for (k in seq_along(ORDEN_CAT)) {
    idx <- which(meta$categoria_principal == ORDEN_CAT[k])
    if (length(idx) < 3) next          # hclust necesita >= 3 filas
    m <- meta[idx, ]
    m$categoria_principal <- droplevels(m$categoria_principal)
    nombre <- paste0(prefijo, "_", k, "_", gsub("[^A-Za-z]+", "_", ORDEN_CAT[k]))
    dibujar_heatmap(z[idx, , drop = FALSE], m, anot[idx, , drop = FALSE],
                    paste0(titulo, ": ", ORDEN_CAT[k]), nombre, tam_fuente)
  }
}

#9. Heatmap Zpq (fosfoproteoma) 
anot_zpq <- data.frame(
  `Sig. imatinib`  = ifelse(fp$sig_imatinib,  "Sí", "No"),
  `Sig. dasatinib` = ifelse(fp$sig_dasatinib, "Sí", "No"),
  `NCBI VIH-1`     = fp$nivel_ncbi,
  check.names = FALSE
)
por_categoria(z_zpq, fp, anot_zpq,
              "Fosfoproteoma (Zpq)", "heatmap_Zpq", tam_fuente = 6)


#10. Heatmap Zq (proteoma)
anot_zq <- data.frame(
  `Fosfo sig. imatinib`  = ifelse(zq_sel$sig_imatinib,  "Sí", "No"),
  `Fosfo sig. dasatinib` = ifelse(zq_sel$sig_dasatinib, "Sí", "No"),
  check.names = FALSE
)
if (tiene_stat_zq) {
  anot_zq[["Zq sig. imatinib"]]  <- ifelse(zq_sel$zq_sig_imatinib,  "Sí", "No")
  anot_zq[["Zq sig. dasatinib"]] <- ifelse(zq_sel$zq_sig_dasatinib, "Sí", "No")
}
anot_zq[["NCBI VIH-1"]] <- zq_sel$nivel_ncbi
por_categoria(z_zq, zq_sel, anot_zq,
              "Proteoma (Zq)", "heatmap_Zq", tam_fuente = 7)


#11. Resumen para Resultados / Métodos
resumen <- tibble(
  proteinas_sig_fosfo      = nrow(prot_sig),
  tras_filtro_NCBI         = nrow(prot_vih),
  en_categorias_funcionales = nrow(seleccion),
  multicategoria           = sum(seleccion$n_categorias > 1),
  filas_heatmap_Zq         = nrow(z_zq),
  filas_heatmap_Zpq        = nrow(z_zpq),
  proteinas_sin_Zq         = nrow(sin_zq)
)
print(resumen)
write_csv(resumen, "tablas/09_resumen_heatmaps_vih.csv")
