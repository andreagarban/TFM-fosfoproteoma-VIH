# ==============================================================================
# 05_cruce_bases_datos.R
#
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con inhibidores
#      de tirosina quinasa como estrategia antiviral frente al VIH

# Objetivo: contrastar el listado de proteínas obtenido tras el filtrado con
#           las interacciones VIH-1 / hospedador documentadas en la HIV-1 Human
#           Interaction Database (NCBI), con el fin de distinguir entre
#           proteínas con relación descrita con el virus y proteínas
#           seleccionadas únicamente por su categoría funcional.
#
# Fuente de datos:
#   hiv_interactions.gz       Interacciones moleculares proteína-proteína
#                             descritas en la literatura
#   hiv_siRNA_interactions.gz Genes humanos cuyo silenciamiento mediante siRNA
#                             modifica la replicación viral
# ==============================================================================

library(tidyverse)
library(org.Hs.eg.db)

select <- dplyr::select # para evitar conflictos en el resto del script
rename <- dplyr::rename
filter <- dplyr::filter

fosfo_imatinib      <- readRDS("data_processed/fosfo_imatinib.rds")
fosfo_dasatinib     <- readRDS("data_processed/fosfo_dasatinib.rds")
fosfo_imatinib_VIH  <- readRDS("data_processed/fosfo_imatinib_VIH.rds")
fosfo_dasatinib_VIH <- readRDS("data_processed/fosfo_dasatinib_VIH.rds")
busqueda_VIH        <- readRDS("data_processed/busqueda_VIH.rds")


#1. LECTURA DE LOS FICHEROS DE INTERACCIÓN

# Primeras líneas del fichero de interacciones
readLines("data_raw/hiv_interactions.gz", n = 5)
readLines("data_raw/hiv_siRNA_interactions.gz", n = 5)


# 1.1. Interacciones proteína-proteína
# Se asignan los nombres de columna de forma explícita.

nombres_ppi <- c("tax_id_1", "gene_id_1", "accesion_1", "producto_1",
                 "interaccion",
                 "tax_id_2", "gene_id_2", "accesion_2", "producto_2",
                 "pmid", "actualizacion", "generif")

hiv_ppi <- read_tsv("data_raw/hiv_interactions.gz",
                    skip = 1,
                    col_names = nombres_ppi,
                    col_types = cols(.default = col_character()),
                    quote = "")

nrow(hiv_ppi) # Interacciones proteína-proteína leídas


# 1.2. Interacciones obtenidas mediante silenciamiento con siRNA

nombres_sirna <- c("tax_id", "gene_id", "interaccion", "accesion",
                   "producto", "pmid", "actualizacion", "generif")

hiv_sirna <- read_tsv("data_raw/hiv_siRNA_interactions.gz",
                      skip = 1,
                      col_names = nombres_sirna,
                      col_types = cols(.default = col_character()),
                      quote = "")

nrow(hiv_sirna) # Interacciones por siRNA leídas


# 2. EXTRACCIÓN DEL COMPONENTE HUMANO DE CADA INTERACCIÓN
# En el fichero de interacciones proteína-proteína, cada registro enfrenta una proteína viral (taxid 11676) con una proteína humana (taxid 9606). 
# El orden de ambos no está garantizado, por lo que se identifica el lado humano en función del identificador de taxón en lugar de asumir una posición fija.

TAXID_HUMANO <- "9606"

ppi_humano <- hiv_ppi %>%
  mutate(
    gene_id_humano = case_when(
      tax_id_1 == TAXID_HUMANO ~ gene_id_1,
      tax_id_2 == TAXID_HUMANO ~ gene_id_2,
      TRUE                     ~ NA_character_
    ),
    proteina_viral = case_when(
      tax_id_1 == TAXID_HUMANO ~ producto_2,
      tax_id_2 == TAXID_HUMANO ~ producto_1,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(gene_id_humano)) %>%
  select(gene_id_humano, proteina_viral, interaccion, pmid, generif)

nrow(ppi_humano) # Registros con componente humano identificado


#3. CONVERSIÓN DE IDENTIFICADORES
# Se convierten los GeneID de NCBI a símbolos génicos (identificador común con el conjunto de datos experimental)

geneids <- unique(c(ppi_humano$gene_id_humano, hiv_sirna$gene_id))
geneids <- geneids[!is.na(geneids)]
length(geneids)

mapeo <- AnnotationDbi::select(org.Hs.eg.db,
                               keys    = geneids,
                               keytype = "ENTREZID",
                               columns = "SYMBOL") %>%
  as_tibble() %>%
  filter(!is.na(SYMBOL)) %>%
  rename(gene_id = ENTREZID, gen = SYMBOL)

#Identificadores convertidos:
nrow(mapeo) # de 4662


#4. Construcción de las tablas de referencia por gen

# 4.1. Interacciones proteína-proteína agregadas por gen humano

ref_ppi <- ppi_humano %>%
  left_join(mapeo, by = c("gene_id_humano" = "gene_id")) %>%
  filter(!is.na(gen)) %>%
  group_by(gen) %>%
  summarise(
    n_interacciones   = n(),
    proteinas_virales = paste(sort(unique(proteina_viral)), collapse = "; "),
    tipos_interaccion = paste(sort(unique(interaccion)),   collapse = "; "),
    n_publicaciones   = n_distinct(unlist(str_split(pmid, ","))),
    .groups = "drop"
  )

#Genes humanos con interacción documentada
nrow(ref_ppi)


# 4.2. Genes identificados en cribados con siRNA

ref_sirna <- hiv_sirna %>%
  left_join(mapeo, by = "gene_id") %>%
  filter(!is.na(gen)) %>%
  group_by(gen) %>%
  summarise(
    n_evidencias_sirna = n(),
    efecto_sirna       = paste(sort(unique(interaccion)), collapse = "; "),
    .groups = "drop"
  )

#Genes humanos con evidencia por siRNA
nrow(ref_sirna)


#5. Cruce con el conjunto de datos experimental

cruzar_con_bd <- function(datos, columna_AvA, columna_p, columna_fdr, etiqueta) {
  datos %>%
    group_by(uniprot = fasta, gen, proteina = protein, descripcion = description) %>%
    summarise(
      n_fosfopeptidos = n(),
      n_aumentan      = sum(direccion == "Aumenta"),
      n_disminuyen    = sum(direccion == "Disminuye"),
      AvA_max         = .data[[columna_AvA]][which.max(abs(.data[[columna_AvA]]))],
      p_valor_min     = min(.data[[columna_p]]),
      FDR_min         = min(.data[[columna_fdr]]),
      n_fdr10         = sum(.data[[columna_fdr]] < 0.10),
      .groups = "drop"
    ) %>%
    left_join(ref_ppi,   by = "gen") %>%
    left_join(ref_sirna, by = "gen") %>%
    mutate(
      interaccion_directa = !is.na(n_interacciones),
      evidencia_sirna     = !is.na(n_evidencias_sirna),
      tratamiento         = etiqueta
    ) %>%
    arrange(desc(interaccion_directa), desc(abs(AvA_max)))
}

cruce_imatinib <- cruzar_con_bd(fosfo_imatinib,
                                "Zpq_AvA_imatinib",
                                "Zpq_Ttest_imatinib",
                                "Imatinib_FDR",
                                "Imatinib")

cruce_dasatinib <- cruzar_con_bd(fosfo_dasatinib,
                                 "Zpq_AvA_dasatinib",
                                 "Zpq_Ttest_dasatinib",
                                 "Dasatinib_FDR",
                                 "Dasatinib")


#6. Resumen de la cobertura del cruce

resumen_cruce <- bind_rows(cruce_imatinib, cruce_dasatinib) %>%
  group_by(tratamiento) %>%
  summarise(
    proteinas_totales   = n(),
    con_interaccion_ppi = sum(interaccion_directa),
    con_evidencia_sirna = sum(evidencia_sirna),
    con_alguna_evidencia = sum(interaccion_directa | evidencia_sirna),
    .groups = "drop"
  ) %>%
  mutate(
    pct_ppi    = round(100 * con_interaccion_ppi / proteinas_totales, 1),
    pct_sirna  = round(100 * con_evidencia_sirna / proteinas_totales, 1),
    pct_alguna = round(100 * con_alguna_evidencia / proteinas_totales, 1)
  )

#Cobertura del cruce con la base de datos
print(resumen_cruce)


#7. Clasificación en niveles de justificación
# Las proteínas se clasifican en tres niveles según la solidez de su vinculación con el VIH-1:
#
#   Nivel 1. Relación directa: interacción proteína-proteína documentada con
#            alguna proteína viral.
#   Nivel 2. Relación indirecta: sin interacción directa documentada, pero con
#            evidencia funcional procedente de cribados con siRNA, o
#            pertenencia a una categoría funcional implicada en el ciclo viral.
#   Nivel 3. Candidatas exploratorias: seleccionadas únicamente por el
#            resultado del filtrado, sin respaldo específico en la literatura
#            de VIH-1.

# Categorías funcionales de interés definidas en el diseño del trabajo
keywords_categorias <- c(
  "nucleoporin", "nuclear pore", "importin", "karyopherin", "transportin",
  "actin", "tubulin", "myosin", "vimentin", "filamin", "cofilin", "profilin",
  "talin", "vinculin", "splicing", "spliceosom", "snRNP",
  "heterogeneous nuclear", "hnRNP", "serine/arginine"
)

patron_categorias <- paste(keywords_categorias, collapse = "|")

clasificar_niveles <- function(cruce) {
  cruce %>%
    mutate(
      en_categoria = str_detect(descripcion,
                                regex(patron_categorias, ignore_case = TRUE)),
      nivel = case_when(
        interaccion_directa            ~ "1. Relación directa con VIH-1",
        evidencia_sirna | en_categoria ~ "2. Relación indirecta o funcional",
        TRUE                           ~ "3. Candidata exploratoria"
      )
    )
}

cruce_imatinib  <- clasificar_niveles(cruce_imatinib)
cruce_dasatinib <- clasificar_niveles(cruce_dasatinib)

#Distribución por nivel de justificación
bind_rows(cruce_imatinib, cruce_dasatinib) %>%
  count(tratamiento, nivel) %>%
  group_by(tratamiento) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  print()


#8. Contraste con el panel de búsqueda dirigida
# Se comprueba qué proteínas del panel definido en el script anterior cuentan
# con interacción documentada, con independencia de que hayan superado los criterios de filtrado

panel_contrastado <- busqueda_VIH %>%
  distinct(gen, categoria) %>%
  left_join(ref_ppi,   by = "gen") %>%
  left_join(ref_sirna, by = "gen") %>%
  mutate(
    interaccion_directa = !is.na(n_interacciones),
    evidencia_sirna     = !is.na(n_evidencias_sirna)
  ) %>%
  select(gen, categoria, interaccion_directa, n_interacciones,
         proteinas_virales, evidencia_sirna, n_evidencias_sirna)

#Panel dirigido contrastado con la base de datos
print(panel_contrastado, n = Inf, width = Inf)


#9. Validación de la selección por palabras clave
#Se evalúa qué proporción de las proteínas así seleccionadas cuenta con respaldo en la base de datos de
#interacciones, y qué proteínas con interacción documentada quedaron fuera de dicha selección.

# 9.1. Cobertura de la lista obtenida por palabras clave
validar_keywords <- function(datos_keywords, cruce_completo, etiqueta) {
  genes_keywords <- datos_keywords %>%
    distinct(uniprot = fasta, gen, proteina = protein, descripcion = description)

  genes_keywords %>%
    left_join(cruce_completo %>%
                select(uniprot, interaccion_directa, evidencia_sirna,
                       n_interacciones, proteinas_virales, nivel),
              by = "uniprot") %>%
    mutate(tratamiento = etiqueta,
           interaccion_directa = replace_na(interaccion_directa, FALSE),
           evidencia_sirna     = replace_na(evidencia_sirna, FALSE))
}

keywords_validado_imatinib <- validar_keywords(fosfo_imatinib_VIH,
                                               cruce_imatinib,
                                               "Imatinib")

keywords_validado_dasatinib <- validar_keywords(fosfo_dasatinib_VIH,
                                                cruce_dasatinib,
                                                "Dasatinib")

resumen_keywords <- bind_rows(keywords_validado_imatinib,
                              keywords_validado_dasatinib) %>%
  group_by(tratamiento) %>%
  summarise(
    seleccionadas_keywords = n(),
    con_interaccion_ppi    = sum(interaccion_directa),
    con_evidencia_sirna    = sum(evidencia_sirna),
    sin_respaldo           = sum(!interaccion_directa & !evidencia_sirna),
    .groups = "drop"
  ) %>%
  mutate(
    pct_con_respaldo = round(100 * (seleccionadas_keywords - sin_respaldo) /
                               seleccionadas_keywords, 1)
  )

#Validación de la selección por palabras clave
print(resumen_keywords)


#9.2. Proteínas con interacción documentada omitidas por las palabras clave
# Identifica las proteínas del conjunto filtrado que presentan interacción
# descrita con el VIH-1 y que, sin embargo, no fueron recuperadas por la
# búsqueda de términos

omitir_keywords <- function(cruce, keywords_validado) {
  cruce %>%
    filter(interaccion_directa,
           !uniprot %in% keywords_validado$uniprot) %>%
    select(uniprot, gen, proteina, descripcion, n_fosfopeptidos,
           AvA_max, p_valor_min, n_interacciones, proteinas_virales) %>%
    arrange(desc(n_interacciones))
}

omitidas_por_keywords_imatinib  <- omitir_keywords(cruce_imatinib,
                                                   keywords_validado_imatinib)
omitidas_por_keywords_dasatinib <- omitir_keywords(cruce_dasatinib,
                                                   keywords_validado_dasatinib)

# Resumen comparable entre fármacos
resumen_omitidas <- tibble(
  tratamiento     = c("Imatinib", "Dasatinib"),
  con_interaccion = c(sum(cruce_imatinib$interaccion_directa),
                      sum(cruce_dasatinib$interaccion_directa)),
  omitidas        = c(nrow(omitidas_por_keywords_imatinib),
                      nrow(omitidas_por_keywords_dasatinib))
) %>%
  mutate(pct_omitidas = round(100 * omitidas / con_interaccion, 1))

print(resumen_omitidas)
print(head(omitidas_por_keywords_imatinib, 20),  width = Inf)
print(head(omitidas_por_keywords_dasatinib, 20), width = Inf)

#9.3. Bloque de señalización de linfocitos T identificado a partir del cruce

bloque_senalizacion <- c("LCK", "NFKB1", "STAT1", "STAT3", "MAPK14", "PIK3CD")

keywords_recuperadas <- bind_rows(keywords_validado_imatinib,
                                  keywords_validado_dasatinib) %>%
  distinct(tratamiento, uniprot) %>%
  mutate(recuperada_por_keywords = TRUE)

bloque_cruce <- bind_rows(cruce_imatinib, cruce_dasatinib) %>%
  filter(gen %in% bloque_senalizacion) %>%
  left_join(keywords_recuperadas, by = c("tratamiento", "uniprot")) %>%
  mutate(
    recuperada_por_keywords = replace_na(recuperada_por_keywords, FALSE),
    categoria_hallazgo = "Señalización de linfocitos T (hallazgo del cruce)"
  ) %>%
  select(tratamiento, gen, uniprot, proteina, descripcion, n_fosfopeptidos,
         AvA_max, p_valor_min, FDR_min, interaccion_directa, n_interacciones,
         proteinas_virales, evidencia_sirna, nivel, recuperada_por_keywords,
         categoria_hallazgo) %>%
  arrange(gen, tratamiento)

#Bloque de señalización de linfocitos T en el conjunto filtrado
print(bloque_cruce, n = Inf, width = Inf)

#Genes del bloque ausentes en cada comparación (no superan los filtros)
bloque_senalizacion %>%
  setdiff(bloque_cruce$gen[bloque_cruce$tratamiento == "Imatinib"])
bloque_senalizacion %>%
  setdiff(bloque_cruce$gen[bloque_cruce$tratamiento == "Dasatinib"])

#Genes con más de una entrada UniProt por tratamiento
bloque_cruce %>% count(tratamiento, gen) %>% filter(n > 1)


# 10. GUARDAR

saveRDS(ref_ppi,           "data_processed/ref_ppi_hiv.rds")
saveRDS(ref_sirna,         "data_processed/ref_sirna_hiv.rds")
saveRDS(cruce_imatinib,    "data_processed/cruce_imatinib.rds")
saveRDS(cruce_dasatinib,   "data_processed/cruce_dasatinib.rds")
saveRDS(resumen_cruce,     "data_processed/resumen_cruce.rds")
saveRDS(panel_contrastado, "data_processed/panel_contrastado.rds")
saveRDS(keywords_validado_imatinib,  "data_processed/keywords_validado_imatinib.rds")
saveRDS(keywords_validado_dasatinib, "data_processed/keywords_validado_dasatinib.rds")
saveRDS(resumen_keywords,       "data_processed/resumen_keywords.rds")
saveRDS(omitidas_por_keywords_imatinib,  "data_processed/omitidas_por_keywords_imatinib.rds")
saveRDS(omitidas_por_keywords_dasatinib, "data_processed/omitidas_por_keywords_dasatinib.rds")
saveRDS(resumen_omitidas,                "data_processed/resumen_omitidas.rds")
saveRDS(bloque_cruce, "data_processed/bloque_senalizacion.rds")
