# =============================================================================
# 14_string_red_global.R
#
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con inhibidores
#      de tirosina quinasa como estrategia antiviral frente al VIH
#
# Objetivo : Construir una red STRING GLOBAL por fármaco (imatinib y dasatinib)
#            con todos los genes asociados a los fosfopéptidos que cumplen los
#            criterios de selección (p < 0.05 y |AvA| > 1.3), SIN restringir
#            STRING a las proteínas de la base de datos de interacciones con VIH.
# =============================================================================

#0. Paquetes y alias

library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(STRINGdb)
library(igraph)


select    <- dplyr::select
rename    <- dplyr::rename
filter    <- dplyr::filter
slice     <- dplyr::slice
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
arrange   <- dplyr::arrange
count     <- dplyr::count

# Criterios de selección
p_umbral   <- 0.05
ava_umbral <- 1.3

# Parámetros de STRING

string_version <- "12.0"
especie        <- 9606            # Homo sapiens
tipo_red       <- "full"          # red funcional + física (opción por defecto)
umbrales       <- c(400, 700)     # confianza media (0.4) y alta (0.7), escala 0-1000

dir.create("data_raw/string_cache", recursive = TRUE, showWarnings = FALSE)

# Recuentos esperados según el pipeline ya verificado (solo para comprobar)
ESPERADO_FOSFOPEPTIDOS <- c(imatinib = 416, dasatinib = 1926)

farmacos <- c("imatinib", "dasatinib")


# 1. Carga de los fosfopéptidos seleccionados

fosfo <- list(
  imatinib  = readRDS("data_processed/fosfo_imatinib.rds"),
  dasatinib = readRDS("data_processed/fosfo_dasatinib.rds")
)


#2. Entrada a nivel de gen (una fila por gen) 
# Varios fosfopéptidos de una misma proteína colapsan en un único gen. Se
# guarda el número de fosfopéptidos por gen (n_fosfopeptidos) para la revisión
# posterior de proteínas con varios fosfopéptidos (no se guarda. la direccion up/down)

entrada_gen   <- list()
sin_simbolo   <- list()

# Tabla puente: nombre de entrada UniProt -> accession (tabla_completa, script 07)
tabla_completa <- readRDS("data_processed/tabla_uniprot_niveles.rds")
acc_lookup <- tabla_completa %>% distinct(proteina, uniprot)
stopifnot(!anyDuplicated(acc_lookup$proteina))

for (f in farmacos) {
  df <- fosfo[[f]] %>%
    mutate(gen = na_if(str_trim(as.character(gen)), ""))

  n_antes <- nrow(df)
  df <- df %>% left_join(acc_lookup, by = c("protein" = "proteina"))
  stopifnot(nrow(df) == n_antes)   # el cruce no debe duplicar filas
  message(sprintf("[%s] fosfopéptidos sin accession tras el cruce: %d", f, sum(is.na(df$uniprot))))

  sin_simbolo[[f]] <- df %>%
    filter(is.na(gen)) %>%
    distinct(protein, uniprot)

  entrada_gen[[f]] <- df %>%
    filter(!is.na(gen)) %>%
    group_by(gen) %>%
    summarise(
      n_fosfopeptidos = n(),
      n_proteinas     = n_distinct(protein),
      proteinas       = paste(sort(unique(protein)), collapse = ";"),
      uniprot         = if (all(is.na(uniprot))) NA_character_
                        else paste(sort(unique(na.omit(uniprot))), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(gen)

  message(sprintf("[%s] %d fosfopéptidos -> %d genes con símbolo (%d proteínas sin símbolo génico excluidas)",
                  f, nrow(fosfo[[f]]), nrow(entrada_gen[[f]]), nrow(sin_simbolo[[f]])))

  saveRDS(entrada_gen[[f]], paste0("data_processed/string_global_entrada_", f, ".rds"))
  writeLines(entrada_gen[[f]]$gen,
             paste0("tablas/14string_global_genes_entrada_", f, ".txt"))
}

message(sprintf("Genes compartidos entre ambas listas de entrada: %d",
                length(intersect(entrada_gen$imatinib$gen, entrada_gen$dasatinib$gen))))

#3. Mapeo de símbolos génicos a identificadores STRING
# Se mapea por símbolo. Los genes sin identificador se reintentan actualizando
# el símbolo a partir del código UniProt (org.Hs.eg.db) y todo caso revisado
# queda registrado en un fichero de log.

string_db <- STRINGdb$new(
  version         = string_version,
  species         = especie,
  score_threshold = min(umbrales),
  network_type    = tipo_red,
  input_directory = "data_raw/string_cache"
)

prot <- string_db$get_proteins()

# Correcciones manuales de símbolos (símbolo del dataset -> símbolo STRING/HGNC vigente).
# Todas quedan registradas en el log de mapeo (columna corregido = TRUE).
correcciones <- c(
  GARS          = "GARS1",   # el alias "GARS" caía en GART
  MARS          = "MARS1",   # el alias "MARS" caía en SLA2 (imatinib)
  HN1           = "JPT1",    # el alias "HN1" caía en MTRNR2L1
  HIST1H4A      = "H4C1",    # unifica la histona H4 entre fármacos (verificar en HGNC)
  DKFZP686J1372 = "TPM3",     # vía accession Q5HYB6 (alias UniProt de STRING)
  FAM21A        = "WASHC2A"  # el alias "FAM21A" caía en WASHC2C (parálogo)
)

mapear_string <- function(genes) {
  genes <- genes %>% mutate(gen_mapeo = coalesce(unname(correcciones[toupper(gen)]), gen))
  
  m <- string_db$map(as.data.frame(genes), "gen_mapeo",
                     removeUnmappedRows = FALSE, takeFirst = FALSE) %>%
    as_tibble() %>%
    left_join(select(prot, STRING_id = protein_external_id, nombre_string = preferred_name),
              by = "STRING_id") %>%
    mutate(exacto     = toupper(gen_mapeo) == toupper(nombre_string),
           gen_string = gen_mapeo,
           rescatado  = FALSE,
           corregido  = toupper(gen_mapeo) != toupper(gen))

  # (a) genes con varios identificadores STRING: se conserva el primero
  multi <- m %>% filter(!is.na(STRING_id)) %>% count(gen) %>% filter(n > 1) %>% pull(gen)
  log_multi <- m %>% filter(gen %in% multi) %>% mutate(incidencia = "varios_ids_STRING")
  m <- m %>% arrange(is.na(STRING_id), desc(exacto)) %>% distinct(gen, .keep_all = TRUE)

  # (b) genes sin identificador: reintento con el símbolo vigente según UniProt
  log_rescate <- tibble()
  sin_id <- m %>% filter(is.na(STRING_id), !is.na(uniprot))

  if (nrow(sin_id) > 0) {
    cand <- sin_id %>%
      select(gen, uniprot) %>%
      separate_rows(uniprot, sep = ";") %>%
      mutate(uniprot = str_remove(uniprot, "-\\d+$")) %>%
      mutate(gen_nuevo = suppressWarnings(suppressMessages(
        unname(mapIds(org.Hs.eg.db, keys = uniprot, keytype = "UNIPROT",
                      column = "SYMBOL", multiVals = "first"))))) %>%
      filter(!is.na(gen_nuevo), gen_nuevo != gen) %>%
      distinct(gen, gen_nuevo)

    if (nrow(cand) > 0) {
      m2 <- string_db$map(as.data.frame(cand), "gen_nuevo", removeUnmappedRows = TRUE) %>%
        as_tibble() %>%
        distinct(gen, .keep_all = TRUE) %>%
        select(gen, STRING_id_rescate = STRING_id, gen_nuevo)

      m <- m %>%
        left_join(m2, by = "gen") %>%
        mutate(
          rescatado  = is.na(STRING_id) & !is.na(STRING_id_rescate),
          STRING_id  = coalesce(STRING_id, STRING_id_rescate),
          gen_string = if_else(rescatado, gen_nuevo, gen_string)
        ) %>%
        select(-STRING_id_rescate, -gen_nuevo)

      log_rescate <- m %>% filter(rescatado) %>% mutate(incidencia = "rescatado_por_UniProt")
    }
  }

  # (c) genes que siguen sin identificador
  log_nomap <- m %>% filter(is.na(STRING_id)) %>% mutate(incidencia = "sin_ID_STRING")

  list(mapeo = m, log = bind_rows(log_multi, log_rescate, log_nomap))
}

mapeo <- list()
nodos <- list()

for (f in farmacos) {
  res <- mapear_string(entrada_gen[[f]])
  mapeo[[f]] <- res$mapeo

  write_csv(res$log, paste0("tablas/14string_global_log_mapeo_", f, ".csv"))

  # Un nodo por identificador STRING (si dos símbolos convergen en el mismo
  # identificador, se fusionan y el nodo conserva ambos símbolos)
  nodos[[f]] <- res$mapeo %>%
    filter(!is.na(STRING_id)) %>%
    group_by(STRING_id) %>%
    summarise(
      gen             = paste(unique(gen), collapse = "/"),
      gen_string      = paste(unique(gen_string), collapse = "/"),
      n_fosfopeptidos = sum(n_fosfopeptidos),
      .groups = "drop"
    ) %>%
    relocate(STRING_id)

  message(sprintf("[%s] genes de entrada: %d | mapeados: %d (rescatados por UniProt: %d) | sin ID: %d | nodos únicos: %d",
                  f, nrow(res$mapeo), sum(!is.na(res$mapeo$STRING_id)),
                  sum(res$mapeo$rescatado), sum(is.na(res$mapeo$STRING_id)),
                  nrow(nodos[[f]])))
}

#3b. Auditoría del mapeo
# STRING no avisa cuando un símbolo antiguo cae en otra proteína (p. ej. GARS -> GART).
# Se listan los genes cuyo símbolo enviado a STRING no coincide con el nombre preferido
# de STRING: los renombres legítimos de HGNC se aceptan; el resto hay que corregirlo.

auditoria <- bind_rows(mapeo, .id = "farmaco") %>%
  filter(!is.na(STRING_id)) %>%
  filter(toupper(gen_string) != toupper(nombre_string)) %>%
  select(farmaco, gen, gen_string, nombre_string, STRING_id, proteinas, corregido)

write_csv(auditoria, "tablas/14string_global_auditoria_mapeo.csv")
message(sprintf("Auditoría del mapeo: %d genes con nombre distinto al de STRING (ver tablas/string_global_auditoria_mapeo.csv)",
                nrow(auditoria)))


#4. Interacciones STRING entre los genes de entrada
# Se descargan todas las interacciones con puntuación >= min(UMBRALES) y después
# se filtra por cada umbral (la red de 700 es un subconjunto de la de 400).

aristas <- list()

for (f in farmacos) {
  inter <- string_db$get_interactions(nodos[[f]]$STRING_id) %>% as_tibble()

  # STRING devuelve cada par en ambos sentidos: se deja un único registro por par
  aristas[[f]] <- inter %>%
    filter(from != to) %>%
    mutate(a = pmin(from, to), b = pmax(from, to)) %>%
    group_by(a, b) %>%
    summarise(combined_score = max(combined_score), .groups = "drop") %>%
    rename(from = a, to = b)

  message(sprintf("[%s] interacciones únicas con puntuación >= %d: %d",
                  f, min(umbrales), nrow(aristas[[f]])))
}

#5. Construcción de las redes globales

redes   <- list()
resumen <- list()

for (f in farmacos) {
  for (u in umbrales) {
    ar <- aristas[[f]] %>% filter(combined_score >= u)
    g  <- graph_from_data_frame(ar, directed = FALSE, vertices = nodos[[f]])

    comp <- components(g)
    V(g)$grado      <- degree(g)
    V(g)$componente <- comp$membership

    nombre <- paste(f, u, sep = "_")
    redes[[nombre]] <- g
    saveRDS(g, paste0("data_processed/string_global_red_", nombre, ".rds"))

    resumen[[nombre]] <- tibble(
      farmaco                    = f,
      umbral_string              = u,
      n_genes_entrada            = nrow(entrada_gen[[f]]),
      n_nodos                    = vcount(g),
      n_aristas                  = ecount(g),
      grado_medio                = round(mean(degree(g)), 2),
      densidad                   = round(edge_density(g), 4),
      n_nodos_aislados           = sum(degree(g) == 0),
      n_componentes_no_triviales = sum(comp$csize > 1),
      tam_componente_gigante     = max(comp$csize),
      coef_clustering_global     = round(transitivity(g, type = "global"), 3)
    )
  }
}

resumen <- bind_rows(resumen)
print(resumen, width = Inf)

#6. Descripción básica: nodos de mayor grado
# (Solo descriptivo; la identificación de módulos se hace en una etapa posterior

for (nombre in names(redes)) {
  g <- redes[[nombre]]
  cat("\n---", nombre, ": 10 nodos de mayor grado ---\n")
  print(
    tibble(gen = V(g)$gen, grado = V(g)$grado, n_fosfopeptidos = V(g)$n_fosfopeptidos) %>%
      arrange(desc(grado)) %>%
      head(10)
  )
}

#7. Exportación

write_csv(resumen, "tablas/14_string_global_resumen_redes.csv")

for (nombre in names(redes)) {
  g <- redes[[nombre]]

  nod <- igraph::as_data_frame(g, what = "vertices") %>%
    as_tibble() %>%
    rename(STRING_id = name)

  ari <- igraph::as_data_frame(g, what = "edges") %>%
    as_tibble() %>%
    left_join(select(nod, STRING_id, gen_origen = gen), by = c("from" = "STRING_id")) %>%
    left_join(select(nod, STRING_id, gen_destino = gen), by = c("to" = "STRING_id")) %>%
    transmute(source = gen_origen, target = gen_destino, combined_score)

  write_csv(nod, paste0("tablas/14string_global_nodos_", nombre, ".csv"))
  write_csv(ari, paste0("tablas/14string_global_aristas_", nombre, ".csv"))
}

writeLines(capture.output(sessionInfo()), "tablas/14string_global_sessioninfo.txt")

