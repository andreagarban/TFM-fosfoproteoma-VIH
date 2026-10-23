# ==============================================================================
# Script: 13_comparacion_proteinas_virales.R
# Descripcion: Comprueba si la distribucion de proteinas virales (gp120, Nef,
#              Tat, Vpr, Gag...) observada en los genes de cada bloque con
#              enriquecimiento real (FDR, ver 12_cruce_go_ncbi.R) es
#              caracteristica del bloque o solo refleja el patron general de la
#              base de datos NCBI HIV-1. Para cada proteina viral se calcula el
#              numero y porcentaje de genes con interaccion documentada en
#              (a) el bloque y (b) el conjunto completo de genes significativos
#              de dasatinib de la misma direccion, y se comparan (descriptivo,
#              sin contraste de hipotesis).
#              Se calculan dos versiones del porcentaje segun el denominador:
#                - sobre genes con relacion directa con VIH-1 (nivel 1)
#                - sobre todos los genes del grupo
# ==============================================================================

library(tidyverse)

select <- dplyr::select
filter <- dplyr::filter

# 1. Bloques con enriquecimiento real (FDR) y conjuntos de referencia

genes_bloques <- readRDS("data_processed/genes_bloques_fdr.rds")

# Se excluye senalizacion celular: 6 genes (3 de nivel 1), porcentajes inestables
genes_bloques <- genes_bloques %>%
  filter(tratamiento == "Dasatinib", bloque != "Dasatinib | up | Senalizacion celular")

# Referencia: genes significativos de dasatinib de la misma direccion
# (se reconstruyen igual que en 10_go_kegg.R)
fosfo_dasatinib <- readRDS("data_processed/fosfo_dasatinib.rds") %>%
  filter(!is.na(gen))

genes_dir <- list(
  up   = fosfo_dasatinib %>% filter(Zpq_AvA_dasatinib >  1.3) %>% distinct(gen) %>% pull(gen),
  down = fosfo_dasatinib %>% filter(Zpq_AvA_dasatinib < -1.3) %>% distinct(gen) %>% pull(gen)
)
map_int(genes_dir, length)   # esperado: up 139, down 915


# 2. Tabla a nivel de gen para todo dasatinib 

tabla_completa <- readRDS("data_processed/tabla_uniprot_niveles.rds")

# Une las proteinas virales de varias entradas en un unico texto sin repetidos
unir_virales <- function(x) {
  v <- trimws(unlist(str_split(x, ";")))
  v <- v[!is.na(v) & v != ""]
  paste(sort(unique(v)), collapse = "; ")
}

dasatinib_gen <- tabla_completa %>%
  filter(tratamiento == "Dasatinib", !is.na(gen)) %>%
  mutate(nivel_num = as.integer(str_extract(nivel, "^\\d"))) %>%
  group_by(gen) %>%
  summarise(
    nivel_num         = min(nivel_num),
    proteinas_virales = unir_virales(proteinas_virales),
    .groups = "drop"
  )


# 3. Comprobaciones de coherencia con los numeros ya reportados
cat("\nGenes de dasatinib (nivel gen, sin NA):", nrow(dasatinib_gen),
    "(esperado: 991)\n")
cat("Referencia up:", length(genes_dir$up), "(esperado: 139); down:",
    length(genes_dir$down), "(esperado: 915)\n")
stopifnot(all(genes_bloques$gen %in% dasatinib_gen$gen))

# Las proteinas virales se cuentan solo en genes de nivel 1. Se comprueba que
# los genes de nivel 2/3 no traen proteinas virales rellenas (si las tuvieran,
# habria que decidir si cuentan como interaccion documentada).
n_viral_fuera_nivel1 <- dasatinib_gen %>%
  filter(nivel_num > 1, proteinas_virales != "") %>%
  nrow()
cat("Genes de nivel 2/3 con proteinas_virales no vacio:", n_viral_fuera_nivel1,
    "(esperado: 0)\n")

# 4. Funcion de conteo por proteina viral
# Para un grupo de genes devuelve, por proteina viral: n de genes de nivel 1
# que interactuan con ella y su porcentaje con los dos denominadores.
conteo_virales <- function(df_gen, grupo) {
  n_total  <- nrow(df_gen)
  n_nivel1 <- sum(df_gen$nivel_num == 1)

  df_gen %>%
    filter(nivel_num == 1) %>%
    select(gen, proteinas_virales) %>%
    separate_rows(proteinas_virales, sep = ";\\s*") %>%
    filter(!is.na(proteinas_virales), proteinas_virales != "") %>%
    distinct() %>%
    count(proteinas_virales, name = "n") %>%
    mutate(
      grupo            = grupo,
      pct_sobre_nivel1 = round(100 * n / n_nivel1, 1),
      pct_sobre_total  = round(100 * n / n_total, 1)
    )
}

bloques <- unique(genes_bloques$bloque)

# Devuelve los genes del bloque y los de su referencia (misma direccion)
preparar <- function(b) {
  dir <- str_split(b, " \\| ")[[1]][2]          # "up" o "down"
  genes_b <- genes_bloques %>% filter(bloque == b) %>% pull(gen)
  stopifnot(all(genes_b %in% genes_dir[[dir]]))
  list(bloque     = dasatinib_gen %>% filter(gen %in% genes_b),
       referencia = dasatinib_gen %>% filter(gen %in% genes_dir[[dir]]))
}

tamanos_grupos <- map_dfr(bloques, function(b) {
  p <- preparar(b)
  tibble(bloque          = b,
         n_total_bloque  = nrow(p$bloque),
         n_nivel1_bloque = sum(p$bloque$nivel_num == 1),
         n_total_ref     = nrow(p$referencia),
         n_nivel1_ref    = sum(p$referencia$nivel_num == 1))
})
cat("\nTamaño de los grupos comparados:\n")
print(tamanos_grupos, width = Inf)

# 5. Comparacion bloque vs conjunto de referencia (misma direccion)
comparacion <- map_dfr(bloques, function(b) {
  p <- preparar(b)
  bind_rows(
    conteo_virales(p$bloque,     "bloque"),
    conteo_virales(p$referencia, "referencia")
  ) %>% mutate(bloque = b)
}) %>%
  pivot_wider(
    names_from  = grupo,
    values_from = c(n, pct_sobre_nivel1, pct_sobre_total),
    values_fill = 0
  ) %>%
  mutate(
    # diferencia en puntos porcentuales (bloque - referencia)
    dif_pp_nivel1 = round(pct_sobre_nivel1_bloque - pct_sobre_nivel1_referencia, 1),
    dif_pp_total  = round(pct_sobre_total_bloque  - pct_sobre_total_referencia, 1)
  ) %>%
  arrange(bloque, desc(n_bloque), desc(n_referencia))

# Version A: porcentaje sobre genes con relacion directa con VIH-1 (nivel 1)
tabla_sobre_nivel1 <- comparacion %>%
  select(bloque, proteinas_virales,
         n_bloque, pct_sobre_nivel1_bloque,
         n_referencia, pct_sobre_nivel1_referencia,
         dif_pp_nivel1)

# Version B: porcentaje sobre todos los genes del grupo
tabla_sobre_total <- comparacion %>%
  select(bloque, proteinas_virales,
         n_bloque, pct_sobre_total_bloque,
         n_referencia, pct_sobre_total_referencia,
         dif_pp_total)

cat("\n--- A) Porcentaje sobre genes con relacion directa (nivel 1) ---\n")
print(tabla_sobre_nivel1, n = Inf, width = Inf)

cat("\n--- B) Porcentaje sobre todos los genes del grupo ---\n")
print(tabla_sobre_total, n = Inf, width = Inf)


# 6. Exportar
write_csv(comparacion,        "tablas/13_comparacion_proteinas_virales_dasatinib.csv")
write_csv(tabla_sobre_nivel1, "tablas/13_comparacion_virales_sobre_nivel1.csv")
write_csv(tabla_sobre_total,  "tablas/13_comparacion_virales_sobre_total.csv")
saveRDS(comparacion,          "data_processed/comparacion_proteinas_virales.rds")
