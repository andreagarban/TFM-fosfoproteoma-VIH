# ==============================================================================
# Script: 12_cruce_go_ncbi.R
#
# # TFM: "Analisis proteomico y fosfoproteomico de PBMCs tratadas con
# inhibidores de tirosina quinasa como estrategia antiviral frente al VIH"
#
# Extrae los genes de los terminos GO/KEGG/Reactome que superan la correccion
# de BH (FDR < 0,05) en cada bloque (farmaco x direccion x categoria de interes)
# con enriquecimiento real, y los cruza con tabla_completa (clasificacion de
# evidencia de interaccion con VIH-1, NCBI HIV-1 Human Interaction Database).
# Se excluye el bloque de preintegracion: su unico termino con FDR
# (Viral mRNA Translation) es un termino viral general, no especifico del PIC.
# ==============================================================================

library(tidyverse)
select <- dplyr::select
filter <- dplyr::filter

#0. Carga de datos 
terminos_interes  <- readRDS("data_processed/terminos_interes.rds")
tabla_completa <- readRDS("data_processed/tabla_uniprot_niveles.rds")


#1. Genes de los terminos con FDR de cada bloque

# Bloques con al menos un termino con FDR (excluido preintegracion)
bloques_fdr <- terminos_interes %>%
  filter(significacion == "FDR (BH < 0,05)",
         categoria_interes != "Complejo de preintegracion") %>%
  distinct(tratamiento, direccion, categoria_interes)

print(bloques_fdr)   # deben salir exactamente 4 bloques, todos de Dasatinib

# Genes de los terminos con FDR de cada bloque
genes_bloques <- terminos_interes %>%
  filter(significacion == "FDR (BH < 0,05)",
         categoria_interes != "Complejo de preintegracion") %>%
  mutate(bloque = paste(tratamiento, direccion, categoria_interes, sep = " | ")) %>%
  select(bloque, tratamiento, gen = geneID) %>%
  separate_rows(gen, sep = "/") %>%
  distinct()

genes_bloques %>% count(bloque, name = "n_genes")   # genes unicos por bloque


#2. Cruce con tabla_completa

cruce_go_ncbi <- genes_bloques %>%
  inner_join(
    tabla_completa %>%
      select(tratamiento, gen, proteina, uniprot, nivel,
             proteinas_virales, tipos_interaccion),
    by = c("tratamiento", "gen"),
    relationship = "many-to-many"
  ) %>%
  distinct()

cruce_go_ncbi %>% count(bloque, name = "n_filas")   # filas encontradas por bloque


#3. Resumen por nivel de evidencia (un nivel por gen: el de menor numero)

resumen_niveles_gen <- cruce_go_ncbi %>%
  mutate(nivel_num = as.integer(str_extract(nivel, "^\\d"))) %>%
  group_by(bloque, gen) %>%
  slice_min(nivel_num, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(bloque, nivel, name = "n_genes") %>%
  group_by(bloque) %>%
  mutate(pct = round(100 * n_genes / sum(n_genes), 1)) %>%
  ungroup() %>%
  arrange(bloque, nivel)

print(resumen_niveles_gen, n = Inf)


#4. Genes de los bloques sin fila en tabla_completa

genes_sin_match <- anti_join(genes_bloques, cruce_go_ncbi, by = c("bloque", "gen"))

nrow(genes_sin_match)
print(genes_sin_match, n = Inf)


# 5. Lista final: genes de cada bloque con relacion directa (nivel 1)

lista_nivel1 <- cruce_go_ncbi %>%
  filter(nivel == "1. Relación directa con VIH-1") %>%
  select(bloque, gen, proteinas_virales, tipos_interaccion) %>%
  distinct() %>%
  arrange(bloque, gen)

count(lista_nivel1, bloque, name = "n_genes_nivel1")
print(lista_nivel1, n = Inf)


#6. Frecuencia de proteinas virales entre los genes de nivel 1 de cada bloque

frecuencia_proteinas_virales <- lista_nivel1 %>%
  separate_rows(proteinas_virales, sep = ";\\s*") %>%
  distinct(bloque, gen, proteinas_virales) %>%
  count(bloque, proteinas_virales, name = "n_genes") %>%
  arrange(bloque, desc(n_genes))

print(frecuencia_proteinas_virales, n = Inf)


#7. Exportar
write_csv(cruce_go_ncbi, "tablas/12_cruce_go_ncbi_bloques_fdr.csv")
write_csv(resumen_niveles_gen, "tablas/12_resumen_niveles_bloques_fdr.csv")
write_csv(lista_nivel1, "tablas/12_nivel1_bloques_fdr.csv")
write_csv(frecuencia_proteinas_virales, "tablas/12_frecuencia_proteinas_virales_bloques_fdr.csv")

saveRDS(genes_bloques, "data_processed/genes_bloques_fdr.rds")
saveRDS(cruce_go_ncbi, "data_processed/cruce_go_ncbi_bloques_fdr.rds")


