# ==============================================================================
# Script: 11_terminos_interes_vih.R
#
# # TFM: "Analisis proteomico y fosfoproteomico de PBMCs tratadas con
# inhibidores de tirosina quinasa como estrategia antiviral frente al VIH"
#
# Filtra, entre los terminos GO/KEGG/Reactome enriquecidos de imatinib y dasatinib, 
# los relacionados con las categorias de interes biologico del TFM: 
# poro nuclear, complejo de preinfeccion/ preintegracion, citoesqueleto, 
# spliceosoma/splicing y señalización por linfocitos T.
# ==============================================================================

library(tidyverse)

select <- dplyr::select
filter <- dplyr::filter

# 1. Cargar datos
# Terminos enriquecidos de 10_go_kegg.R (corrida sin correccion FDR), con la marca
# de significacion: "FDR (BH < 0,05)" o "Solo nominal (p < 0,05)".
# Incluye GO_BP, GO_MF, GO_CC, KEGG y Reactome de los 4 conjuntos.
terminos_enriquecidos <- read_csv("tablas/10_terminos_nivel_significacion.csv",
                                  show_col_types = FALSE) %>%
  separate(conjunto, into = c("tratamiento", "direccion"), sep = "_")


# 2. Palabras clave por categoria de interes 
#NOTA: los terminos GO/KEGG estan en ingles

palabras_clave <- list(
  "Poro nuclear"               = c("nuclear pore", "nucleoporin"),
  "Complejo de preintegracion" = c("preintegration", "\\bviral\\b"),
  "Citoesqueleto"              = c("\\bactin\\b", "\\bmicrotubule",
                                   "cytoskeleton", "cytoskeletal"),
  "Spliceosoma/splicing"       = c("splic", "spliceosome"),
  "Senalizacion celular"  = c("\\bt cell", "\\bt-cell", "nf-kappa ?b",
                                   "jak-stat", "\\bstat\\b", "\\bmapk",
                                   "\\bpi3k", "phosphatidylinositol 3-kinase",
                                   "interleukin")
)


# 3. Funcion de filtrado
buscar_terminos_interes <- function(tabla, palabras_clave) {
  map_dfr(names(palabras_clave), function(cat_interes) {
    patron <- paste(palabras_clave[[cat_interes]], collapse = "|")
    tabla %>%
      filter(str_detect(tolower(Description), patron)) %>%
      transmute(tratamiento, direccion,
                categoria,                        # GO_BP / GO_MF / GO_CC / KEGG / Reactome
                categoria_interes = cat_interes,  # categoria biologica del TFM
                significacion,                    # FDR o solo nominal
                Description,
                p_nominal, p_ajustado,            # p_ajustado: solo en los terminos con FDR
                Count, geneID)
  })
}

# 4. Ejecutar y exportar 
terminos_interes <- buscar_terminos_interes(terminos_enriquecidos, palabras_clave)

print(terminos_interes)
terminos_interes %>%
  count(tratamiento, direccion, categoria_interes, significacion) %>%
  print(n = Inf)

saveRDS(terminos_interes, "data_processed/terminos_interes.rds")
write_csv(terminos_interes, "tablas/11_terminos_interes_vih.csv")


