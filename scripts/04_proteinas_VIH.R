# ==============================================================================
# 04_proteinas_VIH.R
#
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con inhibidores
#     de tirosina quinasa como estrategia antiviral frente al VIH
# 
# Selección de las proteínas de interés para la infección por VIH-1 mediante dos estrategias complementarias:
#             (a) búsqueda por palabras clave sobre las categorías funcionales
#             (b) búsqueda dirigida de proteínas concretas descritas en la
#                 literatura, con independencia del resultado del filtrado
# ==============================================================================

library(tidyverse)


fosfoproteoma   <- readRDS("data_processed/fosfoproteoma.rds")
fosfo_imatinib  <- readRDS("data_processed/fosfo_imatinib.rds")
fosfo_dasatinib <- readRDS("data_processed/fosfo_dasatinib.rds")


# PRIMERA PARTE. BÚSQUEDA POR PALABRAS CLAVE

#1. DEFINICIÓN DE LAS PALABRAS CLAVE POR CATEGORÍA FUNCIONAL

keywords_VIH <- c(
  # Poro nuclear
  "nucleoporin", "nuclear pore", "NPC",
  # Complejo de preintegración / transporte nuclear
  "importin", "karyopherin", "TNPO", "transportin",
  # Citoesqueleto
  "actin", "tubulin", "myosin", "vimentin", "filamin",
  "cofilin", "profilin", "talin", "vinculin",
  # Espliceosoma y factores de splicing
  "splicing", "spliceosom", "snRNP", "heterogeneous nuclear",
  "hnRNP", "serine/arginine", "SR protein"
)

# Construir patrón de búsqueda (insensible a mayúsculas)
patron <- paste(keywords_VIH, collapse = "|")

# Aplicación del filtro sobre los conjuntos ya seleccionados

fosfo_imatinib_VIH <- fosfo_imatinib %>%
  filter(str_detect(description, regex(patron, ignore_case = TRUE)))

fosfo_dasatinib_VIH <- fosfo_dasatinib %>%
  filter(str_detect(description, regex(patron, ignore_case = TRUE)))

# Resultados
# === IMATINIB — proteínas de interés VIH === #
nrow(fosfo_imatinib_VIH) # fosfopéptidos 
n_distinct(fosfo_imatinib_VIH$protein) # proteínas únicas

# === DESATINIB — proteínas de interés VIH === #
nrow(fosfo_dasatinib_VIH) # fosfopéptidos
n_distinct(fosfo_dasatinib_VIH$protein) # proteínas únicas

# Listado de proteínas recuperadas, para su revisión manual.
# La búsqueda por palabras clave puede recuperar proteínas que no pertenecen a ninguna de las categorías de interés
# por lo que el listado requiere revisión y validación posterior.
proteinas_keywords_imatinib <- fosfo_imatinib_VIH %>%
  distinct(uniprot = fasta, gen, protein, description)
View(proteinas_keywords_imatinib)   

proteinas_keywords_dasatinib <- fosfo_dasatinib_VIH %>%
  distinct(uniprot = fasta, gen, protein, description)
View(proteinas_keywords_dasatinib)


# SEGUNDA PARTE. BÚSQUEDA DIRIGIDA DE PROTEÍNAS RELEVANTES EN VIH
# La búsqueda se realiza sobre el fosfoproteoma completo, sin aplicar los
# criterios de filtrado, con el fin de distinguir entre proteínas no
# identificadas en el experimento y proteínas identificadas cuya fosforilación
# no se modifica de forma apreciable.


#1. DEFINICIÓN DEL PANEL DE PROTEINAS

genes_VIH <- tribble(
  ~gen,       ~categoria,
  # Poro nuclear y transporte nucleocitoplasmático
  "RANBP2",   "Poro nuclear / transporte nuclear",
  "NUP62",    "Poro nuclear / transporte nuclear",
  "NUP98",    "Poro nuclear / transporte nuclear",
  "NUP107",   "Poro nuclear / transporte nuclear",
  "NUP160",   "Poro nuclear / transporte nuclear",
  "NUP85",    "Poro nuclear / transporte nuclear",
  "TNPO3",    "Poro nuclear / transporte nuclear",
  "KPNB1",    "Poro nuclear / transporte nuclear",
  "IPO7",     "Poro nuclear / transporte nuclear",
  "RAN",      "Poro nuclear / transporte nuclear",
  # Integración y tráfico intracelular
  "CPSF6",    "Integración / tráfico intracelular",
  "BANF1",    "Integración / tráfico intracelular",
  "HMGA1",    "Integración / tráfico intracelular",
  "SMARCB1",  "Integración / tráfico intracelular",
  # Factores celulares clásicos del ciclo del VIH-1
  "PPIA",     "Factor celular clásico VIH",
  "XRCC5",    "Factor celular clásico VIH"
)


#2. BÚSQUEDA SOBRE EL FOSFOPROTEOMA COMPLETO

busqueda_VIH <- fosfoproteoma %>%
  filter(gen %in% genes_VIH$gen) %>%
  left_join(genes_VIH, by = "gen") %>%
  mutate(
    signif_imatinib  = Zpq_Ttest_imatinib  < 0.05 & abs(Zpq_AvA_imatinib)  > 1.3,
    signif_dasatinib = Zpq_Ttest_dasatinib < 0.05 & abs(Zpq_AvA_dasatinib) > 1.3
  ) %>%
  select(gen, categoria, uniprot = fasta, protein, sequence,
         Zpq_AvA_imatinib,  Zpq_Ttest_imatinib,  Imatinib_FDR,  signif_imatinib,
         Zpq_AvA_dasatinib, Zpq_Ttest_dasatinib, Dasatinib_FDR, signif_dasatinib)



# Clasificación en los tres escenarios posibles

detectados    <- unique(busqueda_VIH$gen)
no_detectados <- setdiff(genes_VIH$gen, detectados)

# Proteínas identificadas con al menos un fosfopéptido significativo
con_cambio <- busqueda_VIH %>%
  filter(signif_imatinib | signif_dasatinib) %>%
  pull(gen) %>%
  unique()

# Proteínas identificadas sin ningún fosfopéptido significativo
sin_cambio <- setdiff(detectados, con_cambio)

resumen_busqueda <- genes_VIH %>%
  mutate(escenario = case_when(
    gen %in% con_cambio ~ "Identificada con cambio significativo",
    gen %in% sin_cambio ~ "Identificada sin cambio significativo",
    TRUE                ~ "No identificada"
  ))

#Resultados: Búsqueda dirigida
print(resumen_busqueda, n = Inf)

length(con_cambio) # Identificadas con cambio significativo
length(sin_cambio) # Identificadas sin cambio significativo
length(no_detectados) # No identificadas

# Detalle de los fosfopéptidos significativos
peptidos_signif_VIH <- busqueda_VIH %>%
  filter(signif_imatinib | signif_dasatinib)

View(peptidos_signif_VIH)


# Comprobación de entradas redundantes
# Una misma proteína puede aparecer asociada a varios identificadores
# (entrada revisada y entradas no revisadas)

busqueda_VIH %>%
  distinct(gen, uniprot) %>%
  count(gen, name = "n_identificadores") %>%
  filter(n_identificadores > 1) %>%
  print()

busqueda_VIH %>%
  distinct(gen, uniprot, protein)

# El gen RAN aparece asociado a tres códigos UniProt: una entrada revisada
# (P62826, RAN_HUMAN, SwissProt) y dos entradas no revisadas (H0YFC6 y J3KQE5, TrEMBL), 
# derivadas del mismo gen. Se descartan para evitar duplicar artificialmente el recuento de proteínas

busqueda_VIH <- busqueda_VIH %>%
  filter(!uniprot %in% c("H0YFC6", "J3KQE5"))

# Verificación: 13 proteínas únicas del panel dirigido
n_distinct(busqueda_VIH$gen)
n_distinct(busqueda_VIH$uniprot)


# GUARDAR

saveRDS(fosfo_imatinib_VIH,           "data_processed/fosfo_imatinib_VIH.rds")
saveRDS(fosfo_dasatinib_VIH,          "data_processed/fosfo_dasatinib_VIH.rds")
saveRDS(proteinas_keywords_imatinib,  "data_processed/proteinas_keywords_imatinib.rds")
saveRDS(proteinas_keywords_dasatinib, "data_processed/proteinas_keywords_dasatinib.rds")
saveRDS(busqueda_VIH,                 "data_processed/busqueda_VIH.rds")
saveRDS(resumen_busqueda,             "data_processed/resumen_busqueda.rds")
