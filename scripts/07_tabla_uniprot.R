# ==============================================================================
# 07_tabla_uniprot.R
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con ITKs
#      como estrategia antiviral frente al VIH
# ==============================================================================

library(tidyverse)


#1. Cargar y combinar ambas condiciones
cruce_imatinib  <- readRDS("data_processed/cruce_imatinib.rds")
cruce_dasatinib <- readRDS("data_processed/cruce_dasatinib.rds")

tabla_completa <- bind_rows(cruce_imatinib, cruce_dasatinib)

# Comprobación: debería dar 318 + 1002 = 1320
stopifnot(nrow(tabla_completa) == nrow(cruce_imatinib) + nrow(cruce_dasatinib))

# Anotación del bloque de señalización de linfocitos T (script 04, sección 9.3)
# Se reduce a una fila por gen para no duplicar MAPK14 ni PIK3CD al unir
bloque <- readRDS("data_processed/bloque_senalizacion.rds") %>%
  distinct(gen, categoria_hallazgo)

n_antes <- nrow(tabla_completa)

tabla_completa <- tabla_completa %>%
  left_join(bloque, by = "gen")

# Comprobaciones: el join no debe alterar el número de filas
stopifnot(nrow(tabla_completa) == n_antes)

# El bloque debe quedar anotado en 8 filas (6 en dasatinib + 2 en imatinib)
sum(!is.na(tabla_completa$categoria_hallazgo))
tabla_completa %>%
  filter(!is.na(categoria_hallazgo)) %>%
  count(tratamiento, gen)

#Filas totales
nrow(tabla_completa)

nrow(cruce_imatinib) #imatinib
nrow(cruce_dasatinib) #dasatinib

# Orden: tratamiento -> nivel -> gen (alfabético)
tabla_completa <- tabla_completa %>%
  arrange(tratamiento, nivel, gen)


#2. Guardar tabla completa (todas las columnas) en .rds y .csv

saveRDS(tabla_completa, "data_processed/tabla_uniprot_niveles.rds")

write_csv(tabla_completa, "tablas/07_tabla_uniprot_niveles.csv")

sessionInfo()

