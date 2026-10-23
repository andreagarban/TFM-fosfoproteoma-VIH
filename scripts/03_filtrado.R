# ==============================================================================
# 03_filtrado.R
#
#TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con inhibidores
#      de tirosina quinasa como estrategia antiviral frente al VIH
#
# Criterios de filtrado: p-valor (t-test) < 0.05 y |AvA| > 1.3
# ==============================================================================

library(tidyverse)

zpq_clean <- readRDS("data_processed/zpq_clean.rds")


#1. SEPARACIÓN PROTEOMA Y FOSFOPROTEOMA

proteoma      <- zpq_clean %>% filter(!is_phospho)
fosfoproteoma <- zpq_clean %>% filter(is_phospho)

# Péptidos totales
nrow(zpq_clean)
# Péptidos no fosforilados
nrow(proteoma)
# Fosfopéptidos
nrow(fosfoproteoma)


## FDR (Benjamini-Hochberg) 
#El FDR se calcula sobre los fosfopéptidos con p-valor calculable de cada comparación 
# y se reporta como capa adicional de robustez estadística, no como filtro.

fosfoproteoma <- fosfoproteoma %>%
  mutate(
    Imatinib_FDR  = p.adjust(Zpq_Ttest_imatinib,  method = "BH"),
    Dasatinib_FDR = p.adjust(Zpq_Ttest_dasatinib, method = "BH")
  )

#2. FILTRADO

# Comparación imatinib vs control
fosfo_imatinib <- fosfoproteoma %>%
  filter(Zpq_Ttest_imatinib < 0.05,
         abs(Zpq_AvA_imatinib) > 1.3)

# Comparación desatinib vs control
fosfo_dasatinib <- fosfoproteoma %>%
  filter(Zpq_Ttest_dasatinib < 0.05,
         abs(Zpq_AvA_dasatinib) > 1.3)

#Resultados del filtrado
nrow(fosfo_imatinib) # Pasan filtro imatinib
n_distinct(fosfo_imatinib$protein) # Proteínas únicas imatinib

nrow(fosfo_dasatinib) # Pasan filtro dasatinib 
n_distinct(fosfo_dasatinib$protein) # Proteínaas únicas desatinib 


#2b. ROBUSTEZ FDR
# De los fosfopéptidos que YA cumplen p<0.05 + |AvA|>1.3, ¿cuántos
# mantienen además FDR<0.10 y FDR<0.05?

robustez_fdr <- tibble(
  Comparación = c("Imatinib", "Dasatinib"),
  `Pasan filtro p<0.05 + |AvA|>1.3` = c(nrow(fosfo_imatinib), nrow(fosfo_dasatinib)),
  `de ellos, FDR < 0.10` = c(
    sum(fosfo_imatinib$Imatinib_FDR   < 0.10, na.rm = TRUE),
    sum(fosfo_dasatinib$Dasatinib_FDR < 0.10, na.rm = TRUE)
  ),
  `de ellos, FDR < 0.05` = c(
    sum(fosfo_imatinib$Imatinib_FDR   < 0.05, na.rm = TRUE),
    sum(fosfo_dasatinib$Dasatinib_FDR < 0.05, na.rm = TRUE)
  )
)
robustez_fdr


# Panel dirigido VIH: FDR por proteína (se mantienen las 4 - búsqueda dirigida por relevancia biológica)
resumen_panel <- fosfoproteoma %>%
  filter(gen %in% c("RANBP2", "HMGA1", "PPIA", "XRCC5")) %>%
  group_by(gen) %>%
  summarise(
    n_fosfopeptidos      = n(),
    Imatinib_sig_actual  = sum(Zpq_Ttest_imatinib  < 0.05 & abs(Zpq_AvA_imatinib)  > 1.3, na.rm = TRUE),
    Imatinib_FDR_min     = min(Imatinib_FDR, na.rm = TRUE),
    Imatinib_sig_FDR10   = sum(Imatinib_FDR  < 0.10, na.rm = TRUE),
    Dasatinib_sig_actual = sum(Zpq_Ttest_dasatinib < 0.05 & abs(Zpq_AvA_dasatinib) > 1.3, na.rm = TRUE),
    Dasatinib_FDR_min    = min(Dasatinib_FDR, na.rm = TRUE),
    Dasatinib_sig_FDR10  = sum(Dasatinib_FDR < 0.10, na.rm = TRUE)
  )

as.data.frame(resumen_panel)


#3. EXPLORACIÓN DE UMBRALES ALTERNATIVOS PARA DESATINIB
# El número de proteínas obtenido con dasatinib es sustancialmente mayor que
# con imatinib. Se comprueba si dicha diferencia responde al umbral empleado

umbrales <- c(1.3, 1.4, 1.5, 1.6, 1.7, 2.0)

exploracion_umbrales <- tibble(umbral = numeric(), n_pep = integer(), n_prot = integer()) 
for (u in umbrales) {
  n_pep <- fosfoproteoma %>%
    filter(Zpq_Ttest_dasatinib < 0.05, abs(Zpq_AvA_dasatinib) > u) %>%
    nrow()
  n_prot <- fosfoproteoma %>%
    filter(Zpq_Ttest_dasatinib < 0.05, abs(Zpq_AvA_dasatinib) > u) %>%
    pull(protein) %>%
    n_distinct()
  cat("Umbral", u, "→", n_pep, "péptidos,", n_prot, "proteínas\n")
  exploracion_umbrales <- bind_rows(
    exploracion_umbrales,
    tibble(umbral = u, n_pep = n_pep, n_prot = n_prot)
  )
}

# Se conserva el umbral de 1.3 en ambas comparaciones y la diferencia se interpreta como un efecto biológico.


#4. DIRECCIÓN DEL CAMBIO

fosfo_imatinib <- fosfo_imatinib %>%
  mutate(direccion = ifelse(Zpq_AvA_imatinib > 0, "Aumenta", "Disminuye"))

fosfo_dasatinib <- fosfo_dasatinib %>%
  mutate(direccion = ifelse(Zpq_AvA_dasatinib > 0, "Aumenta", "Disminuye"))

# Ver cuántos péptidos suben y cuántos bajan en cada comparación
# === IMATINIB === #
table(fosfo_imatinib$direccion)

# === DASATINIB === #
table(fosfo_dasatinib$direccion)

# Dirección del cambio expresada en porcentaje

# Nivel fosfopéptido
direccion_peptidos <- bind_rows(
  fosfo_imatinib  %>% count(direccion) %>% mutate(tratamiento = "Imatinib"),
  fosfo_dasatinib %>% count(direccion) %>% mutate(tratamiento = "Dasatinib")
) %>%
  group_by(tratamiento) %>%
  mutate(total      = sum(n),
         porcentaje = round(100 * n / total, 1)) %>%
  ungroup() %>%
  select(tratamiento, direccion, n, total, porcentaje)

print(direccion_peptidos)

# Nivel proteína
# Una misma proteína puede presentar fosfopéptidos en ambas direcciones, así que la clasificamos en tres grupos.
direccion_proteinas <- bind_rows(
  fosfo_imatinib  %>% mutate(tratamiento = "Imatinib"),
  fosfo_dasatinib %>% mutate(tratamiento = "Dasatinib")
) %>%
  group_by(tratamiento, protein) %>%
  summarise(sube = any(direccion == "Aumenta"),
            baja = any(direccion == "Disminuye"),
            .groups = "drop") %>%
  mutate(patron = case_when(
    sube & baja ~ "Mixto",
    sube        ~ "Solo aumenta",
    baja        ~ "Solo disminuye"
  )) %>%
  count(tratamiento, patron) %>%
  group_by(tratamiento) %>%
  mutate(total      = sum(n),
         porcentaje = round(100 * n / total, 1)) %>%
  ungroup()

print(direccion_proteinas)


# 6. GUARDAR

saveRDS(fosfoproteoma,        "data_processed/fosfoproteoma.rds")
saveRDS(fosfo_imatinib,       "data_processed/fosfo_imatinib.rds")
saveRDS(fosfo_dasatinib,      "data_processed/fosfo_dasatinib.rds")
saveRDS(exploracion_umbrales, "data_processed/exploracion_umbrales.rds")
saveRDS(direccion_peptidos,   "data_processed/direccion_peptidos.rds")
saveRDS(direccion_proteinas,  "data_processed/direccion_proteinas.rds")
saveRDS(robustez_fdr,  "data_processed/robustez_fdr.rds")
saveRDS(resumen_panel, "data_processed/resumen_panel_fdr.rds")
