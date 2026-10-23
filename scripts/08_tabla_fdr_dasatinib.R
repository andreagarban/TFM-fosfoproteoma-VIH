# ==============================================================================
# 08_tabla_fdr_dasatinib.R
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con ITKs
#      como estrategia antiviral frente al VIH
#
# Tabla de fosfopéptidos de dasatinib que, cumpliendo el filtro p<0.05 +
# |AvA|>1.3, mantienen además FDR<0.10
# ==============================================================================

library(tidyverse)

fosfo_dasatinib <- readRDS("data_processed/fosfo_dasatinib.rds")

tabla_fdr_dasatinib <- fosfo_dasatinib %>%
  filter(Dasatinib_FDR < 0.10) %>%
  transmute(
    Gen          = gen,
    UniProt      = fasta,
    Fosfopeptido = sequence,
    AvA          = round(Zpq_AvA_dasatinib, 2),
    p_valor      = signif(Zpq_Ttest_dasatinib, 3),
    FDR          = signif(Dasatinib_FDR, 3)
  ) %>%
  arrange(FDR)

nrow(tabla_fdr_dasatinib)  # debería dar 95

saveRDS(tabla_fdr_dasatinib, "data_processed/tabla_fdr_dasatinib.rds")
write_csv(tabla_fdr_dasatinib, "tablas/tabla_fdr_dasatinib.csv")

