# ==============================================================================
# 06_volcano.R
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con ITKs
#      como estrategia antiviral frente al VIH
#
# Objetivo: Volcano plots (imatinib-vs-control, dasatinib-vs-control) sobre
#           el conjunto completo de fosfopéptidos, marcando el panel dirigido de proteínas VIH-
#           relevantes con cambios significativos: RANBP2, HMGA1, PPIA, XRCC5.
# ==============================================================================

library(tidyverse)
library(ggrepel)

zpq_clean <- readRDS("data_processed/zpq_clean.rds")

#1. Panel dirigido: proteínas VIH-relevantes con cambios significativos

panel_dirigido <- tibble(
  uniprot_entry = c("RBP2_HUMAN", "HMGA1_HUMAN", "PPIA_HUMAN", "XRCC5_HUMAN"),
  gen_panel     = c("RANBP2", "HMGA1", "PPIA", "XRCC5")
)


#2. Filtrar fosfopéptidos

phospho <- zpq_clean %>% filter(is_phospho == TRUE)


# 3. Función: preparar datos de volcano para una condición y clasificar
#    significancia según los filtros estándar (p < 0.05, |AvA| > 1.3)

prep_volcano <- function(df, condicion) {

  col_ava   <- paste0("Zpq_AvA_", condicion)
  col_ttest <- paste0("Zpq_Ttest_", condicion)

  df %>%
    select(RowProtein, protein, gen, all_of(c(col_ava, col_ttest))) %>%
    rename(AvA = all_of(col_ava), pval = all_of(col_ttest)) %>%
    filter(!is.na(AvA), !is.na(pval)) %>%
    mutate(
      neg_log10_p = -log10(pval),
      categoria = case_when(
        pval < 0.05 & AvA >  1.3 ~ "Aumenta",
        pval < 0.05 & AvA < -1.3 ~ "Disminuye",
        TRUE                     ~ "No significativo"
      ),
      panel = protein %in% panel_dirigido$uniprot_entry
    )
}

volcano_imatinib  <- prep_volcano(phospho, "imatinib")
volcano_dasatinib <- prep_volcano(phospho, "dasatinib")


# 4. Función de plotting

plot_volcano <- function(datos, titulo) {
  
  etiquetas <- datos %>%
    filter(panel) %>%
    left_join(panel_dirigido, by = c("protein" = "uniprot_entry"))
  
  etiquetas_texto <- etiquetas %>%
    group_by(gen_panel) %>%
    slice_max(abs(AvA), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ggplot(datos, aes(x = AvA, y = neg_log10_p)) +
    geom_point(aes(color = categoria), alpha = 0.5, size = 1.2) +
    geom_point(data = etiquetas, color = "black", size = 2.8,
               shape = 21, fill = "yellow") +
    geom_text_repel(data = etiquetas_texto, aes(label = gen_panel),
                    fontface = "bold", size = 4, box.padding = 0.6) +
    geom_vline(xintercept = c(-1.3, 1.3), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Aumenta"           = "#D55E00",
                                   "Disminuye"         = "#0072B2",
                                   "No significativo"  = "grey70")) +
    labs(
      title = titulo,
      x     = "AvA",
      y     = expression(-log[10](italic(p)*"-value")),
      color = "Significancia"
    ) +
    theme_minimal(base_size = 13)
}

p_imatinib  <- plot_volcano(volcano_imatinib,  "Volcano plot - Imatinib vs Control (fosfoproteoma)")
p_dasatinib <- plot_volcano(volcano_dasatinib, "Volcano plot - Dasatinib vs Control (fosfoproteoma)")

p_imatinib
p_dasatinib


#5. Guardar figuras

ggsave("figuras/volcano_imatinib.png",  p_imatinib,  width = 8, height = 6, dpi = 300)
ggsave("figuras/volcano_dasatinib.png", p_dasatinib, width = 8, height = 6, dpi = 300)
