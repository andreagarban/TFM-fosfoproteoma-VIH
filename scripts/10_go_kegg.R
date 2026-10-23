# ==============================================================================
# Script: 10_go_kegg.R
#
# # TFM: "Analisis proteomico y fosfoproteomico de PBMCs tratadas con
# inhibidores de tirosina quinasa como estrategia antiviral frente al VIH"
#
# Enriquecimiento funcional GO (BP, MF, CC), KEGG y Reactome para los
# fosfopeptidos significativos de imatinib y dasatinib.
# ==============================================================================


library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)

# Alias explicitos para evitar conflictos 
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
slice  <- dplyr::slice


# 1. Cargar datos
fosfo_imatinib  <- readRDS("data_processed/fosfo_imatinib.rds")
fosfo_dasatinib <- readRDS("data_processed/fosfo_dasatinib.rds")
fosfoproteoma   <- readRDS("data_processed/fosfoproteoma.rds") 

# 2. Universo de fondo 
# Se usa como background el conjunto de genes identificados en el experimento
# (fosfoproteoma completo, sin filtrar por significancia)

universo_genes <- fosfoproteoma %>%
  filter(!is.na(gen)) %>%
  filter(!is.na(Zpq_Ttest_imatinib) | !is.na(Zpq_Ttest_dasatinib)) %>%
  distinct(gen) %>%
  pull(gen)

universo_entrez <- bitr(universo_genes, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db, drop = TRUE) %>%
  distinct(ENTREZID, .keep_all = TRUE)

length(universo_genes) #Universo de fondo (genes unicos identificados)
nrow(universo_entrez) # mapeados a Entrez ID 
round(100 * nrow(universo_entrez) / length(universo_genes), 1) #porcentaje de mapeo
setdiff(universo_genes, bitr(universo_genes, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$SYMBOL) #qué símbolos se pierden

# 3. Conjuntos de genes por farmaco y direccion del cambio
# fosfo_imatinib / fosfo_dasatinib ya cumplen p < 0.05 y |AvA| > 1.3;
# aqui solo se separan por el signo del AvA correspondiente a cada farmaco.

UMBRAL_AVA <- 1.3
MIN_GS <- 10    # tamano minimo de conjunto de genes (valor por defecto de clusterProfiler)
MAX_GS <- 500   # tamano maximo (valor por defecto)

genes_direccion <- function(datos, col_AvA, sentido) {
  datos <- datos %>% filter(!is.na(gen))
  if (sentido == "up") {
    datos <- datos %>% filter(.data[[col_AvA]] >  UMBRAL_AVA)
  } else {
    datos <- datos %>% filter(.data[[col_AvA]] < -UMBRAL_AVA)
  }
  datos %>% distinct(gen) %>% pull(gen)
}

conjuntos_genes <- list(
  Imatinib_up    = genes_direccion(fosfo_imatinib,  "Zpq_AvA_imatinib",  "up"),
  Imatinib_down  = genes_direccion(fosfo_imatinib,  "Zpq_AvA_imatinib",  "down"),
  Dasatinib_up   = genes_direccion(fosfo_dasatinib, "Zpq_AvA_dasatinib", "up"),
  Dasatinib_down = genes_direccion(fosfo_dasatinib, "Zpq_AvA_dasatinib", "down")
)

# Diagnostico: tamano de cada conjunto y genes presentes en ambas direcciones
diag_conjuntos <- tibble(
  conjunto = names(conjuntos_genes),
  n_genes  = map_int(conjuntos_genes, length)
)
solape <- tibble(
  farmaco = c("Imatinib", "Dasatinib"),
  n_genes_en_ambas_direcciones = c(
    length(intersect(conjuntos_genes$Imatinib_up,  conjuntos_genes$Imatinib_down)),
    length(intersect(conjuntos_genes$Dasatinib_up, conjuntos_genes$Dasatinib_down))
  )
)
print(diag_conjuntos)
print(solape)
write_csv(diag_conjuntos, "tablas/10_tamano_conjuntos.csv")
write_csv(solape, "tablas/10_solapamiento_up_down.csv")

#Un gen con fosfopéptidos que suben y otros que bajan aparece en los dos conjuntos.
#Anotar en Métodos cuántos genes se solapan y por qué. 


# 4. Funcion de enriquecimiento GO + KEGG + Reactome
# ajuste: "BH" (con correccion FDR) o "none" (sin correccion)

enriquecer <- function(genes, universo_entrez, ajuste, nombre = "") {
  
  genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db, drop = TRUE) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  stopifnot(all(genes_entrez$ENTREZID %in% universo_entrez$ENTREZID))
  
  # Genes de entrada y mapeados a Entrez (se imprime una sola vez por conjunto)
  if (ajuste == "BH") {
    cat(nombre, "- Genes:", length(genes), "; mapeados a Entrez:",
        nrow(genes_entrez), "(",
        round(100 * nrow(genes_entrez) / length(genes), 1), "%)\n")
  }
  
  resultados <- list()
  
  # GO: Biological Process, Molecular Function, Cellular Component
  for (onto in c("BP", "MF", "CC")) {
    go_res <- enrichGO(
      gene          = genes_entrez$ENTREZID,
      universe      = universo_entrez$ENTREZID,
      OrgDb         = org.Hs.eg.db,
      ont           = onto,
      pAdjustMethod = ajuste,
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 1,
      minGSSize     = MIN_GS,
      maxGSSize     = MAX_GS,
      readable      = TRUE
    )
    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      go_res <- simplify(go_res, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }
    resultados[[paste0("GO_", onto)]] <- go_res
  }
  
  # KEGG
  kegg_res <- enrichKEGG(
    gene          = genes_entrez$ENTREZID,
    universe      = universo_entrez$ENTREZID,
    organism      = "hsa",
    pAdjustMethod = ajuste,
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 1,
    minGSSize     = MIN_GS,
    maxGSSize     = MAX_GS
  )
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  resultados[["KEGG"]] <- kegg_res
  
  # Reactome
  resultados[["Reactome"]] <- enrichPathway(
    gene          = genes_entrez$ENTREZID,
    universe      = universo_entrez$ENTREZID,
    organism      = "human",
    pAdjustMethod = ajuste,
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 1,
    minGSSize     = MIN_GS,
    maxGSSize     = MAX_GS,
    readable      = TRUE
  )
  
  resultados
}


# 5. Ejecutar: 4 conjuntos x 2 modos (con FDR = BH, sin FDR = none)
ajustes <- c(BH = "BH", sin_FDR = "none")

resultados <- map(ajustes, function(aj) {
  imap(conjuntos_genes, function(g, nombre) {
    enriquecer(g, universo_entrez, aj, nombre)
  })
})
# Acceso: resultados[["BH"]][["Dasatinib_down"]][["GO_BP"]]


# Diagnóstico de la discrepancia 1.002 proteínas vs 992 genes en dasatinib
sum(is.na(fosfo_dasatinib$gen)) #NA en gen

fosfo_dasatinib %>% distinct(protein) %>% nrow() # nº proteínas (UniProt)
fosfo_dasatinib %>% distinct(gen) %>% nrow() # nº genes (símbolo)

fosfo_dasatinib %>%
  filter(!is.na(gen)) %>%
  group_by(gen) %>%
  summarise(n_proteinas = n_distinct(protein), .groups = "drop") %>%
  filter(n_proteinas > 1) # genes con >1 proteína asociada


#6. Exportar tablas y objetos
for (aj in names(resultados)) {
  for (conj in names(resultados[[aj]])) {
    for (cat in names(resultados[[aj]][[conj]])) {
      res <- resultados[[aj]][[conj]][[cat]]
      if (is.null(res)) next
      nombre <- paste(tolower(cat), tolower(conj), tolower(aj), sep = "_")
      tabla <- as.data.frame(res)
      if (nrow(tabla) > 0) write_csv(tabla, paste0("tablas/", nombre, ".csv"))
      saveRDS(res, paste0("data_processed/", nombre, ".rds"))
    }
  }
}


#7. Resumen de terminos significativos (con vs sin FDR)
resumen <- map_dfr(names(resultados), function(aj) {
  map_dfr(names(resultados[[aj]]), function(conj) {
    tibble(
      correccion     = aj,
      conjunto       = conj,
      n_genes        = length(conjuntos_genes[[conj]]),
      categoria      = names(resultados[[aj]][[conj]]),
      n_terminos_sig = map_int(resultados[[aj]][[conj]], function(x) {
        if (is.null(x)) return(0L)
        tab <- as.data.frame(x)
        if (nrow(tab) == 0) return(0L)
        sum(tab$p.adjust < 0.05)
      })
    )
  })
})

resumen_ancho <- resumen %>%
  pivot_wider(names_from = correccion, values_from = n_terminos_sig)

print(resumen_ancho, n = Inf)
write_csv(resumen_ancho, "tablas/10_resumen_enriquecimiento_direccion.csv")

resumen %>%
  filter(n_terminos_sig == 0) %>%
  arrange(conjunto, correccion, categoria) %>%
  print(n = Inf)

# 7b. Nivel de significacion de cada termino: FDR vs solo nominal
etiquetar <- function(x, etiqueta, excluir = character(0)) {
  if (is.null(x)) return(NULL)
  d <- as_tibble(as.data.frame(x))
  if (nrow(d) == 0) return(NULL)
  d %>% filter(!ID %in% excluir) %>% mutate(significacion = etiqueta)
}

# Terminos con FDR: de la ejecucion con correccion (BH).
# Solo nominales: de la ejecucion sin correccion, los que no estan en la de BH.
tabla_significacion <- map_dfr(names(conjuntos_genes), function(conj) {
  map_dfr(names(resultados$sin_FDR[[conj]]), function(cat) {
    d_bh   <- etiquetar(resultados$BH[[conj]][[cat]], "FDR (BH < 0,05)")
    ids_bh <- if (is.null(d_bh)) character(0) else d_bh$ID
    d_nom  <- etiquetar(resultados$sin_FDR[[conj]][[cat]],
                        "Solo nominal (p < 0,05)", excluir = ids_bh)
    bind_rows(d_bh, d_nom) %>% mutate(conjunto = conj, categoria = cat)
  })
}) %>%
  mutate(p_nominal  = pvalue,
         p_ajustado = if_else(significacion == "FDR (BH < 0,05)", p.adjust, NA_real_))

write_csv(tabla_significacion, "tablas/10_terminos_nivel_significacion.csv")
tabla_significacion %>%
  count(conjunto, categoria, significacion) %>%
  as_tibble() %>%
  print(n = Inf)


#8. Dotplots
guardar_dotplot <- function(resultado, titulo, archivo) {
  if (!is.null(resultado) && nrow(as.data.frame(resultado)) > 0) {
    p <- dotplot(resultado, showCategory = 15) + ggtitle(titulo)
    ggsave(paste0("figuras/", archivo), p, width = 8, height = 6, dpi = 300)
  }
}

for (aj in names(resultados)) {
  for (conj in names(resultados[[aj]])) {
    for (cat in names(resultados[[aj]][[conj]])) {
      guardar_dotplot(
        resultados[[aj]][[conj]][[cat]],
        paste(cat, "-", gsub("_", " ", conj), "-", aj),
        paste0("dotplot_", tolower(cat), "_", tolower(conj), "_", tolower(aj), ".png")
      )
    }
  }
}


#9. Redundancia entre terminos Reactome (solapamiento de genes, Jaccard)
agrupar_terminos <- function(res, h = 0.5) {
  tab <- as.data.frame(res)
  if (nrow(tab) < 2) return(NULL)
  g <- strsplit(tab$geneID, "/")
  jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  m <- outer(seq_along(g), seq_along(g),
             Vectorize(function(i, j) jac(g[[i]], g[[j]])))
  hc <- hclust(as.dist(1 - m), method = "average")
  tab$grupo <- cutree(hc, h = h)
  tab
}

cortes <- c(0.3, 0.5, 0.7)

redundancia <- map_dfr(names(conjuntos_genes), function(conj) {
  map_dfr(cortes, function(h) {
    tab <- agrupar_terminos(resultados$BH[[conj]]$Reactome, h)
    if (is.null(tab)) return(NULL)
    tibble(conjunto = conj, corte = h,
           n_terminos = nrow(tab),
           n_genes_distintos = n_distinct(unlist(strsplit(tab$geneID, "/"))),
           n_grupos = n_distinct(tab$grupo))
  })
})
print(redundancia, n = Inf)
write_csv(redundancia, "tablas/10_redundancia_reactome.csv")

# Terminos de Reactome de dasatinib (aumento) con su grupo (corte 0.5)
tab_grupos <- agrupar_terminos(resultados$BH$Dasatinib_up$Reactome, 0.5)
write_csv(tab_grupos, "tablas/10_reactome_dasatinib_up_grupos.csv")


# 10. Versiones de paquetes y fecha de consulta (KEGG se consulta en linea)
paquetes <- c("clusterProfiler", "ReactomePA", "org.Hs.eg.db", "enrichplot")
writeLines(
  c(paste("Fecha de ejecucion:", Sys.Date()),
    paste(paquetes, sapply(paquetes, function(p) as.character(packageVersion(p))))),
  "tablas/10_versiones.txt"
)

