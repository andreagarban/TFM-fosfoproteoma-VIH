# ==============================================================================
# 01_lectura_datos.R
#
# TFM: Análisis proteómico y fosfoproteómico de PBMCs tratadas con inhibidores
#      de tirosina quinasa como estrategia antiviral frente al VIH
# ==============================================================================


library(readxl)
library(tidyverse)

#1. LECTURA
zpq_raw <- read_excel("data_raw/Results_Coiras_ProtAndPhosphoCorrect.xlsx", 
                      sheet = "Zpq",
                      skip = 3)

#2. RENOMBRAR LAS COLUMNAS 
colnames(zpq_raw) <- c(
  
  # Columnas de identificación (1-6)
  "RowProtein", "fasta", "protein", "description", "sequence", "nP",
  
  # Zq por muestra - 11 muestras (7-17)
  "Zq_126t", "Zq_127N", "Zq_127C", "Zq_128N", "Zq_128C",
  "Zq_129N", "Zq_129C", "Zq_130N", "Zq_130C", "Zq_131N", "Zq_131C",
  
  # Zq estadísticas: Imatinib vs Control (18-20)
  "Zq_AvA_imatinib", "Zq_logA_imatinib", "Zq_Ttest_imatinib",
  
  # Zq estadísticas: Dasatinib vs Control (21-23)
  "Zq_AvA_dasatinib", "Zq_logA_dasatinib", "Zq_Ttest_dasatinib",
  
  # Zp por muestra - 11 muestras (24-34)
  "Zp_126t", "Zp_127N", "Zp_127C", "Zp_128N", "Zp_128C",
  "Zp_129N", "Zp_129C", "Zp_130N", "Zp_130C", "Zp_131N", "Zp_131C",
  
  # Zp estadísticas: Imatinib vs Control (35-37)
  "Zp_AvA_imatinib", "Zp_logA_imatinib", "Zp_Ttest_imatinib",
  
  # Zp estadísticas: Dasatinib vs Control (38-40)
  "Zp_AvA_dasatinib", "Zp_logA_dasatinib", "Zp_Ttest_dasatinib",
  
  # Zpq por muestra - 11 muestras (41-51)
  "Zpq_126t", "Zpq_127N", "Zpq_127C", "Zpq_128N", "Zpq_128C",
  "Zpq_129N", "Zpq_129C", "Zpq_130N", "Zpq_130C", "Zpq_131N", "Zpq_131C",
  
  # Zpq estadísticas: Imatinib vs Control (52-54)
  "Zpq_AvA_imatinib", "Zpq_logA_imatinib", "Zpq_Ttest_imatinib",
  
  # Zpq estadísticas: Dasatinib vs Control (55-57)
  "Zpq_AvA_dasatinib", "Zpq_logA_dasatinib", "Zpq_Ttest_dasatinib",
  
  # FDRq por muestra (58-68)
  "FDRq_126t", "FDRq_127N", "FDRq_127C", "FDRq_128N", "FDRq_128C",
  "FDRq_129N", "FDRq_129C", "FDRq_130N", "FDRq_130C", "FDRq_131N", "FDRq_131C",
  
  # FDRp por muestra (69-79)
  "FDRp_126t", "FDRp_127N", "FDRp_127C", "FDRp_128N", "FDRp_128C",
  "FDRp_129N", "FDRp_129C", "FDRp_130N", "FDRp_130C", "FDRp_131N", "FDRp_131C",
  
  # FDRpq por muestra (80-90)
  "FDRpq_126t", "FDRpq_127N", "FDRpq_127C", "FDRpq_128N", "FDRpq_128C",
  "FDRpq_129N", "FDRpq_129C", "FDRpq_130N", "FDRpq_130C", "FDRpq_131N", "FDRpq_131C"
)

# Verificar que quedó bien
glimpse(zpq_raw)


#3. LIMPIEZA Y SEPARACIÓN
#   is_phospho : el sufijo "_Phospho" en la secuencia identifica los
#                fosfopéptidos (fosforilación en serina, treonina o tirosina)
#   gen        : símbolo génico extraído del campo GN= de la descripción.
#                Se emplea porque es la nomenclatura estándar en la literatura
#                y en las bases de datos de interacción virus-hospedador.
#   ant_       : prefijo de la base de datos de contaminantes incorporada a la
#                búsqueda (proteínas séricas bovinas del medio de cultivo y
#                queratinas de contaminación ambiental). Se eliminan.


zpq_clean <- zpq_raw %>%
  mutate(across(Zq_126t:FDRpq_131C, as.numeric),
         is_phospho = str_detect(sequence, "_Phospho"),
         gen        = str_extract(description, "(?<=GN=)[^ ]+")) %>%
  filter(!str_detect(protein, "^ant_"))


#4. COMPROBACIONES

dim(zpq_clean) # Dimensiones del conjunto de datos

sum(!zpq_clean$is_phospho) # Péptidos no fosforilados

sum(zpq_clean$is_phospho) # Fosfopéptidos

n_distinct(zpq_clean$protein) #Proteínas únicas:

sum(is.na(zpq_clean$gen)) # Entradas sin símbolo génico (descripciones sin campo GN=).

# Fosfopéptidos sin valor de Zpq calculable.
# Solo disponen de valor de Zpq los péptidos cuya proteína fue identificada en el proteoma con más de un péptido.
zpq_clean %>%
  filter(is_phospho) %>%
  summarise(total       = n(),
            sin_Zpq     = sum(is.na(Zpq_AvA_dasatinib)),
            pct_sin_Zpq = round(100 * sum(is.na(Zpq_AvA_dasatinib)) / n(), 1)) %>%
  print()

# 5. GUARDAR

saveRDS(zpq_clean, "data_processed/zpq_clean.rds")

