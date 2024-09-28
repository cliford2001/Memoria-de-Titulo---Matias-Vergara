library(vcfR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

# Cargar el archivo VCF fusionado BRUTOSSSSSS
vcf1 <- read.vcfR("~/vep_data/VCFs/merged.vcf.gz")

# Extraer la matriz de genotipos
gt <- extract.gt(vcf1)

# Crear un dataframe con las posiciones y la matriz de genotipos
vcf1_df <- cbind(data.frame(CHROM=vcf1@fix[,1], POS=vcf1@fix[,2]), as.data.frame(gt))

# Crear un dataframe largo donde cada fila es una muestra
vcf1_df_long <- vcf1_df %>%
  pivot_longer(cols = c(A1, A2, A3, G1, G2, G3), 
               names_to = "Sample", 
               values_to = "Genotype")

# Cargar el archivo VCF fusionado FILTRO PASS
vcf2 <- read.vcfR("/home/administrador/vep_data/VCFs/anotados/Filtradas/merged_rediportal_filtered.vcf.gz")

# Extraer la matriz de genotipos
gt <- extract.gt(vcf2)

# Crear un dataframe con las posiciones y la matriz de genotipos
vcf2_df <- cbind(data.frame(CHROM=vcf2@fix[,1], POS=vcf2@fix[,2]), as.data.frame(gt))

# Crear un dataframe largo donde cada fila es una muestra
vcf2_df_long <- vcf2_df %>%
  pivot_longer(cols = c(A1, A2, A3, G1, G2, G3), 
               names_to = "Sample", 
               values_to = "Genotype")

# Cargar el archivo VCF fusionado
vcf3 <- read.vcfR("/home/administrador/vep_data/VCFs/anotados/Filtradas/VEP/SNV/merged_rediportal_filtered_annotated_SNV.vcf.gz")

# Extraer la matriz de genotipos
gt <- extract.gt(vcf3)

# Crear un dataframe con las posiciones y la matriz de genotipos
vcf3_df <- cbind(data.frame(CHROM=vcf3@fix[,1], POS=vcf3@fix[,2]), as.data.frame(gt))

# Crear un dataframe largo donde cada fila es una muestra
vcf3_df_long <- vcf3_df %>%
  pivot_longer(cols = c(A1, A2, A3, G1, G2, G3), 
               names_to = "Sample", 
               values_to = "Genotype")

# Cargar el archivo VCF fusionado
vcf4 <- read.vcfR("/home/administrador/vep_data/REDIPORTAL/Filtradas_SNPs/filtered_vcfs/merged_rediportal_filtered.vcf.gz")

# Extraer la matriz de genotipos
gt <- extract.gt(vcf4)

# Crear un dataframe con las posiciones y la matriz de genotipos
vcf4_df <- cbind(data.frame(CHROM=vcf4@fix[,1], POS=vcf4@fix[,2]), as.data.frame(gt))

# Crear un dataframe largo donde cada fila es una muestra
vcf4_df_long <- vcf4_df %>%
  pivot_longer(cols = c(SHADAR1, SHADAR2, SHC1, SHC2), 
               names_to = "Sample", 
               values_to = "Genotype")



# Definir colores pastel personalizados para cada muestra usando códigos hexadecimales
colors <- c("A1" = "#762a83",   # Pastel green
            "A2" = "#af8dc3",   # Light pastel green
            "A3" = "#e7d4e8",      # Light pastel purple
            "G1" = "#d9f0d3",
            "G2" = "#7fbf7b",
            "G3" = "#1b7837")      # Lavender



# Generar gráfico con colores personalizados
plot1 <-  ggplot(vcf1_df_long %>% filter(!is.na(Genotype)), aes(x = Sample, fill = Sample)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Distribución de variantes por muestra",
       x = "Muestra",
       y = "Número de variantes")

# Generar gráfico con colores personalizados
plot2 <- ggplot(vcf2_df_long %>% filter(!is.na(Genotype)), aes(x = Sample, fill = Sample)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Distribución de variantes por muestra",
       x = "Muestra",
       y = "Número de variantes")

# Generar gráfico con colores personalizados
plot3 <- ggplot(vcf3_df_long %>% filter(!is.na(Genotype)), aes(x = Sample, fill = Sample)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Distribución de variantes por muestra",
       x = "Muestra",
       y = "Número de variantes")

# Generar gráfico con colores personalizados
plot4 <- ggplot(vcf4_df_long %>% filter(!is.na(Genotype)), aes(x = Sample, fill = Sample)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Distribución de variantes por muestra",
       x = "Muestra",
       y = "Número de variantes")


# Combinar los tres gráficos en una sola imagen
titulo_grid <- textGrob("Cantidad de eventos de splicing brutos: 78707", gp = gpar(fontsize = 16, fontface = "bold"))

grid.arrange(
  titulo_grid, 
  arrangeGrob(plot1, plot2, plot3, plot4, ncol = 2),
  ncol = 1,
  heights = unit.c(unit(1, "lines"), unit(1, "npc") - unit(1, "lines"))
)

grid.arrange(
  titulo_grid, 
  arrangeGrob(plot1, plot2, plot3, ncol = 2),
  ncol = 1,
  heights = unit.c(unit(1, "lines"), unit(1, "npc") - unit(1, "lines"))
)




library(vcfR)
library(tidyverse)

# Función para cargar y procesar un archivo VCF
procesar_vcf <- function(file_path) {
  vcf <- read.vcfR(file_path)
  gt <- extract.gt(vcf)
  vcf_df <- cbind(data.frame(CHROM=vcf@fix[,1], POS=vcf@fix[,2]), as.data.frame(gt))
  
  # Crear un dataframe largo donde cada fila es una muestra
  vcf_df_long <- vcf_df %>%
    pivot_longer(cols = c(A1, A2, A3, G1, G2, G3), 
                 names_to = "Sample", 
                 values_to = "Genotype") %>%
    filter(!is.na(Genotype) & Genotype != "./.")  # Filtrar para eliminar valores missing y sin llamada de genotipo
  
  return(vcf_df_long)
}
vcf6_df_long <- procesar_vcf("/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/Merge/merged_filtered.vcf") # A_G T_C


# Procesar cada archivo VCF
vcf1_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/merged.vcf") # Brutos
vcf2_df_long <- procesar_vcf("/home/administrador/Escritorio/Variantes/Filtered/merged_Snp_PASS.vcf.gz") # Anotados por VEP #antoados por Rediportal / filtro PASS
vcf3_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/anotados/Filtradas/hola/VEP/merged_rediportal_filtered_annotated.vcf.gz") # Solamente 1:1
vcf4_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/anotados/Filtradas/hola/VEP/SNV/ADAR/merged_rediportal_filtered_annotated_SNV_ADAR.vcf.gz") # A_G T_C
vcf5_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/anotados/Filtradas/hola/VEP/SNV/ADAR/AD/merged_filtered.vcf.gz") # Filtro Por DP
vcf6_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/anotados/Filtradas/hola/VEP/SNV/ADAR/AD/COSMIC/merged_COSMIC.vcf.gz") 
vcf7_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/anotados/Filtradas/hola/VEP/SNV/ADAR/AD/COSMIC/GENO/final_filtered_COSMIC.vcf.gz") #GENO
vcf8_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/anotados/Filtradas/hola/VEP/SNV/ADAR/AD/COSMIC/GENO/merged_filtered_A_G_MAX_AF_0.01.vcf") # AF


# Procesar cada archivo VCF
vcf1_df_long <- procesar_vcf("/home/administrador/vep_data/VCFs/merged.vcf") # Brutos
vcf2_df_long <- procesar_vcf("/home/administrador/Escritorio/Variantes/Filtered/merged_Snp_PASS.vcf.gz") # Solamente 1:1
vcf3_df_long <- procesar_vcf("/home/administrador/Escritorio/Variantes/Filtered/ADAR/merged_A_G.vcf") # A_G 
vcf4_df_long <- procesar_vcf("/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/anotado_dbsbp.vcf.gz") # A_G T_C
vcf5_df_long <- procesar_vcf("/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/merge.vcf") # A_G T_C
vcf6_df_long <- procesar_vcf("/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/Merge/merged_filtered.vcf") # A_G T_C



contar_variantes <- function(df_long, filter_name) {
  df_long %>%
    group_by(Sample) %>%
    summarize(Count = n()) %>%
    mutate(Filter = filter_name)
}

# Generar conteos para cada archivo VCF
vcf1_variant_counts <- contar_variantes(vcf1_df_long, "Variantes Brutas")
vcf2_variant_counts <- contar_variantes(vcf2_df_long, "Filtro 1")
vcf3_variant_counts <- contar_variantes(vcf3_df_long, "Filtro 2")
vcf4_variant_counts <- contar_variantes(vcf4_df_long, "Filtro 3")
vcf5_variant_counts <- contar_variantes(vcf5_df_long, "Filtro 4")
vcf6_variant_counts <- contar_variantes(vcf6_df_long, "Ediciones Finales")
vcf7_variant_counts <- contar_variantes(vcf7_df_long, "Filtro 6")
vcf8_variant_counts <- contar_variantes(vcf8_df_long, "Variantes Finales")

# Definir el orden de los niveles en la columna 'Filter'
# Combinar todos los conteos en un solo dataframe
combined_variant_counts <- bind_rows(
  vcf1_variant_counts,
  vcf2_variant_counts,
  vcf3_variant_counts,
  vcf4_variant_counts,
  vcf5_variant_counts,
  vcf6_variant_counts,
  vcf7_variant_counts,
  vcf8_variant_counts
)


library(ggplot2)

# Crear un vector de colores para cada filtro
filtros_colores <- c("Variantes Brutas" = "#c7e9c0",  # Color para Variantes Brutas
                     "Filtro 1" = "#a1d99b",          # Color para Filtro 1
                     "Filtro 2" = "#74c476",          # Color para Filtro 2
                     "Filtro 3" = "#238b45",          # Color para Filtro 3
                     "Filtro 4" = "#005a32",          # Color para Filtro 4
                     "Ediciones Finales" = "#756bb1",          # Color para Filtro 5
                     "Filtro 6" = "#6a51a3",          # Color para Filtro 6
                     "Filtro 7" = "#fc9272",          # Color para Filtro 7
                     "Filtro 8" = "#de2d26",          # Color para Filtro 8
                     "Filtro 9" = "#756bb1",          # Color para Filtro 9
                     "Filtro 10" = "#fdae61",         # Color para Filtro 10
                     "Final" = "#d73027")  

# Asegurar que los filtros se ordenen en el orden deseado
combined_variant_counts$Filter <- factor(combined_variant_counts$Filter, 
                                         levels = names(filtros_colores))

# Crear el gráfico con colores para cada filtro
ggplot(combined_variant_counts, aes(x = Sample, y = Count, fill = Filter)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +  # Ajustar el ancho de las barras
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Conteo de Variantes por Muestra y Filtro", 
       x = "Muestra", 
       y = "Conteo de Variantes", 
       fill = "Filtro") +
  theme_minimal() +
  scale_fill_manual(values = filtros_colores) +  # Asignar colores a los filtros
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ajustar el eje Y para comenzar en 0
  theme(
    legend.position = "right",  # Posición de la leyenda: puede ser "left", "right", "bottom", "top"
    axis.line = element_line(color = "black"),  # Color negro para los ejes
    axis.title = element_text(color = "black"),  # Color negro para los títulos de los ejes
    axis.text = element_text(color = "black"),  # Color negro para los textos de los ejes
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Ajustar los márgenes del gráfico
  )


ggplot(combined_variant_counts, aes(x = Sample, y = Count, fill = Filter)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +  
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Conteo de Variantes por Muestra y Filtro", 
       x = "Muestra", 
       y = "Conteo de Variantes",  # Cambiar el título del eje Y
       fill = "Filtro") +
  theme_minimal() +
  scale_fill_manual(values = filtros_colores) +  
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks = NULL) +  # Ajustar el eje Y para comenzar en 0
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),  
    axis.title = element_text(color = "black"),  
    axis.text = element_text(color = "black"),  
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )



library(vcfR)
library(dplyr)

# Cargar el archivo VCF anotado
vcf <- read.vcfR("/home/administrador/vep_data/VCFs/anotados/Filtradas/VEP/SNV/ADAR/Cualite/AD/Rediportal/COSMIC/merged_filtered_COSMIC.vcf")
head(vcf@meta)
head(vcf@fix)
head(vcf@gt)









# Cargar las librerías necesarias
library(vcfR)
library(tidyverse)

# Cargar el archivo VCF fusionado
vcf <- read.vcfR("/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/anotado_dbsbp.vcf.gz")

# Extraer la información de profundidad de cobertura (DP)
dp <- extract.info(vcf, element = "DP")

# Crear un dataframe con las posiciones y la profundidad de cobertura
coverage_df <- data.frame(CHROM=vcf@fix[,1], POS=vcf@fix[,2], DP=as.numeric(dp))

# Filtrar las variantes con DP mayor a 10
coverage_df_above_10 <- coverage_df %>%
  filter(DP > 5)

# Filtrar las variantes con DP menor o igual a 10
coverage_df_below_or_equal_10 <- coverage_df %>%
  filter(DP <= 5)

# Contar cuántas variantes quedan en cada grupo
count_above_10 <- nrow(coverage_df_above_10)
count_below_or_equal_10 <- nrow(coverage_df_below_or_equal_10)

# Mostrar los resultados
cat("Número de variantes con DP mayor a 5:", count_above_10, "\n")
cat("Número de variantes con DP menor o igual a 5:", count_below_or_equal_10, "\n")

# Cargar las librerías necesarias
library(vcfR)
library(tidyverse)

# Cargar el archivo VCF fusionado
vcf <- read.vcfR("/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/filtered_output.vcf")

# Extraer la información de profundidad de cobertura (DP)
dp <- extract.info(vcf, element = "DP")

# Crear un dataframe con las posiciones y la profundidad de cobertura
coverage_df <- data.frame(CHROM=vcf@fix[,1], POS=vcf@fix[,2], DP=as.numeric(dp))

# Generar un histograma de las profundidades de cobertura (DP) con un límite superior de 50
ggplot(coverage_df, aes(x = DP)) +
  geom_histogram(binwidth = 1, fill = "#82c294FF", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribución de las profundidades de cobertura (DP)",
       x = "Profundidad de cobertura (DP)",
       y = "Número de variantes") +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ajustar el eje Y para comenzar en 0
  theme(
    axis.line = element_line(color = "black"),  # Color negro para los ejes
    axis.title = element_text(color = "black"),  # Color negro para los títulos de los ejes
    axis.text = element_text(color = "black"),  # Color negro para los textos de los ejes
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Ajustar los márgenes del gráfico
  )




# Cargar las librerías necesarias
library(vcfR)
library(dplyr)

# Cargar el archivo VCF
vcf <- read.vcfR("/home/administrador/vep_data/VCFs/merged.vcf.gz")

# Extraer la información de la profundidad de cobertura del alelo alternativo (AD)
ad <- extract.gt(vcf, element = "AD")

# El campo AD generalmente contiene dos valores separados por coma: ref_reads,alt_reads.
# Necesitamos extraer el segundo valor (alt_reads) para todas las variantes.

# Convertir el valor de AD en un dataframe
ad_df <- data.frame(
  CHROM = vcf@fix[,1],
  POS = vcf@fix[,2],
  ALT_COV = as.numeric(sapply(strsplit(ad, ","), `[`, 2))  # Extraer la cobertura del alelo alternativo
)


# Filtrar las variantes sin datos o con valores nulos
ad_df <- ad_df %>% filter(!is.na(ALT_COV))

# Generar un histograma de las profundidades de cobertura del alelo alternativo (ALT_COV) con un límite superior de 50
ggplot(ad_df, aes(x = ALT_COV)) +
  geom_histogram(binwidth = 1, fill = "#82c294FF", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribución de la profundidad de cobertura del alelo alternativo",
       x = "Profundidad de cobertura del alelo alternativo (ALT_COV)",
       y = "Número de variantes") +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ajustar el eje Y para comenzar en 0
  theme(
    axis.line = element_line(color = "black"),  # Color negro para los ejes
    axis.title = element_text(color = "black"),  # Color negro para los títulos de los ejes
    axis.text = element_text(color = "black"),  # Color negro para los textos de los ejes
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Ajustar los márgenes del gráfico
  )







# Cargar librerías necesarias
library(dplyr)
library(ggplot2)
library(readr)

# Definir la ruta base
base_path <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/rmats/SHADAR-SHC/rmats_post"

# Leer archivos de eventos de splicing
a3ss <- read_delim(file.path(base_path, "A3SS.MATS.JCEC.txt"), delim = "\t") %>%
  mutate(Tipo = "A3SS")
a5ss <- read_delim(file.path(base_path, "A5SS.MATS.JCEC.txt"), delim = "\t") %>%
  mutate(Tipo = "A5SS")
mxe <- read_delim(file.path(base_path, "MXE.MATS.JCEC.txt"), delim = "\t") %>%
  mutate(Tipo = "MXE")
ri <- read_delim(file.path(base_path, "RI.MATS.JCEC.txt"), delim = "\t") %>%
  mutate(Tipo = "RI")
se <- read_delim(file.path(base_path, "SE.MATS.JCEC.txt"), delim = "\t") %>%
  mutate(Tipo = "SE")

# Combinar los dataframes en uno solo
datos_jcec_comb <- bind_rows(a3ss, a5ss, mxe, ri, se)

# Verificar que se ha añadido la columna Tipo correctamente
print(head(datos_jcec_comb))

# Filtrar los datos para incluir solo los valores de p-values menores a 0.05, FDR < 0.1 y |PSI| > 0.8
datos_pvalue_filtrados_jcec <- datos_jcec_comb %>%
  rename(PSI = IncLevelDifference) %>%
  filter(PValue < 0.05 & FDR < 0.1 & abs(PSI) > 0.1)

# Exportar los datos filtrados a un archivo TSV
ruta_salida <- file.path(base_path, "datos_filtrados.tsv")
write_tsv(datos_pvalue_filtrados_jcec, ruta_salida)

# Confirmar la exportación
cat("Datos filtrados exportados a:", ruta_salida, "\n")

# Mostrar las primeras filas del dataframe para verificar el contenido
print(head(datos_pvalue_filtrados_jcec$Tipo))







# Cargar librerías necesarias
library(dplyr)
library(readr)
library(VariantAnnotation)
library(vcfR)
library(seqinr)
library(data.table)


# Cargar los datos de splicing filtrados
splicing_data_path <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/rmats/SHADAR-SHC/rmats_post/datos_filtrados.tsv"
splicing_data <- read_delim(splicing_data_path, delim = "\t")
splicing_data_subset <- splicing_data %>% dplyr::select(1, 2)
splicing_data_subset <- splicing_data_subset %>%
  dplyr::select(ID = ID...1, everything())

upEE <- splicing_data[, c("chr", "upstreamEE", "upstreamEE","geneSymbol", "ID...1" ,"strand")]
star <- splicing_data[, c("chr", "exonStart_0base","exonStart_0base","geneSymbol","ID...1" ,"strand")]
end <- splicing_data[, c("chr",  "exonEnd","exonEnd","geneSymbol","ID...1" ,"strand")]
downEE <- splicing_data[, c("chr", "downstreamES", "downstreamES" ,"geneSymbol","ID...1" ,"strand")]
# Definir la ruta de salida
output_dir <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/"

# Guardar los datos en nuevos archivos
write_delim(upEE, file.path(output_dir, "upEE.tsv"), delim = "\t")
write_delim(star, file.path(output_dir, "star.tsv"), delim = "\t")
write_delim(end, file.path(output_dir, "end.tsv"), delim = "\t")
write_delim(downEE, file.path(output_dir, "downEE.tsv"), delim = "\t")

ventanaA <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/ventanaA.tsv"
ventanaA <- read_delim(ventanaA, delim = "\t", col_names = FALSE)
ventanaB <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/ventanaB.tsv"
ventanaB <- read_delim(ventanaB, delim = "\t", col_names = FALSE)
ventanaB <- ventanaB %>% dplyr::select(-c(1, 4, 6))
ventanaC <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/ventanaC.tsv"
ventanaC <- read_delim(ventanaC, delim = "\t", col_names = FALSE)
ventanaC <- ventanaC %>% dplyr::select(-c(1, 4, 6))
ventanaD <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/ventanaD.tsv"
ventanaD <- read_delim(ventanaD, delim = "\t", col_names = FALSE)
ventanaD <- ventanaD %>% dplyr::select(-c(1, 4, 6))

ventanaS <- left_join(ventanaA  , ventanaD, by = "X5")
# Reorganizar las columnas y agregar encabezado
ventanaS_ordenado <- ventanaS %>%
  dplyr::select(ID = X5, cromosomas = X1, gen_ID = X4, strand = X6, everything())
# Renombrar las columnas en el dataset ventanaS_ordenado
ventanaS_ordenado <- ventanaS_ordenado %>%
  dplyr::rename(
    Start_ventanaA = X2.x,
    End_ventanaA = X3.x,
    Start_ventanaB = X2.y,
    End_ventanaB = X3.y,
    Start_ventanaC = X2.x.x,
    End_ventanaC = X3.x.x,
    Start_ventanaD = X2.y.y,
    End_ventanaD = X3.y.y
  )
# Mostrar el resultado
head(ventanaS_ordenado)
write_delim(ventanaS_ordenado, file.path(output_dir, "Ventanas.tsv"), delim = "\t")


##################################################################################################################################################################
ventanaS <- left_join(splicing_data  , splicing_data_subset, by = "ID")
# Supongamos que tu dataset se llama ventanaS_ordenado
ventanaS_ordenado <- ventanaS %>%
  dplyr::select(ID, cromosomas, GeneID, gen_ID, strand, Start_ventanaA, End_ventanaA, Start_ventanaB, End_ventanaB, Start_ventanaC, End_ventanaC, Start_ventanaD, End_ventanaD)


write_delim(ventanaS_ordenado, file.path(output_dir, "Ventanas.tsv"), delim = "\t")

#######################################################
# Cargar los datos de splicing filtrados
splicing_data_path <- "/home/administrador/Escritorio/tesis/RESULTADOS/RNASPLICE/Ventanas.tsv"
splicing_data <- read_delim(splicing_data_path, delim = "\t")

vcf_path <- "/home/administrador/vep_data/SHADAR1_rediportal_filtered_final_filtered_annotated.vcf"
vcf <- read.vcfR(vcf_path, verbose = FALSE)
# Convertir a un dataframe
vcf_df <- as.data.frame(vcf@fix)












# Definir la ruta de salida del archivo VCF
vcf_output_path <- "/home/administrador/vep_data/matches_output.vcf"

# Escribir el archivo VCF
write.vcf(vcf_new, file = vcf_output_path)

# Confirmar la exportación
cat("VCF con coincidencias exportado a:", vcf_output_path, "\n")

info <- matched_variants$INFO
# Convertir la matriz resultante en un DataFrame
info_df <- as_tibble(info)
# Crear una nueva columna separando desde el primer "|"
info_split <- info_df %>%
  mutate(Part2 = str_extract(value, "(?<=\\|).*")) %>%
  dplyr::select(Part2)  # Seleccionar solo la columna Part2

matched_variants_with_part2 <- matched_variants %>%
  bind_cols(info_split)

# Mostrar las primeras filas para verificar el resultado
print(head(info_split))
info_split_separated <- info_split %>%
  separate(Part2, into = paste0("Segment_", 1:10000), sep = "\\|", fill = "right")













library(vcfR)
library(tidyverse)

vcf <- read.vcfR("/home/administrador/vep_data/VCFs/anotados/output_filtered.vcf")
vcf_df <- as.data.frame(vcf@fix)
info_column <- vcf_df$INFO
head(info_column)

# Extraer solo el campo CSQ de cada elemento de la columna
csq_info <- sapply(info_column, function(x) {
  # Buscar la posición del campo CSQ dentro de cada cadena
  matches <- regmatches(x, regexpr("CSQ=[^;]*", x))
  
  # Retornar el contenido de CSQ sin el prefijo "CSQ="
  if(length(matches) > 0) {
    return(sub("CSQ=", "", matches))
  } else {
    return(NA)
  }
})

# Convertir en un dataframe
csq_df <- as.data.frame(csq_info, stringsAsFactors = FALSE)

# Crear un dataframe a partir de la columna 'csq_info'
csq_df <- as.data.frame(csq_info, stringsAsFactors = FALSE)

# Separar las columnas de 'csq_info' usando el delimitador "|"
csq_df_separado <- do.call(rbind, strsplit(csq_df$csq_info, "\\|"))

# Convertir el resultado a un data.frame
csq_df_separado <- as.data.frame(csq_df_separado, stringsAsFactors = FALSE)

# Mostrar el dataframe resultante
print(csq_df_separado)



# Suponiendo que csq_df_separado es el dataframe resultante del paso anterior y V4 es la columna con los nombres de los genes

# Extraer los nombres de los genes de la columna V4
genes <- csq_df_separado$V4

# Contar la frecuencia de cada gen
gene_counts <- table(genes)

# Convertir la tabla de frecuencias en un dataframe
gene_counts_df <- as.data.frame(gene_counts)

# Renombrar las columnas para mayor claridad
colnames(gene_counts_df) <- c("Gene", "Ediciones")

# Mostrar el dataframe resultante
print(gene_counts_df)



dea <- read.csv("/home/administrador/Escritorio/rnaseq/DEA_DTU_DEU_Splice.csv")
dea <- dea[ , -1]

# Verificar el resultado
head(dea)
head(gene_counts_df)


# Filtrar las filas sin nombre de gen
gene_counts_df <- gene_counts_df[gene_counts_df$Gene != "", ]

# Renombrar la columna 'Gene' en gene_counts_df para que coincida con 'gene_symbol' en dea
colnames(gene_counts_df)[colnames(gene_counts_df) == "Gene"] <- "gene_symbol"

# Realizar el left join entre dea y gene_counts_df usando la columna 'gene_symbol'
merged_df <- merge(dea, gene_counts_df, by = "gene_symbol", all.x = TRUE)

# Convertir los valores NA en 0 en la columna 'Ediciones'
merged_df$Ediciones[is.na(merged_df$Ediciones)] <- 0

# Verificar el resultado
head(merged_df)
write.csv(merged_df, "/home/administrador/Escritorio/rnaseq/Tabla_maestra_2.csv")





library(vcfR)
library(tidyverse)

vcf <- read.vcfR("/home/administrador/vep_data/VCFs/anotados/Solo_A_annotated.vcf.gz")
vcf_df <- as.data.frame(vcf@fix)
info_column <- vcf_df$INFO
head(info_column)

# Extraer solo el campo CSQ de cada elemento de la columna
csq_info <- sapply(info_column, function(x) {
  # Buscar la posición del campo CSQ dentro de cada cadena
  matches <- regmatches(x, regexpr("CSQ=[^;]*", x))
  
  # Retornar el contenido de CSQ sin el prefijo "CSQ="
  if(length(matches) > 0) {
    return(sub("CSQ=", "", matches))
  } else {
    return(NA)
  }
})

# Convertir en un dataframe
csq_df <- as.data.frame(csq_info, stringsAsFactors = FALSE)

# Crear un dataframe a partir de la columna 'csq_info'
csq_df <- as.data.frame(csq_info, stringsAsFactors = FALSE)

# Separar las columnas de 'csq_info' usando el delimitador "|"
csq_df_separado <- do.call(rbind, strsplit(csq_df$csq_info, "\\|"))

# Convertir el resultado a un data.frame
csq_df_separado <- as.data.frame(csq_df_separado, stringsAsFactors = FALSE)

# Mostrar el dataframe resultante
print(csq_df_separado)

# Suponiendo que csq_df_separado es el dataframe resultante del paso anterior y V4 es la columna con los nombres de los genes

# Extraer los nombres de los genes de la columna V4
genes <- csq_df_separado$V4

# Contar la frecuencia de cada gen
gene_counts <- table(genes)

# Convertir la tabla de frecuencias en un dataframe
gene_counts_df <- as.data.frame(gene_counts)

# Renombrar las columnas para mayor claridad
colnames(gene_counts_df) <- c("Gene", "Ediciones")

# Mostrar el dataframe resultante
print(gene_counts_df)

dea <- read.csv("/home/administrador/Escritorio/rnaseq/Tabla_maestra.csv")
dea <- dea[ , -1]

# Verificar el resultado
head(dea)
head(gene_counts_df)


# Filtrar las filas sin nombre de gen
gene_counts_df <- gene_counts_df[gene_counts_df$Gene != "", ]

# Renombrar la columna 'Gene' en gene_counts_df para que coincida con 'gene_symbol' en dea
colnames(gene_counts_df)[colnames(gene_counts_df) == "Gene"] <- "gene_symbol"

# Realizar el left join entre dea y gene_counts_df usando la columna 'gene_symbol'
merged_df <- merge(dea, gene_counts_df, by = "gene_symbol", all.x = TRUE)

# Convertir los valores NA en 0 en la columna 'Ediciones'
merged_df$Cruce[is.na(merged_df$Cruce)] <- 0

# Verificar el resultado
head(merged_df)
write.csv(merged_df, "/home/administrador/Escritorio/rnaseq/Tabla_maestra_2.csv")







library(vcfR)

# Cargar el archivo VCF
vcf <- read.vcfR("/home/administrador/vep_data/VCFs/anotados/Solo_A.vcf")

# Extraer la columna INFO del VCF
info_column <- vcf@fix[,"INFO"]

# Identificar si las variantes están en dbSNP, Rediportal, o COSMIC
en_dbSNP <- grepl("rs", info_column, ignore.case = TRUE)
en_Rediportal <- grepl("Rediportal=T", info_column, ignore.case = TRUE)
en_COSMIC <- grepl("COSMIC", info_column, ignore.case = TRUE)

# Contar variantes conocidas y noveles
n_conocidas <- sum(en_dbSNP | en_Rediportal | en_COSMIC)
n_noveles <- sum(!(en_dbSNP | en_Rediportal | en_COSMIC))

# Mostrar resultados
cat("Número de variantes conocidas: ", n_conocidas, "\n")
cat("Número de variantes noveles: ", n_noveles, "\n")

# Contar cuántas variantes están en Rediportal
n_rediportal <- sum(en_Rediportal)
n_no_rediportal <- sum(!en_Rediportal)

# Mostrar resultados
cat("Número de variantes en Rediportal: ", n_rediportal, "\n")
cat("Número de variantes NO en Rediportal: ", n_no_rediportal, "\n")


# Extraer el valor de QUAL
qual_values <- as.numeric(vcf@fix[,"QUAL"])

# Resumen estadístico de la calidad
summary(qual_values)

# Histograma de la calidad de las variantes
hist(qual_values, breaks = 50, main = "Distribución de la Calidad de las Variantes", xlab = "QUAL", col = "lightgreen")



# Extraer la información de GT (genotipo) de cada muestra
genotipos <- extract.gt(vcf)

# Ver el head de los genotipos
head(genotipos)





library(ggplot2)
library(tidyverse)

# Datos
variantes <- data.frame(
  Tipo = c("Conocidas", "Noveles"),
  Cantidad = c(6264, 15647)
)

# Gráfico de barras
ggplot(variantes, aes(x = Tipo, y = Cantidad, fill = Tipo)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribución de Variantes Conocidas vs. Noveles", x = "Tipo de Variante", y = "Cantidad") +
  theme_minimal()

# Gráfico de pastel
ggplot(variantes, aes(x = "", y = Cantidad, fill = Tipo)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Distribución de Variantes Conocidas vs. Noveles") +
  theme_minimal()

# Datos
rediportal <- data.frame(
  Tipo = c("En Rediportal", "No en Rediportal"),
  Cantidad = c(5772, 16139)
)

# Gráfico de pastel
ggplot(rediportal, aes(x = "", y = Cantidad, fill = Tipo)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Distribución de Variantes en Rediportal vs. No en Rediportal") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_text(aes(label = paste0(Tipo, "\n", Cantidad, " (", round(Cantidad/sum(Cantidad) * 100, 1), "%)")), 
            position = position_stack(vjust = 0.5))

library(ggplot2)

# Datos
consecuencias <- data.frame(
  Consecuencia = c("splice_acceptor_variant", "splice_donor_variant", "stop_lost", 
                   "start_lost", "missense_variant", "splice_donor_5th_base_variant", 
                   "splice_region_variant", "splice_donor_region_variant", 
                   "splice_polypyrimidine_tract_variant", "stop_retained_variant", 
                   "synonymous_variant", "mature_miRNA_variant", "5_prime_UTR_variant", 
                   "3_prime_UTR_variant", "non_coding_transcript_exon_variant", 
                   "intron_variant", "upstream_gene_variant", "downstream_gene_variant", 
                   "regulatory_region_variant"),
  Cantidad = c(22, 16, 24, 2, 532, 4, 82, 25, 71, 4, 545, 2, 361, 9024, 11139, 35, 17, 4, 2)
)

# Calcular porcentajes
consecuencias$Porcentaje <- round(consecuencias$Cantidad / sum(consecuencias$Cantidad) * 100, 1)
consecuencias$Etiqueta <- paste0(consecuencias$Consecuencia, " (", consecuencias$Porcentaje, "%)")

# Paleta de colores
colores <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
             "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", 
             "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f")

# Gráfico de pastel
ggplot(consecuencias, aes(x = "", y = Cantidad, fill = Etiqueta)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colores) +
  labs(title = "Distribución de las Consecuencias Funcionales de las Variantes") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = "Consecuencia"))



# Gráfico de pastel
library(ggplot2)
ggplot(consecuencias, aes(x = "", y = Cantidad, fill = Etiqueta)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colores) +
  labs(title = "Distribución de las Consecuencias Funcionales de las Variantes") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)  # Ajuste de tamaño de la leyenda a 14
  ) +
  guides(fill = guide_legend(title = "Consecuencia"))








library(ggplot2)
library(tidyverse)

# Datos
variantes <- data.frame(
  Tipo = c("Conocidas", "Noveles"),
  Cantidad = c(6264, 15647)
)

# Gráfico de barras
ggplot(variantes, aes(x = Tipo, y = Cantidad, fill = Tipo)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribución de Variantes Conocidas vs. Noveles", x = "Tipo de Variante", y = "Cantidad") +
  theme_minimal()

# Gráfico de pastel
ggplot(variantes, aes(x = "", y = Cantidad, fill = Tipo)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Distribución de Variantes Conocidas vs. Noveles") +
  theme_minimal()

# Datos
rediportal <- data.frame(
  Tipo = c("En Rediportal", "No en Rediportal"),
  Cantidad = c(5772, 16139)
)

# Gráfico de pastel
ggplot(rediportal, aes(x = "", y = Cantidad, fill = Tipo)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Distribución de Variantes en Rediportal vs. No en Rediportal") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_text(aes(label = paste0(Tipo, "\n", Cantidad, " (", round(Cantidad/sum(Cantidad) * 100, 1), "%)")), 
            position = position_stack(vjust = 0.5))








library(ggplot2)

# Datos
consecuencias <- data.frame(
  Consecuencia = c("missense_variant", 
                   "splice_polypyrimidine_tract_variant", 
                   "synonymous_variant", 
                   "5_prime_UTR_variant", 
                   "3_prime_UTR_variant", 
                   "non_coding_transcript_exon_variant"),
  Cantidad = c(5, 2, 2, 19, 10, 32)
)

# Calcular porcentajes
consecuencias$Porcentaje <- round(consecuencias$Cantidad / sum(consecuencias$Cantidad) * 100, 1)
consecuencias$Etiqueta <- paste0(consecuencias$Consecuencia, " (", consecuencias$Porcentaje, "%)")

# Paleta de colores
colores <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")

# Gráfico de pastel
ggplot(consecuencias, aes(x = "", y = Cantidad, fill = Etiqueta)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colores) +
  labs(title = "Distribución de las Consecuencias Funcionales de las Variantes") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  ) +
  guides(fill = guide_legend(title = "Consecuencia"))














