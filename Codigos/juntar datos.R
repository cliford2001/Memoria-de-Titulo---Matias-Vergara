library(edgeR)
library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(VennDiagram)


# Leer el archivo tx2gene
tx2gene <- read_tsv("/home/administrador/Escritorio/rnaseq/tx2gene.tsv", col_names = FALSE)
colnames(tx2gene) <- c("transcript_id", "gene_id", "gene_symbol")
tx2gene_symbols_unique <- tx2gene %>%
  dplyr::select(gene_id, gene_symbol) %>%
  distinct(gene_id, .keep_all = TRUE)

tx2gene_symbols_unique_t <- tx2gene %>%
  dplyr::select(transcript_id, gene_symbol) %>%
  distinct(transcript_id, .keep_all = TRUE)

sig_dtu <- read.csv("/home/administrador/Escritorio/rnaseq/DTU_SIG.csv")
sig_deu <- read.csv("/home/administrador/Escritorio/rnaseq/SIG_DEU.csv")
sig_dea <- read.csv("/home/administrador/Escritorio/rnaseq/SIG_DEA.csv")
sig_dea <- sig_dea[, -c(1:5)]
sig_splicing <- read.csv("/home/administrador/Escritorio/rnaseq/SIG_SPLICING.csv")

dtu <- read.csv("/home/administrador/Escritorio/rnaseq/DTU.csv")
deu <- read.csv("/home/administrador/Escritorio/rnaseq/contrast_A_G.usage.exon.csv")
dea <- read.csv("/home/administrador/Escritorio/rnaseq/DEA_no_pseudo.csv")
x
# Usar left_join para agregar la columna gene_symbol a sig_deu
sig_deu <- sig_deu %>%
  left_join(tx2gene_symbols_unique, by = c("Geneid" = "gene_id"))
deu <- deu %>%
  left_join(tx2gene_symbols_unique, by = c("Geneid" = "gene_id"))
dtu <- dtu %>%
  left_join(tx2gene_symbols_unique_t, by = c("X" = "transcript_id"))

# Extraer los símbolos de genes de cada conjunto de datos
genes_dtu <- sig_dtu$gene_symbol
genes_dea <- sig_dea$Name_genes  # Suponiendo que la columna de símbolos en sig_dea es Name_genes
genes_deu <- sig_deu$gene_symbol
genes_splicing <- sig_splicing$geneSymbol

write.csv(sig_deu, "/home/administrador/Escritorio/rnaseq/SIG_DEU.csv")
write.csv(genes_dtu, "genes_dtu.txt", row.names = FALSE, quote  = FALSE)
write.csv(genes_dea, "genes_dea.txt", row.names = FALSE, quote  = FALSE)
write.csv(genes_deu, "genes_deu.txt", row.names = FALSE, quote  = FALSE)
write.csv(genes_splicing, "genes_splicing.txt", row.names = FALSE, quote  = FALSE)


# Crear una lista con los genes de cada análisis
gene_lists <- list(DTU = genes_dtu, DEA = genes_dea, DEU = genes_deu)

# Generar el diagrama de Venn con colores
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("DTU", "DEA", "DEU"),
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  filename = NULL, # Para generar el gráfico en el entorno de R y no como un archivo
  output = TRUE
)

# Mostrar el diagrama de Venn
grid.draw(venn.plot)

# Encontrar los genes que son comunes a los tres análisis
common_genes <- intersect(intersect(genes_dtu, genes_dea), genes_deu)

# Mostrar los genes comunes
length(common_genes)
common_genes
write.csv(common_genes,"genes_todos.txt", row.names = FALSE, quote  = FALSE)


# Diagrama VENN sub y sobre expresados.


# Crear una lista con los genes de cada análisis
gene_lists <- list(DTU = genes_dtu, DEA = genes_dea, DEU = genes_deu, SPLICING = genes_splicing)

# Generar el diagrama de Venn con colores
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("DTU", "DEA", "DEU","SPLICING"),
  fill = c("red", "green", "blue", "purple"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  filename = NULL, # Para generar el gráfico en el entorno de R y no como un archivo
  output = TRUE
)

# Mostrar el diagrama de Venn
grid.draw(venn.plot)

# Encontrar los genes que son comunes a los tres análisis
common_genes <- intersect(intersect(intersect(genes_dtu, genes_dea), genes_deu), genes_splicing)

# Mostrar los genes comunes
length(common_genes)
common_genes















# Contar el número de transcritos asociados a cada gen en el conjunto DTU
transcripts_per_gene_dtu <- dtu %>%
  group_by(gene_symbol) %>%
  summarise(n_transcripts = n())

# Visualizar los resultados
head(transcripts_per_gene_dtu)

# Contar el número de transcritos diferenciales en DTU (usando FDR < 0.05 como criterio de significancia)
differential_transcripts_dtu <- sig_dtu %>%
  filter(FDR < 0.05) %>%
  group_by(gene_symbol) %>%
  summarise(n_differential_transcripts = n())
head(differential_transcripts_dtu)

# Filtrar los transcritos diferenciales que están sobreexpresados (logFC > 0)
overexpressed_transcripts_dtu <- sig_dtu %>%
  filter(FDR < 0.05, logFC > 0) %>%
  group_by(gene_symbol) %>%
  summarise(n_overexpressed_transcripts = n())

# Filtrar los transcritos diferenciales que están subexpresados (logFC < 0)
underexpressed_transcripts_dtu <- sig_dtu %>%
  filter(FDR < 0.05, logFC < 0) %>%
  group_by(gene_symbol) %>%
  summarise(n_underexpressed_transcripts = n())

# Visualizar las primeras filas de ambos resultados
head(overexpressed_transcripts_dtu)
head(underexpressed_transcripts_dtu)


library(dplyr)

# Combinar los resultados de transcritos totales y transcritos diferenciales
combined_dtu <- left_join(transcripts_per_gene_dtu, differential_transcripts_dtu, by = "gene_symbol")

# Combinar con los sobreexpresados
combined_dtu <- left_join(combined_dtu, overexpressed_transcripts_dtu, by = "gene_symbol")

# Combinar con los subexpresados
combined_dtu <- left_join(combined_dtu, underexpressed_transcripts_dtu, by = "gene_symbol")

# Rellenar los NA con 0 para las columnas que faltan valores
combined_dtu <- combined_dtu %>%
  mutate(
    n_differential_transcripts = ifelse(is.na(n_differential_transcripts), 0, n_differential_transcripts),
    n_overexpressed_transcripts = ifelse(is.na(n_overexpressed_transcripts), 0, n_overexpressed_transcripts),
    n_underexpressed_transcripts = ifelse(is.na(n_underexpressed_transcripts), 0, n_underexpressed_transcripts)
  )

# Visualizar el resultado final
head(combined_dtu)

head(sig_dea)

# Realizar el merge entre combined_dtu y sig_dea, conservando solo las filas que están en ambas
final_combined_data <- merge(combined_dtu, sig_dea, by.x = "gene_symbol", by.y = "Name_genes", all = FALSE)

# Visualizar las primeras filas del resultado final
head(final_combined_data)

library(dplyr)


# Visualizar las primeras filas para confirmar la eliminación
head(final_combined_data)



library(dplyr)




# Realiza
library(dplyr)

# Reemplazar NA por 0 en todas las columnas que provienen de combined_deu
final_combined_data_with_deu <- final_combined_data_with_deu %>%
  mutate(across(starts_with("n_"), ~ ifelse(is.na(.), 0, .)))

# Visualizar las primeras filas del resultado final
head(final_combined_data_with_deu)





# Contar el número de exones asociados a cada gen en el conjunto DEU
exons_per_gene_deu <- deu %>%
  group_by(gene_symbol) %>%
  summarise(n_exons = n())

# Visualizar los resultados
head(exons_per_gene_deu)

# Contar el número de exones diferenciales en DEU (usando FDR < 0.05 como criterio de significancia)
differential_exons_deu <- sig_deu %>%
  filter(FDR < 0.05) %>%
  group_by(gene_symbol) %>%
  summarise(n_differential_exons = n())

# Visualizar los resultados
head(differential_exons_deu)

# Filtrar los exones diferenciales que están sobreexpresados (logFC > 0)
overexpressed_exons_deu <- sig_deu %>%
  filter(FDR < 0.05, logFC > 0) %>%
  group_by(gene_symbol) %>%
  summarise(n_overexpressed_exons = n())

# Visualizar los resultados
head(overexpressed_exons_deu)
# Filtrar los exones diferenciales que están subexpresados (logFC < 0)
underexpressed_exons_deu <- sig_deu %>%
  filter(FDR < 0.05, logFC < 0) %>%
  group_by(gene_symbol) %>%
  summarise(n_underexpressed_exons = n())

# Visualizar los resultados
head(underexpressed_exons_deu)

genes_SU <- overexpressed_events$geneSymbol  
genes_SO <- underexpressed_events$geneSymbol  

write.csv(sig_deu, "/home/administrador/Escritorio/rnaseq/SIG_DEU.csv")
write.csv(genes_dtu, "genes_dtu.txt", row.names = FALSE, quote  = FALSE)

# Combinar el número total de exones y el número de exones diferenciales
combined_deu <- left_join(exons_per_gene_deu, differential_exons_deu, by = "gene_symbol")

# Combinar con los exones sobreexpresados
combined_deu <- left_join(combined_deu, overexpressed_exons_deu, by = "gene_symbol")

# Combinar con los exones subexpresados
combined_deu <- left_join(combined_deu, underexpressed_exons_deu, by = "gene_symbol")

# Rellenar los NA con 0 para las columnas que faltan valores
combined_deu <- combined_deu %>%
  mutate(
    n_differential_exons = ifelse(is.na(n_differential_exons), 0, n_differential_exons),
    n_overexpressed_exons = ifelse(is.na(n_overexpressed_exons), 0, n_overexpressed_exons),
    n_underexpressed_exons = ifelse(is.na(n_underexpressed_exons), 0, n_underexpressed_exons)
  )

# Visualizar el resultado final
head(combined_deu)

final_combined_data_with_deu <- left_join(final_combined_data, combined_deu, by = "gene_symbol")
# Reemplazar NA por 0 en todas las columnas que provienen de combined_deu
final_combined_data_with_deu <- final_combined_data_with_deu %>%
  mutate(across(starts_with("n_"), ~ ifelse(is.na(.), 0, .)))

# Visualizar las primeras filas del resultado final
head(final_combined_data_with_deu)

# Visualizar las primeras filas del resultado final
head(final_combined_data)


# # Reorganizar las columnas dejando la información de DEA primero, seguida por la información de transcritos y exones
final_combined_data_reordered <- final_combined_data_with_deu %>%
  dplyr::select(
    gene_symbol,         # Mantener el identificador del gen al principio
    ensembl_id,          # Colocar la información de DEA inmediatamente después
    baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,  # Información de DEA
    n_transcripts, n_differential_transcripts, n_overexpressed_transcripts, n_underexpressed_transcripts, # Luego la info de transcritos
    n_exons, n_differential_exons, n_overexpressed_exons, n_underexpressed_exons                          # Finalmente la info de exones
  )

# Verifica que la reorganización se haya realizado correctamente
head(final_combined_data_reordered)


# Asignar un nombre al archivo CSV
output_file <- "/home/administrador/Escritorio/rnaseq/DEA_DTU_DEU_combined_data.csv"
# Agregar un nombre a la primera fila del CSV
write.csv(
  final_combined_data_reordered, 
  file = output_file, 
  row.names = FALSE
)

# Verificar la primera fila del resultado exportado para asegurarse de que todo esté correcto
final_combined_data_reordered <- read.csv(output_file)
head(final_combined_data_reordered)












# Leer la tabla
sig_splicing <- read.csv("/home/administrador/Escritorio/rnaseq/SIG_SPLICING.csv")
splicing <- read.csv("/home/administrador/Escritorio/rnaseq/SPLICING.csv")






# Contar el número de eventos de splicing por gen
events_per_gene <- splicing %>%
  group_by(geneSymbol) %>%
  summarise(n_events = n())

# Visualizar los resultados
head(events_per_gene)
sum(events_per_gene$n_events, na.rm = TRUE)


# Contar el número de eventos de splicing diferenciales (FDR < 0.05)
differential_events <- sig_splicing %>%
  group_by(geneSymbol) %>%
  summarise(n_differential_events = n())
sum(differential_events$n_differential_events, na.rm = TRUE)
# Visualizar los resultados
head(differential_events)

# Filtrar los eventos de splicing sobreexpresados (IncLevelDifference > 0)
overexpressed_events <- sig_splicing %>%
  filter(IncLevelDifference < 0)
sum(overexpressed_events$n_overexpressed_events, na.rm = TRUE)

# Visualizar los resultados
head(overexpressed_events)

# Filtrar los eventos de splicing subexpresados (IncLevelDifference < 0)
underexpressed_events <- sig_splicing %>%
  filter(IncLevelDifference > 0)
sum(underexpressed_events$n_underexpressed_events, na.rm = TRUE)

# Visualizar los resultados
head(underexpressed_events)



genes_SO <- overexpressed_events$geneSymbol  
genes_SU <- underexpressed_events$geneSymbol  
write.csv(genes_SU, "genes_SU.txt", row.names = FALSE, quote  = FALSE)
write.csv(genes_SO, "genes_SO.txt", row.names = FALSE, quote  = FALSE)

#ontología Splciing overexpresed
library(igraph)

# Cargar el dataframe desde el archivo
go_mf <- read.table("/home/administrador/Escritorio/ontologias/splicing/GO_Molecular_Function_2023_table.txt", 
                    header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Lista de términos de ontología de interés relacionados con ADAR
interesantes <- c(
  "Nuclear Receptor Binding",
  "DNA-directed RNA Polymerase Activity",
  "RNA Polymerase Activity",
  "Poly(A) Binding",
  "Histone Methyltransferase Activity",
  "GTPase Regulator Activity",
  "DNA-binding Transcription Activator Activity",
  "Transcription Coregulator Binding"
)

# Crear un patrón de búsqueda a partir de la lista de términos
pattern <- paste(interesantes, collapse = "|")

# Filtrar los términos en el dataframe usando el patrón
go_mf_filtrado <- go_mf %>%
  filter(grepl(pattern, Term, ignore.case = TRUE))

# Crear una lista de términos de interés específicos
interesantes <- c(
  "Double-Stranded RNA Binding",
  "RNA Exonuclease Activity, Producing 5'-Phosphomonoesters",
  "Exonuclease Activity",
  "GTPase Regulator Activity",
  "GTPase Activator Activity",
  "Phosphatidylinositol Binding",
  "3'-5'-RNA Exonuclease Activity",
  "3'-5'-DNA Exonuclease Activity"
)

# Filtrar los términos en el dataframe usando coincidencias parciales

pattern <- paste(interesantes, collapse = "|")
go_mf_filtrado <- go_mf %>%
  filter(`P.value` < 0.05) %>%
  filter(grepl(pattern, Term, ignore.case = TRUE)) 


ggplot(go_mf_filtrado, aes(x = Odds.Ratio, y = reorder(Term, Odds.Ratio))) +
  geom_point(aes(size = -log10(Adjusted.P.value), color = P.value)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "Dot Plot de Ontología",
    x = "Odds Ratio",
    y = "Término de Ontología"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),    # Ejes en color negro
    axis.text = element_text(color = "black", size = 14),    # Texto de los ejes en negro y mayor tamaño
    axis.title = element_text(color = "black", size = 16),   # Títulos de los ejes en negro y mayor tamaño
    plot.title = element_text(hjust = 0.5, size = 18),       # Centrar el título y aumentar tamaño
    legend.text = element_text(size = 12),                   # Tamaño de la leyenda
    legend.title = element_text(size = 14)                   # Tamaño del título de la leyenda
  )

sig_splicing <- read.csv("/home/administrador/Escritorio/rnaseq/SIG_SPLICING.csv")
# Filtrar los eventos con una diferencia significativa en la inclusión
sig_splicing_filtered <- sig_splicing %>%
  filter(abs(IncLevelDifference) >= 0.1)  # Cambia el umbral según lo que consideres significativo
# Crear una nueva columna combinando geneSymbol y ID para crear un identificador único
sig_splicing_filtered <- sig_splicing_filtered %>%
  mutate(eventID = paste(geneSymbol, ID...1, sep = "_"))
# Ordenar los eventos por PValue (o FDR si prefieres) de manera ascendente (los más pequeños primero)
sig_splicing_sorted <- sig_splicing_filtered %>%
  arrange(FDR)
# Seleccionar los primeros 40 eventos más significativos
top_40_splicing <- sig_splicing_sorted %>%
  head(40)
# Crear la matriz de PSI usando eventID para asegurar unicidad
psi_matrix_top40 <- top_40_splicing %>%
  select(eventID, PSI_MDA_A = IncLevel1, PSI_MDA_G = IncLevel2) %>%
  # Convertir las columnas de PSI a valores numéricos, tomando el promedio si hay varios valores
  mutate(
    PSI_MDA_A = sapply(strsplit(PSI_MDA_A, ","), function(x) mean(as.numeric(x), na.rm = TRUE)),
    PSI_MDA_G = sapply(strsplit(PSI_MDA_G, ","), function(x) mean(as.numeric(x), na.rm = TRUE))
  ) %>%
  column_to_rownames("eventID") %>%
  as.matrix()


# Generar el heatmap con la nueva gama de colores
pheatmap(
  psi_matrix_top40,
  cluster_rows = TRUE,     # Clusterizar las filas (eventos de splicing)
  cluster_cols = TRUE,     # Clusterizar las columnas (condiciones MDA A y MDA G)
  scale = "none",          # No escalar los datos
  color = colorRampPalette(c("skyblue", "white", "darkblue"))(100),  # Nueva paleta de colores
  main = "Heatmap de Eventos de Splicing Más Significativos (PSI MDA A vs MDA G)",
  labels_row = sig_splicing_sorted$geneSymbol[1:40]  # Usar solo geneSymbol como etiquetas de filas
)



# Combinar el número total de eventos y el número de eventos diferenciales
combined_splicing <- left_join(events_per_gene, differential_events, by = "geneSymbol")

# Combinar con los eventos sobreexpresados
combined_splicing <- left_join(combined_splicing, overexpressed_events, by = "geneSymbol")

# Combinar con los eventos subexpresados
combined_splicing <- left_join(combined_splicing, underexpressed_events, by = "geneSymbol")

# Rellenar los NA con 0 para las columnas que faltan valores
combined_splicing <- combined_splicing %>%
  mutate(
    n_differential_events = ifelse(is.na(n_differential_events), 0, n_differential_events),
    n_overexpressed_events = ifelse(is.na(n_overexpressed_events), 0, n_overexpressed_events),
    n_underexpressed_events = ifelse(is.na(n_underexpressed_events), 0, n_underexpressed_events)
  )

# Visualizar el resultado final
head(combined_splicing)

# Sumar los valores de cada columna
total_n_events <- sum(combined_splicing$n_events, na.rm = TRUE)
total_n_differential_events <- sum(combined_splicing$n_differential_events, na.rm = TRUE)
total_n_overexpressed_events <- sum(combined_splicing$n_overexpressed_events, na.rm = TRUE)
total_n_underexpressed_events <- sum(combined_splicing$n_underexpressed_events, na.rm = TRUE)

# Imprimir los resultados
cat("Suma total de n_events:", total_n_events, "\n")
cat("Suma total de n_differential_events:", total_n_differential_events, "\n")
cat("Suma total de n_overexpressed_events:", total_n_overexpressed_events, "\n")
cat("Suma total de n_underexpressed_events:", total_n_underexpressed_events, "\n")





taba_maestra <- read.csv("/home/administrador/Escritorio/rnaseq/DEA_DTU_DEU_combined_data.csv")



#Agregar las columnas de la tabla combined_splicing a taba_maestra
taba_maestra_extended <- left_join(taba_maestra, combined_splicing, by = c("gene_symbol" = "geneSymbol"))
# Reemplazar los valores NA por 0 en las columnas de eventos de splicing
taba_maestra_extended <- taba_maestra_extended %>%
  mutate(
    n_events = ifelse(is.na(n_events), 0, n_events),
    n_differential_events = ifelse(is.na(n_differential_events), 0, n_differential_events),
    n_overexpressed_events = ifelse(is.na(n_overexpressed_events), 0, n_overexpressed_events),
    n_underexpressed_events = ifelse(is.na(n_underexpressed_events), 0, n_underexpressed_events)
  )
# Visualizar las primeras filas de la tabla extendida
head(taba_maestra_extended)





# Sumar los valores de cada columna que comienza con "n_"
total_n_transcripts <- sum(taba_maestra_extended$n_transcripts, na.rm = TRUE)
total_n_differential_transcripts <- sum(taba_maestra_extended$n_differential_transcripts, na.rm = TRUE)
total_n_overexpressed_transcripts <- sum(taba_maestra_extended$n_overexpressed_transcripts, na.rm = TRUE)
total_n_underexpressed_transcripts <- sum(taba_maestra_extended$n_underexpressed_transcripts, na.rm = TRUE)
total_n_exons <- sum(taba_maestra_extended$n_exons, na.rm = TRUE)
total_n_differential_exons <- sum(taba_maestra_extended$n_differential_exons, na.rm = TRUE)
total_n_overexpressed_exons <- sum(taba_maestra_extended$n_overexpressed_exons, na.rm = TRUE)
total_n_underexpressed_exons <- sum(taba_maestra_extended$n_underexpressed_exons, na.rm = TRUE)
total_n_events <- sum(taba_maestra_extended$n_events, na.rm = TRUE)
total_n_differential_events <- sum(taba_maestra_extended$n_differential_events, na.rm = TRUE)
total_n_overexpressed_events <- sum(taba_maestra_extended$n_overexpressed_events, na.rm = TRUE)
total_n_underexpressed_events <- sum(taba_maestra_extended$n_underexpressed_events, na.rm = TRUE)

# Imprimir los resultados
cat("Suma total de n_transcripts:", total_n_transcripts, "\n")
cat("Suma total de n_differential_transcripts:", total_n_differential_transcripts, "\n")
cat("Suma total de n_overexpressed_transcripts:", total_n_overexpressed_transcripts, "\n")
cat("Suma total de n_underexpressed_transcripts:", total_n_underexpressed_transcripts, "\n")
cat("Suma total de n_exons:", total_n_exons, "\n")
cat("Suma total de n_differential_exons:", total_n_differential_exons, "\n")
cat("Suma total de n_overexpressed_exons:", total_n_overexpressed_exons, "\n")
cat("Suma total de n_underexpressed_exons:", total_n_underexpressed_exons, "\n")
cat("Suma total de n_events:", total_n_events, "\n")
cat("Suma total de n_differential_events:", total_n_differential_events, "\n")
cat("Suma total de n_overexpressed_events:", total_n_overexpressed_events, "\n")
cat("Suma total de n_underexpressed_events:", total_n_underexpressed_events, "\n")


write.csv(taba_maestra_extended, "/home/administrador/Escritorio/rnaseq/DEA_DTU_DEU_Splice.csv")










# Cargar los datos
tabla_maestra <- read.csv("/home/administrador/Escritorio/rnaseq/Tabla_maestra.csv")

# Eliminar la primera columna si es necesario
tabla_maestra <- tabla_maestra[, -1]

# Renombrar las columnas
colnames(tabla_maestra)[which(names(tabla_maestra) == "n_differential_transcripts")] <- "DTU"
colnames(tabla_maestra)[which(names(tabla_maestra) == "n_differential_exons")] <- "DEU"
colnames(tabla_maestra)[which(names(tabla_maestra) == "n_differential_events")] <- "Eventos_Splicing"
colnames(tabla_maestra)[which(names(tabla_maestra) == "Ediciones")] <- "Ediciones"

# Crear la columna para DEA (sin intersección)
tabla_maestra$DEA <- ifelse(
  tabla_maestra$DTU == 0 &
    tabla_maestra$DEU == 0 &
    tabla_maestra$Eventos_Splicing == 0 &
    tabla_maestra$Ediciones == 0,
  1, 0
)

# Verificar los nombres de las columnas
print(colnames(tabla_maestra))

# Seleccionar las columnas necesarias para el UpSet
upset_data <- tabla_maestra[, c("DTU", "DEU", "Eventos_Splicing", "Ediciones", "DEA")]

# Convertir a binario para el UpSet plot
upset_data_bin <- as.data.frame(lapply(upset_data, function(x) as.integer(x > 0)))

# Generar el UpSet plot
upset(
  upset_data_bin,
  sets = c("DTU", "DEU", "Eventos_Splicing", "Ediciones", "DEA"),
  keep.order = TRUE,
  order.by = "freq",
  sets.bar.color = "#56B4E9",
  main.bar.color = "#7fbf7b",
  matrix.color = "#af8dc3",
  text.scale = 2
)



# Identificar los genes que están en todas las categorías excepto DEA
interseccion_menos_dea <- tabla_maestra[
  upset_data_bin$DTU == 1 & 
    upset_data_bin$DEU == 1 & 
    upset_data_bin$Eventos_Splicing == 1 & 
    upset_data_bin$Ediciones == 1 & 
    upset_data_bin$DEA == 0, 
  "gene_symbol"
]

# Mostrar los nombres de los genes en esta intersección
print(interseccion_menos_dea)




# Cargar los datos
tabla_maestra <- read.csv("/home/administrador/Escritorio/rnaseq/Tabla_maestra.csv")

# Eliminar la primera columna si es necesario
tabla_maestra <- tabla_maestra[, -1]

# Renombrar las columnas
colnames(tabla_maestra)[which(names(tabla_maestra) == "n_differential_transcripts")] <- "DTU"
colnames(tabla_maestra)[which(names(tabla_maestra) == "n_differential_exons")] <- "DEU"
colnames(tabla_maestra)[which(names(tabla_maestra) == "n_differential_events")] <- "Eventos_Splicing"
colnames(tabla_maestra)[which(names(tabla_maestra) == "Ediciones")] <- "Ediciones"

# Seleccionar las columnas necesarias para el UpSet (sin DEA)
upset_data <- tabla_maestra[, c("DTU", "DEU", "Eventos_Splicing", "Ediciones")]

# Convertir a binario para el UpSet plot
upset_data_bin <- as.data.frame(lapply(upset_data, function(x) as.integer(x > 0)))

# Generar el UpSet plot en blanco y negro
upset(
  upset_data_bin,
  sets = c("DTU", "DEU", "Eventos_Splicing", "Ediciones"),
  keep.order = TRUE,
  order.by = "freq",
  sets.bar.color = "black",
  main.bar.color = "black",
  matrix.color = "black",
  text.scale = 2
)

# Identificar los genes que están en todas las categorías
interseccion_genes <- tabla_maestra[
  upset_data_bin$DTU == 1 & 
    upset_data_bin$DEU == 1 & 
    upset_data_bin$Eventos_Splicing == 1 & 
    upset_data_bin$Ediciones == 1, 
  "gene_symbol"
]

# Mostrar los nombres de los genes en esta intersección
print(interseccion_genes)
