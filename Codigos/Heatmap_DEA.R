# Cargar las librerías necesarias
library(DESeq2)
library(tidyverse) # para manipulación de datos
library(sva)       # para corrección de efectos de batch
library(dplyr)
library(org.Hs.eg.db)
library(ggpubr)

# Leer y preparar los datos
countData <- read.table("/home/administrador/Escritorio/tesis_MDA/RNASeq/counts.csv", header = TRUE, row.names = 1, sep = ",")
# Extraer la primera columna y agregar los rownames como una nueva columna
Name_genes <- countData[, c("gene_name")]
Name_genes <- cbind(ensembl_id = rownames(countData), Name_genes)
# Eliminar la columna gene_name
countData <- countData[, -1]

# Crear el data frame colData con la información de las muestras y condiciones
sample_names <- colnames(countData)
countData[sample_names] <- sapply(countData[sample_names], function(x) round(as.numeric(x)))

colData <- data.frame(
  sampleID = sample_names,
  condition = rep(c("G", "A"), each = 3)
)

# Crear el objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData[, sample_names],
                              colData = colData,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, "G")
#A Tratadas G control

# Realizar el análisis DESeq2
dds <- DESeq(dds)
# Suponiendo que ya tienes el objeto `dds` de DESeq2

# Realizar la transformación rlog
rlog_data <- rlog(dds, blind = FALSE)

# Obtener y ajustar los resultados usando la función results de DESeq2 explícitamente
results_deseq2 <- DESeq2::results(dds)
results_deseq2 <- DESeq2::results(dds, alpha = 0.05, pAdjustMethod = "BH")

# Convertir el objeto DESeqResults a un data.frame
results_df <- as.data.frame(results_deseq2)
# Agregar una columna con los nombres de los genes
results_df$ensembl_id <- rownames(results_df)
results_df <- merge(results_df, Name_genes, by = "ensembl_id")


genes_dea <- results_df$Name_genes  # Suponiendo que la columna de símbolos en sig_dea es Name_genes

# Conecta a la base de datos de Ensembl

# Obtener la información sobre los biotipos
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = genes_dea, 
                    columns = c("SYMBOL", "GENENAME", "GENETYPE"),
                    keytype = "SYMBOL")

# Filtrar y excluir los pseudogenes
non_pseudogenes <- gene_info[gene_info$GENETYPE != "pseudogene" & gene_info$GENETYPE != "pseudo" & !is.na(gene_info$GENETYPE), ]

# Mostrar la data sin pseudogenes
print(non_pseudogenes)

# Crear un vector con los nombres de los pseudogenes
pseudogene_names <- non_pseudogenes$SYMBOL[non_pseudogenes$GENETYPE != "pseudogene" & non_pseudogenes$GENETYPE != "pseudo"]



# Filtrar sig_dea para que solo incluya los genes que no son pseudogenes
filtered_dea <- results_df[results_df$Name_genes %in% pseudogene_names, ]

# Mostrar el resultado filtrado
head(filtered_dea)


significant_genes_clean <- na.omit(filtered_dea)
significant_genes_clean <- significant_genes_clean[significant_genes_clean$padj <= 0.05, ]


# Filtrar genes subexpresados (log2FoldChange <= -1)
subexpresados <- significant_genes_clean %>%
  filter(log2FoldChange <= -1)

# Filtrar genes sobreexpresados (log2FoldChange >= 1)
sobreexpresados <- significant_genes_clean %>%
  filter(log2FoldChange >= 1)


# Obtener la matriz transformada con rlog
rlog_matrix <- assay(rlog_data)

# Filtrar los genes rlog para los sobreexpresados y subexpresados
rlog_sobreexpresados <- rlog_matrix[rownames(rlog_matrix) %in% sobreexpresados$ensembl_id, ]
rlog_subexpresados <- rlog_matrix[rownames(rlog_matrix) %in% subexpresados$ensembl_id, ]
rlog_sobreexpresados
# Calcular la varianza para cada gen en los conjuntos sobreexpresados y subexpresados
var_sobreexpresados <- apply(rlog_sobreexpresados, 1, var)
var_subexpresados <- apply(rlog_subexpresados, 1, var)


# Seleccionar los 10 genes más variables de sobreexpresión
top10_sobreexpresados <- names(sort(var_sobreexpresados, decreasing = TRUE)[1:17])

# Seleccionar los 10 genes más variables de subexpresión
top10_subexpresados <- names(sort(var_subexpresados, decreasing = TRUE)[1:30])

# Confirmar la extracción de los genes seleccionados
cat("Genes seleccionados para sobreexpresión:\n")
print(top10_sobreexpresados)

cat("\nGenes seleccionados para subexpresión:\n")
print(top10_subexpresados)

# Combinar los datos de los genes seleccionados
genes_seleccionados <- c(top10_sobreexpresados, top10_subexpresados)
heatmap_data <- rlog_matrix[rownames(rlog_matrix) %in% genes_seleccionados, ]

# Verificar que los genes seleccionados están en la matriz de datos
cat("\nGenes presentes en heatmap_data:\n")
print(rownames(heatmap_data))

# Asegurarse de que los nombres de filas en heatmap_data coincidan con los nombres de genes en lugar de los IDs de Ensembl
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]  # Eliminar filas con NA
rownames(heatmap_data) <- results_df$Name_genes[match(rownames(heatmap_data), results_df$ensembl_id)]

# Verificar los nombres de las filas después de la conversión
cat("\nNombres de las filas en heatmap_data después de la conversión:\n")
print(rownames(heatmap_data))

# Crear anotaciones para diferenciar sobreexpresados y subexpresados, confirmando que coinciden con los genes seleccionados
annotation_row <- data.frame(
  `Gene Type` = factor(ifelse(rownames(heatmap_data) %in% results_df$Name_genes[results_df$ensembl_id %in% top10_sobreexpresados], 
                              "Overexpressed", 
                              "Underexpressed"))
)
rownames(annotation_row) <- rownames(heatmap_data)

# Verificar las anotaciones
cat("\nAnotaciones para cada gen:\n")
print(annotation_row)


# Crear anotaciones para las condiciones de las muestras
annotation_col <- data.frame(
  `Condition` = colData$condition,
  row.names = colnames(heatmap_data)
)

# Definir 11 puntos de ruptura (breaks) para que coincidan con los 11 colores
breaks <- seq(min(heatmap_data), max(heatmap_data), length.out = 11)

# Definir la paleta de colores usando la paleta personalizada con los puntos de ruptura
heatmap_colors <- colorRamp2(breaks,
                             c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8",
                               "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", 
                               "#00441b"))

# Verificar la cantidad de grupos en `row_split`
num_slices <- length(unique(annotation_row$`Gene Type`))


annotation_row$Gene.Type <- as.factor(annotation_row$Gene.Type)

# Cambiar los colores para las barras de "Gene Type"
row_ha <- rowAnnotation(
  "Gene Type" = annotation_row$Gene.Type,
  col = list(Gene.Type = c("Overexpressed" = "#40004b", "Underexpressed" = "#1b7837")),  # Aquí cambias los colores
  annotation_legend_param = list(title = "Gene Type")
)

head(row_ha)
ht <- Heatmap(heatmap_data,
              name = "Expression",
              col = heatmap_colors,
              top_annotation = HeatmapAnnotation(df = annotation_col,
                                                 col = list(Condition = c("G" = "blue", "A" = "red"))),
              left_annotation = row_ha,  # Añadir la anotación lateral
              row_split = annotation_row$Gene.Type,  # Separar filas por tipo de gen (clustering en dos grupos)
              cluster_rows = TRUE,  # Clusterizar filas dentro de cada grupo
              cluster_columns = TRUE,  # Clusterizar columnas
              row_title = c("Overexpressed Genes", "Underexpressed Genes"),  # Títulos para cada grupo
              column_names_gp = gpar(fontsize = 6),  # Tamaño de la fuente para los nombres de los genes
              row_names_gp = gpar(fontsize = 8),  # Tamaño de la fuente para los nombres de las muestras
              border = TRUE,  # Añadir borde alrededor de las celdas
              heatmap_legend_param = list(title = "Expression (rlog values)"),
              column_title = "Heatmap of Gene Expression Levels (rlog)",  # Aquí agregas el título
              column_title_gp = gpar(fontsize = 16, fontface = "bold")  # Ajustes de estilo para el título
)
ht
# Ajustar el tamaño de la figura y guardar como PDF o PNG
png("clustered_heatmap_with_column_title.png", width = 1200, height = 1600, res = 150)
draw(ht)
dev.off()



# Cambiar los colores para las barras de "Gene Type"
row_ha <- rowAnnotation(
  "Gene Type" = annotation_row$Gene.Type,
  col = list(Gene.Type = c("Overexpressed" = "#40004b", "Underexpressed" = "#1b7837")),
  annotation_legend_param = list(title = "Gene Type")
)

# Crear el heatmap
ht <- Heatmap(heatmap_data,
              name = "Expression",
              col = heatmap_colors,
              top_annotation = HeatmapAnnotation(df = annotation_col,
                                                 col = list(Condition = c("G" = "blue", "A" = "red"))),
              left_annotation = row_ha,
              row_split = annotation_row$Gene.Type,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              row_title = c("Overexpressed Genes", "Underexpressed Genes"),
              column_names_gp = gpar(fontsize = 6),  # Incrementar tamaño de la fuente para columnas
              row_names_gp = gpar(fontsize = 8),  # Incrementar tamaño de la fuente para filas
              border = TRUE,
              heatmap_legend_param = list(title = "Expression (rlog values)"),
              column_title = "Heatmap of Gene Expression Levels (rlog)",
              column_title_gp = gpar(fontsize = 18, fontface = "bold")  # Ajustar el tamaño del título
)

# Ajustar el tamaño de la figura y guardar como PDF o PNG
png("larger_clustered_heatmap.png", width = 1800, height = 2400, res = 200)
draw(ht)
dev.off()
# Ajustar el tamaño de la figura y guardar como PDF o PNG
png("extra_large_clustered_heatmap.png", width = 1800, height = 3000, res = 200)
draw(ht)
dev.off()


# Filtrar los datos para los genes sobreexpresados
heatmap_sobreexpresados <- rlog_sobreexpresados[rownames(rlog_sobreexpresados) %in% top10_sobreexpresados, ]

# Asegurarse de que los nombres de filas en heatmap_sobreexpresados coincidan con los nombres de genes en lugar de los IDs de Ensembl
rownames(heatmap_sobreexpresados) <- results_df$Name_genes[match(rownames(heatmap_sobreexpresados), results_df$ensembl_id)]

# Crear el heatmap para los genes sobreexpresados
ht_sobreexpresados <- Heatmap(heatmap_sobreexpresados,
                              name = "Expression",
                              col = heatmap_colors,
                              top_annotation = HeatmapAnnotation(df = annotation_col,
                                                                 col = list(Condition = c("G" = "blue", "A" = "red"))),
                              cluster_rows = TRUE,  # Clusterizar filas
                              cluster_columns = TRUE,  # Clusterizar columnas
                              row_names_gp = gpar(fontsize = 8),  # Tamaño de la fuente para los nombres de las muestras
                              border = TRUE,  # Añadir borde alrededor de las celdas
                              heatmap_legend_param = list(title = "Expression (rlog values)"),
                              column_title = "Heatmap of Overexpressed Genes (rlog)",  # Título del heatmap
                              column_title_gp = gpar(fontsize = 16, fontface = "bold")  # Estilo del título
)

# Guardar el heatmap para los genes sobreexpresados
png("heatmap_sobreexpresados.png", width = 1200, height = 1600, res = 150)
draw(ht_sobreexpresados)
dev.off()



# Filtrar los datos para los genes subexpresados
heatmap_subexpresados <- rlog_subexpresados[rownames(rlog_subexpresados) %in% top10_subexpresados, ]

# Asegurarse de que los nombres de filas en heatmap_subexpresados coincidan con los nombres de genes en lugar de los IDs de Ensembl
rownames(heatmap_subexpresados) <- results_df$Name_genes[match(rownames(heatmap_subexpresados), results_df$ensembl_id)]

# Crear el heatmap para los genes subexpresados
ht_subexpresados <- Heatmap(heatmap_subexpresados,
                            name = "Expression",
                            col = heatmap_colors,
                            top_annotation = HeatmapAnnotation(df = annotation_col,
                                                               col = list(Condition = c("G" = "blue", "A" = "red"))),
                            cluster_rows = TRUE,  # Clusterizar filas
                            cluster_columns = TRUE,  # Clusterizar columnas
                            row_names_gp = gpar(fontsize = 8),  # Tamaño de la fuente para los nombres de las muestras
                            border = TRUE,  # Añadir borde alrededor de las celdas
                            heatmap_legend_param = list(title = "Expression (rlog values)"),
                            column_title = "Heatmap of Underexpressed Genes (rlog)",  # Título del heatmap
                            column_title_gp = gpar(fontsize = 16, fontface = "bold")  # Estilo del título
)

# Guardar el heatmap para los genes subexpresados
png("heatmap_subexpresados.png", width = 1200, height = 1600, res = 150)
draw(ht_subexpresados)
dev.off()












library(ComplexHeatmap)
library(circlize)

# Seleccionar los 10 genes más variables de sobreexpresión
top10_sobreexpresados <- names(sort(var_sobreexpresados, decreasing = TRUE)[1:17])

# Seleccionar los 10 genes más variables de subexpresión
top10_subexpresados <- names(sort(var_subexpresados, decreasing = TRUE)[1:60])

# Combinar los datos de los genes seleccionados
genes_seleccionados <- c(top10_sobreexpresados, top10_subexpresados)
heatmap_data <- rlog_matrix[rownames(rlog_matrix) %in% genes_seleccionados, ]
heatmap_data
# Asegurarse de que los nombres de filas en heatmap_data coincidan con los nombres de genes en lugar de los IDs de Ensembl
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]  # Eliminar filas con NA
rownames(heatmap_data) <- results_df$Name_genes[match(rownames(heatmap_data), results_df$ensembl_id)]

heatmap_data# Crear anotaciones para diferenciar sobreexpresados y subexpresados

annotation_row <- data.frame(
  `Gene Type` = factor(c(rep("Overexpressed", length(top10_sobreexpresados)), 
                         rep("Underexpressed", length(top10_subexpresados))))
)
rownames(annotation_row) <- rownames(heatmap_data)
annotation_row




# Seleccionar los 10 genes más variables de sobreexpresión
top10_sobreexpresados <- names(sort(var_sobreexpresados, decreasing = TRUE)[1:10])

# Seleccionar los 10 genes más variables de subexpresión
top10_subexpresados <- names(sort(var_subexpresados, decreasing = TRUE)[1:60])

# Confirmar la extracción de los genes seleccionados
cat("Genes seleccionados para sobreexpresión:\n")
print(top10_sobreexpresados)

cat("\nGenes seleccionados para subexpresión:\n")
print(top10_subexpresados)

# Combinar los datos de los genes seleccionados
genes_seleccionados <- c(top10_sobreexpresados, top10_subexpresados)
heatmap_data <- rlog_matrix[rownames(rlog_matrix) %in% genes_seleccionados, ]

# Verificar que los genes seleccionados están en la matriz de datos
cat("\nGenes presentes en heatmap_data:\n")
print(rownames(heatmap_data))

# Asegurarse de que los nombres de filas en heatmap_data coincidan con los nombres de genes en lugar de los IDs de Ensembl
heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]  # Eliminar filas con NA
rownames(heatmap_data) <- results_df$Name_genes[match(rownames(heatmap_data), results_df$ensembl_id)]

# Verificar los nombres de las filas después de la conversión
cat("\nNombres de las filas en heatmap_data después de la conversión:\n")
print(rownames(heatmap_data))

# Crear anotaciones para diferenciar sobreexpresados y subexpresados, confirmando que coinciden con los genes seleccionados
annotation_row <- data.frame(
  `Gene Type` = factor(ifelse(rownames(heatmap_data) %in% results_df$Name_genes[results_df$ensembl_id %in% top10_sobreexpresados], 
                              "Overexpressed", 
                              "Underexpressed"))
)
rownames(annotation_row) <- rownames(heatmap_data)

# Verificar las anotaciones
cat("\nAnotaciones para cada gen:\n")
print(annotation_row)












# Crear anotaciones para las condiciones de las muestras
annotation_col <- data.frame(
  `Condition` = colData$condition,
  row.names = colnames(heatmap_data)
)

# Definir 11 puntos de ruptura (breaks) para que coincidan con los 11 colores
breaks <- seq(min(heatmap_data), max(heatmap_data), length.out = 11)

# Definir la paleta de colores usando la paleta personalizada con los puntos de ruptura
heatmap_colors <- colorRamp2(breaks,
                             c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8",
                               "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", 
                               "#00441b"))

# Cambiar los colores para las barras de "Gene Type"
row_ha <- rowAnnotation(
  "Gene Type" = annotation_row$Gene.Type,
  col = list(Gene.Type = c("Overexpressed" = "#1b7837", "Underexpressed" = "#40004b")),  # Cambiar colores
  annotation_legend_param = list(title = "Gene Type")
)

# Crear el heatmap
ht <- Heatmap(heatmap_data,
              name = "Expression",
              col = heatmap_colors,
              top_annotation = HeatmapAnnotation(df = annotation_col,
                                                 col = list(Condition = c("G" = "blue", "A" = "red"))),
              left_annotation = row_ha,  # Añadir la anotación lateral
              row_split = annotation_row$Gene.Type,  # Separar filas por tipo de gen (clustering en dos grupos)
              cluster_rows = TRUE,  # Clusterizar filas dentro de cada grupo
              cluster_columns = TRUE,  # Clusterizar columnas
              row_title = c( "Overexpressed Genes", "Underexpressed Genes"),  # Títulos para cada grupo
              column_names_gp = gpar(fontsize = 6),  # Tamaño de la fuente para los nombres de los genes
              row_names_gp = gpar(fontsize = 8),  # Tamaño de la fuente para los nombres de las muestras
              border = TRUE,  # Añadir borde alrededor de las celdas
              heatmap_legend_param = list(title = "Expression (rlog values)"),
              column_title = "Heatmap of Gene Expression Levels (rlog)",  # Aquí agregas el título
              column_title_gp = gpar(fontsize = 16, fontface = "bold")  # Ajustes de estilo para el título
)
ht

# Ajustar el tamaño de la figura y guardar como PDF o PNG
png("sub_over_expressed_genes_heatmap.png", width = 1200, height = 1600, res = 150)
draw(ht)
dev.off()



