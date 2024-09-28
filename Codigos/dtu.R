library(edgeR)
library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)
# Definir las rutas a los archivos
counts_file <- "/home/administrador/Escritorio/rnaseq/salmon.merged.transcript_counts.tsv"
head()
# Leer los archivos
counts <- read_delim(counts_file, delim = "\t")
tx2gene <- read_delim(tx2gene_file, delim = "\t", col_names = c("transcript_id", "gene_id", "gene_symbol"))
# Extraer la matriz de datos de conteos, excluyendo las columnas de identificadores
counts_matrix <- as.matrix(counts[, 3:ncol(counts)])

# Asignar nombres de filas a counts_matrix usando 'transcript_id'
rownames(counts_matrix) <- counts$tx

# Añadir una columna extra con gene_symbol
counts_matrix_with_genes <- cbind(GeneSymbol = tx2gene$gene_symbol, counts_matrix)

# Asegurarse de que todos los datos en counts_matrix_with_genes[, -1] son numéricos
# Aquí conservamos los rownames
numeric_data <- counts_matrix  # Simplemente usamos counts_matrix ya que es numérica y tiene rownames

# Verificar que los rownames se conservan
print(head(numeric_data))

# Crear el objeto DGEList con los datos numéricos de conteo
dge <- DGEList(counts = numeric_data)
# Normalización de los factores
dge <- calcNormFactors(dge)
# Definir las condiciones experimentales, con MDA_G como el primer nivel (control)
group <- factor(c("MDA_A", "MDA_A", "MDA_A", "MDA_G", "MDA_G", "MDA_G"), levels = c("MDA_G", "MDA_A"))
# Asignar el grupo al objeto DGEList
dge$samples$group <- group
# Estimación de la dispersión
dge <- estimateDisp(dge)
# Ajuste del modelo lineal generalizado (GLM)
fit <- glmFit(dge, design = model.matrix(~group))
# Realizar el test de razón de verosimilitud
lrt <- glmLRT(fit)
dtu_results <- topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "PValue")
sig_dtu <- dtu_results$table[dtu_results$table$FDR < 0.05, ]

# Agregar transcript_id como rownames
sig_dtu$transcript_id <- rownames(sig_dtu)

# Agregar gene_symbol basado en transcript_id
sig_dtu$gene_symbol <- tx2gene$gene_symbol[match(sig_dtu$transcript_id, tx2gene$transcript_id)]

# Verificar los resultados
head(sig_dtu)
head(dtu_results)
write.csv(sig_dtu, "/home/administrador/Escritorio/rnaseq/DTU_SIG.csv")
write.csv(dtu_results, "/home/administrador/Escritorio/rnaseq/DTU.csv")

# Agregar gene_symbol basado en transcript_id
sig_dtu$gene_symbol <- tx2gene$gene_symbol[match(sig_dtu$transcript_id, tx2gene$transcript_id)]

genes_of_interest <- c("U2", "NPM1P7", "RPL32P29", "HCN2", "FAM86MP", "LncRNA", "ADAR")
filtered_sig_dtu <- sig_dtu %>% filter(gene_symbol %in% genes_of_interest)

head(dtu_results)
sig_dtu <- dtu_results$table[dtu_results$table$FDR < 0.05, ]

# Crear una tabla resumen
# Contar el número total de transcritos antes del filtrado
total_transcripts_before <- nrow(dtu_results$table)

# Mostrar el número total de transcritos antes del filtrado
total_transcripts_before
# Contar el número de transcritos después del filtrado
total_transcripts_after <- nrow(sig_dtu)

# Mostrar el número de transcritos después del filtrado
total_transcripts_after

# Crear una tabla resumen
filter_summary <- data.frame(
  Stage = c("Antes", "Despues"),
  Number_of_Transcripts = c(total_transcripts_before, total_transcripts_after)
)

# Mostrar la tabla resumen
filter_summary


head(sig_dtu)
library(dplyr)

# Añadir transcript_id como columna en sig_dtu si aún no lo has hecho
sig_dtu <- sig_dtu %>%
  mutate(transcript_id = rownames(sig_dtu))

# Unir las tablas por transcript_id
sig_dtu_annotated <- sig_dtu %>%
  left_join(tx2gene, by = "transcript_id")

# Ver los primeros registros para verificar la unión
head(sig_dtu_annotated)


# Agrupar por gene_symbol y sumarizar los datos
gene_summary <- sig_dtu_annotated %>%
  group_by(gene_symbol) %>%
  summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    min_pvalue = min(PValue, na.rm = TRUE),
    min_fdr = min(FDR, na.rm = TRUE),
    n_transcripts = n()
  )

library(ggplot2)

# Filtrar los top 15 genes por número de transcritos
top_genes <- gene_summary %>%
  top_n(15, n_transcripts)

# Crear un barplot para los top 15 genes
ggplot(top_genes, aes(x = reorder(gene_symbol, -n_transcripts), y = n_transcripts, fill = gene_symbol)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Para girar el gráfico y facilitar la lectura
  theme_minimal() +
  labs(title = "Top 15 Genes por Número de Transcritos",
       x = "Gen",
       y = "Número de Transcritos",
       fill = "Gen") +
  theme(legend.position = "none")  # Quitar la leyenda si no es necesaria


library(ggplot2)

library(ggplot2)
library(dplyr)

# Filtrar los genes más sobreexpresados y subexpresados
top_sobreexpresados <- sig_dtu_annotated %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 15) %>%
  mutate(category = "Sobreexpresados")

top_subexpresados <- sig_dtu_annotated %>%
  arrange(logFC) %>%
  slice_head(n = 15) %>%
  mutate(category = "Subexpresados")

# Combinar los datos en un solo dataframe
combined_data <- bind_rows(top_subexpresados, top_sobreexpresados)

# Crear el gráfico con facetas
ggplot(combined_data, aes(x = reorder(gene_symbol, logFC), y = logFC, fill = logFC)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Top 15 Genes por logFC",
       x = "Gene",
       y = "logFC",
       fill = "logFC") +
  facet_wrap(~category, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10))






# Filtrar los top 15 genes por logFC promedio
top_genes <- gene_summary %>%
  top_n(abs(mean_logFC))

# Crear un scatter plot para los top 15 genes
ggplot(top_genes, aes(x = mean_logFC, y = -log10(min_pvalue), label = gene_symbol)) +
  geom_point(aes(size = n_transcripts, color = min_fdr)) +
  geom_text(vjust = 1.5, hjust = 1.5) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Top 15 Genes por logFC Promedio vs. P-Value Mínimo",
       x = "logFC Promedio",
       y = "-log10(P-Value Mínimo)",
       color = "FDR Mínimo",
       size = "Número de Transcritos")


library(ggplot2)

# Crear un scatter plot para los top 15 genes
ggplot(top_genes, aes(x = mean_logFC, y = -log10(min_pvalue), label = gene_symbol)) +
  geom_point(aes(size = n_transcripts, color = min_fdr)) +
  geom_text(vjust = 1.5, hjust = 1.5) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(min(top_genes$mean_logFC) - 1, max(top_genes$mean_logFC) + 1)) +
  scale_y_continuous(limits = c(0, max(-log10(top_genes$min_pvalue)) + 1)) +
  theme_minimal() +
  labs(title = "Top 15 Genes por logFC Promedio vs. P-Value Mínimo",
       x = "logFC Promedio",
       y = "-log10(P-Value Mínimo)",
       color = "FDR Mínimo",
       size = "Número de Transcritos")
# Ver la tabla resumen por gene_symbol
head(gene_summary)
# Filtrar los top 15 genes por logFC promedio
top_genes <- gene_summary %>%
  top_n(15, abs(mean_logFC))

library(tibble)
library(pheatmap)



# Crear un gráfico de barras para comparar el número de transcritos antes y después del filtrado
ggplot(filter_summary, aes(x = Stage, y = Number_of_Transcripts, fill = Stage)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Número de Transcritos Antes y Después del Filtrado",
       x = "Etapa",
       y = "Número de Transcritos") +
  theme(legend.position = "none")


summary_stats <- sig_dtu %>%
  summarise(
    mean_logFC = mean(logFC),
    median_logFC = median(logFC),
    min_logFC = min(logFC),
    max_logFC = max(logFC),
    mean_pvalue = mean(PValue),
    median_pvalue = median(PValue),
    mean_FDR = mean(FDR)
  )
summary_stats
ggplot(sig_dtu, aes(x = logFC, y = -log10(PValue), color = FDR < 0.05)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano Plot de Transcritos Filtrados",
       x = "Log2 Fold Change", y = "-Log10 P-Value") +
  scale_color_manual(values = c("grey", "red"))

library(tidyr)

# Preparar los datos
heatmap_data <- sig_dtu %>%
  dplyr::select(transcript_id, gene_symbol, logFC) %>%
  spread(key = transcript_id, value = logFC)

# Generar el heatmap
pheatmap::pheatmap(as.matrix(heatmap_data[-1]),
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   annotation_row = heatmap_data[1],
                   main = "Heatmap de Log2 Fold Change en Transcritos Filtrados")


# Filtrar por transcritos significativos
significant_transcripts <- sig_dtu %>%
  filter(FDR < 0.05) %>%
  group_by(gene_symbol) %>%
  summarise(
    num_transcripts = n(),
    avg_logFC = mean(logFC),
    avg_logCPM = mean(logCPM),
    min_pvalue = min(PValue),
    max_pvalue = max(PValue),
    min_FDR = min(FDR),
    max_FDR = max(FDR)
  ) %>%
  arrange(desc(num_transcripts))

significant_transcripts


# Lista de genes de interés
genes_of_interest <- c("U2", "NPM1P7", "RPL32P29", "HCN2", "FAM86MP", "LncRNA", "ADAR")

# Filtrar los resultados significativos para los genes de interés
specific_genes_results <- sig_dtu %>%
  filter(gene_symbol %in% genes_of_interest)

los genes de interés
specific_genes_results_sig <- significant_transcripts %>%
  filter(gene_symbol %in% genes_of_interest)
# Ver los resultados filtrados
print(specific_genes_results)

library(ggplot2)
library(reshape2)

# Supongamos que quieres ver las proporciones relativas de isoformas para un gen específico
gene_of_interest <- "ADAR"  # Cambia esto al gene_id de tu interés

# Filtrar los datos para obtener solo las isoformas del gen de interés
isoform_counts_gene <- counts_matrix_with_genes[tx2gene$gene_symbol == gene_of_interest, ]

# Asegurarse de que los datos sean numéricos
numeric_isoform_counts_gene <- apply(isoform_counts_gene[, -1], 2, as.numeric)
rownames(numeric_isoform_counts_gene) <- isoform_counts_gene[, 1]

# Calcular las proporciones de isoformas para cada condición
isoform_proportions_gene <- prop.table(numeric_isoform_counts_gene, margin = 2)

# Convertir a formato largo para ggplot2
proportions_df_gene <- melt(isoform_proportions_gene, varnames = c("feature_id", "condition"), value.name = "proportion")

# Agregar el gene_id a la tabla para que pueda ser filtrado en el gráfico
proportions_df_gene$gene_id <- gene_of_interest

# Graficar las proporciones relativas de las isoformas del gen de interés
ggplot(proportions_df_gene, aes(x = feature_id, y = proportion, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = paste("Proporciones Relativas de Transcritos para el Gen", gene_of_interest), 
       x = "Transcrito", y = "Proporción")





library(ggplot2)

# Definir umbrales para significancia
logFC_threshold <- 1  # Cambia según tu criterio
FDR_threshold <- 0.05  # FDR menor que 0.05 se considera significativo

# Crear una columna para resaltar los transcritos significativos
sig_dtu$significant <- ifelse(sig_dtu$FDR < FDR_threshold & abs(sig_dtu$logFC) > logFC_threshold, "Significativo", "No significativo")

# Resaltar específicamente los transcritos relacionados con ADAR
adar_related_genes <- c("ADAR", "NDRG2","ATP7A")  # Cambia esto con los genes de interés relacionados con ADAR
sig_dtu$highlight <- ifelse(sig_dtu$gene_symbol %in% adar_related_genes, "ADAR-Related", sig_dtu$significant)

# Crear el Gráfico de Volcán
ggplot(sig_dtu, aes(x = logFC, y = -log10(FDR), color = highlight)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  scale_color_manual(values = c("green", "grey", "red")) +
  labs(title = "Gráfico de Volcán para DTU",
       x = "Log2 Fold Change",
       y = "-Log10(FDR)",
       color = "Tipo de Transcrito") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = ifelse(highlight == "ADAR-Related", gene_symbol, "")),
                  size = 3, box.padding = 0.3)



library(pheatmap)

# Filtrar los transcritos significativos
significant_isoforms <- rownames(sig_dtu)
heatmap_matrix <- counts_matrix[significant_isoforms, ]

# Normalizar los datos (log-transformación opcional)
log_heatmap_matrix <- log2(heatmap_matrix + 1)

# Crear el heatmap
pheatmap(log_heatmap_matrix, 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = TRUE, 
         main = "Heatmap de Transcritos Significativos",
         annotation_row = data.frame(Gene = sig_dtu$gene_symbol[significant_isoforms]))




library(ggplot2)

# Definir los genes de interés
genes_interes <- c("ADAR", "NDRG2","ATP7A")  

isoform_counts <- counts_matrix_with_genes[tx2gene$gene_symbol %in% genes_interes, ]

# Asegurarse de que los nombres de columnas y filas sean correctos
rownames(isoform_counts) <- isoform_counts$GeneSymbol  # Asignar GeneSymbol como rownames
isoform_counts <- isoform_counts[, -1]  # Remover la columna GeneSymbol para que no interfiera en melt

# Convertir a formato largo
isoform_counts_long <- melt(isoform_counts, variable.name = "Condition", value.name = "Expression")

# Agregar la columna GeneSymbol nuevamente
isoform_counts_long$GeneSymbol <- rownames(isoform_counts_long)

# Graficar el perfil de expresión de los transcritos individuales
ggplot(isoform_counts_long, aes(x = Condition, y = Expression, color = GeneSymbol, group = GeneSymbol)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Perfil de Expresión de Transcritos Individuales",
       x = "Condición", y = "Expresión Relativa")





library(igraph)
library(ggraph)

# Crear un dataframe que conecte genes con transcritos
gene_transcript_edges <- data.frame(
  from = tx2gene$gene_symbol[match(rownames(counts_matrix), tx2gene$transcript_id)],
  to = rownames(counts_matrix)
)

# Crear el grafo
g <- graph_from_data_frame(gene_transcript_edges)

# Visualizar el grafo
ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5)) +
  geom_node_point(color = "blue", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle("Mapa de Conexión entre Genes y Transcritos")
 
