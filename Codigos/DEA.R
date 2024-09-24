# Cargar las librerías necesarias
library(DESeq2)
library(sva)     
library(dplyr)
library(org.Hs.eg.db)
library(ggpubr)
library(tidyverse) # para manipulación de datos


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

head(results_df)
  ################## PCA ##############################
# Transformar los datos de conteo para el análisis PCA
vsd <- vst(dds, blind = FALSE)


# Cargar las librerías necesarias
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)

vsd <- vst(dds, blind = FALSE)

# Realizar el PCA
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colores_condiciones <- c("G" = "#af8dc3", "A" = "#7fbf7b")

ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  xlim(-15, 15) +  # Adjust x-axis
  ylim(-15, 15) +  # Adjust y-axis
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # Vertical line at x=0
  scale_color_manual(values = colores_condiciones) +  # Apply specific colors for G and A conditions
  ggtitle("PCA Plot") +  # Add title in English
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.border = element_rect(colour = "black", fill=NA, size=1),  # Add a border to the plot
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center and style the title
  )







############### Volcano ####################################
# Supongamos que tu dataframe se llama results_df y tiene las columnas "padj" y "logfc"
results_df <- read.csv("/home/administrador/Escritorio/rnaseq/DEA_no_pseudo.csv")
significant_genes <- results_df[results_df$padj <= 0.05, ]
significant_genes_clean <- na.omit(significant_genes)
significant_genes <- significant_genes[significant_genes$padj != 0, ]

results_df <- na.omit(results_df)
write.csv(results_df, "/home/administrador/Escritorio/rnaseq/DEA.csv")

# Gen específico para destacar
specific_gene <- results_df[results_df$ensembl_id == "ENSG00000160710.18", ]
specific_gene 


# Filtrar genes subexpresados (log2FoldChange <= -1)
subexpresados <- significant_genes %>%
  filter(log2FoldChange < 0)

# Filtrar genes sobreexpresados (log2FoldChange >= 1)
sobreexpresados <- significant_genes %>%
  filter(log2FoldChange > 0)


# Contar las cantidades de genes significativos, subexpresados y sobreexpresados
num_significant_genes <- nrow(significant_genes_clean)
num_subexpresados <- sum(significant_genes_clean$log2FoldChange < 0)
num_sobreexpresados <- sum(significant_genes_clean$log2FoldChange > 0)

# Mostrar las cantidades
print(paste("Número de genes significativos:", num_significant_genes))
print(paste("Número de genes subexpresados:", num_subexpresados))
print(paste("Número de genes sobreexpresados:", num_sobreexpresados))

# Guardar los nombres de los genes subexpresados en un archivo
write.table(subexpresados$Name_genes, file = "subexpresados_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Guardar los nombres de los genes sobreexpresados en un archivo
write.table(sobreexpresados$Name_genes, file = "sobreexpresados_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Combinar ambos dataframes
combined_results <- bind_rows(subexpresados, sobreexpresados)

# Exportar el dataframe combinado a un archivo TSV
write.table(combined_results, 
            file = "/home/administrador/Escritorio/tesis_MDA/RNASeq/condition_control_overexpression.deseq2.results_filtered.tsv", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)



# Supongamos que 'subexpresados' y 'sobreexpresados' son dataframes con los datos relevantes.


library(pheatmap)
library(RColorBrewer)


# Obtener los 30 genes más sobreexpresados
top30_sobreexpresados <- sobreexpresados %>% arrange(desc(log2FoldChange)) %>% head(20)


# Obtener los 50 o 100 genes más subexpresados
top50_subexpresados <- subexpresados %>% arrange(log2FoldChange) %>% head(50)
# O para 100 genes subexpresados:

# Combinar los genes subexpresados y sobreexpresados
genes_interes <- rbind(top50_subexpresados, top30_sobreexpresados)




# Función para calcular genes significativos
calculate_significant_genes <- function(df) {
  df %>% filter(padj < 0.05) %>% nrow()
}

# Función para identificar los top 3 genes sobreexpresados y subexpresados
identify_top_genes <- function(results_df) {
  top_overexpressed <- results_df %>%
    filter(padj < 0.05, log2FoldChange > 1) %>%
    arrange(desc(log2FoldChange)) %>%
    slice_head(n = 3)
  
  top_underexpressed <- results_df %>%
    filter(padj < 0.05, log2FoldChange < -1) %>%
    arrange(log2FoldChange) %>%
    slice_head(n = 3)
  
  list(overexpressed = top_overexpressed, underexpressed = top_underexpressed)
}

# Función para crear el Volcano Plot con líneas de corte
create_volcano_plot <- function(results_df, top_overexpressed, top_underexpressed, specific_gene, num_significant_genes) {
  # Definir el color según la significancia y el cambio de pliegue
  results_df$color <- ifelse(
    results_df$log2FoldChange > 0 & results_df$padj < 0.05, "Significativo",
    ifelse(results_df$log2FoldChange < 0 & results_df$padj < 0.05, "Significativo", "No Significativo")
  )
  
  # Crear el gráfico Volcano Plot
  ggplot(results_df, aes(x = log2FoldChange, y = -log10(pvalue), color = color)) +
    geom_point(alpha = 0.7, size = 2) +  # Colorear puntos basados en la nueva columna 'color'
    scale_color_manual(values = c("Significativo" = "red", "No Significativo" = "grey")) +
    # Añadir los top sobreexpresados y subexpresados con colores específicos
    geom_point(data = top_overexpressed, 
               aes(x = log2FoldChange, y = -log10(pvalue)), 
               color = "#7fbf7b", size = 5) +  # Verde para sobreexpresados
    geom_point(data = top_underexpressed, 
               aes(x = log2FoldChange, y = -log10(pvalue)), 
               color = "#af8dc3", size = 5) +  # Morado para subexpresados
    geom_point(data = specific_gene, 
               aes(x = log2FoldChange, y = -log10(pvalue)), 
               color = "blue", size = 5) +  # Azul para genes específicos
    geom_text_repel(data = top_overexpressed, 
                    aes(x = log2FoldChange, y = -log10(pvalue), label = Name_genes), 
                    color = "#7fbf7b", size = 3, fontface = "bold") +
    geom_text_repel(data = top_underexpressed, 
                    aes(x = log2FoldChange, y = -log10(pvalue), label = Name_genes), 
                    color = "#af8dc3", size = 3, fontface = "bold") +
    geom_text_repel(data = specific_gene, 
                    aes(x = log2FoldChange, y = -log10(pvalue), label = Name_genes), 
                    color = "blue", size = 3, fontface = "bold") +
    labs(title = paste("Volcano Plot - Significant genes:", num_significant_genes), 
         x = "Log2FoldChange", 
         y = "-Log10(p-value)",
         color = "Significancia") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),  # Centrar el título
      axis.line = element_line(color = "black"), 
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    ) +
    coord_cartesian(xlim = c(-13, 13), ylim = c(0, 300))  # Ajustar los límites del gráfico
}

# Función para crear la tabla de anotación
create_annotation_table <- function(top_overexpressed, top_underexpressed, specific_gene) {
  annotation_data <- data.frame(
    Categoria = c("Sobreexpresados:", "Subexpresados:", "Gen Específico:"),
    Genes = c(paste(top_overexpressed$Name_genes, collapse = "\n"), 
              paste(top_underexpressed$Name_genes, collapse = "\n"), 
              paste(specific_gene$Name_genes, collapse = "\n"))
  )
  
  ggtexttable(
    annotation_data, rows = NULL,
    theme = ttheme("minimal", base_size = 9)  # Reducir el tamaño de la fuente para hacer la tabla más pequeña
  ) +
    theme(plot.background = element_rect(color = "grey", size = 1))
}

# Función para combinar el Volcano Plot y la tabla de anotación
combine_plots <- function(volcano_plot, annotation_plot) {
  ggdraw() +
    draw_plot(volcano_plot) +
    draw_plot(annotation_plot, x = 0.75, y = 0.35, width = 0.2, height = 0.2)  # Ajustar tamaño y posición
}

# Identificar los top 3 sobreexpresados y subexpresados
top_genes <- identify_top_genes(results_df)
top_overexpressed <- top_genes$overexpressed
top_underexpressed <- top_genes$underexpressed

# Cálculo de genes significativos
num_significant_genes <- calculate_significant_genes(results_df)

# Creación de los gráficos
volcano_plot <- create_volcano_plot(results_df, top_overexpressed, top_underexpressed, specific_gene, num_significant_genes)
annotation_plot <- create_annotation_table(top_overexpressed, top_underexpressed, specific_gene)
combined_plot <- combine_plots(volcano_plot, annotation_plot)

# Guardar el gráfico
ggsave("volcano_plot_with_annotation.png", plot = combined_plot, width = 12, height = 10)

# Mostrar el gráfico
print(combined_plot)


##############################################################
library(ggplot2)
library(dplyr)
library(ggpubr)

# Leer los archivos
go_bp <- read.table("/home/administrador/Escritorio/ontologias/sub/GO_Biological_Process_2023_table.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
go_cc <- read.table("/home/administrador/Escritorio/ontologias/sub/GO_Cellular_Component_2023_table.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
go_mf <- read.table("/home/administrador/Escritorio/ontologias/sub/GO_Molecular_Function_2023_table.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Filtrar términos relevantes y ordenar por p-valor
top_terms_bp <- go_bp %>%
  arrange(`P.value`) %>%
  head(10)

top_terms_cc <- go_cc %>%
  arrange(`P.value`) %>%
  head(10)

top_terms_mf <- go_mf %>%
  arrange(`P.value`) %>%
  head(10)

# Filtrar por p-value ajustado menor a 0.01 y por términos relevantes para ADAR
top_terms_bp <- go_bp %>%
  filter(`Adjusted.P.value` < 0.01) %>%
  filter(grepl("Inflammatory Response|Cytokine-Mediated Signaling Pathway|Regulation Of Interleukin-6 Production|Positive Regulation Of T Cell Proliferation", Term, ignore.case = TRUE)) %>%
  arrange(`Adjusted.P.value`)

# Mostrar los términos filtrados y ordenados
print(top_terms_bp)


# Filtrar por p-value ajustado menor a 0.15 y por términos relevantes para ADAR
top_terms_cc <- go_cc %>%
  filter(`Adjusted.P.value` < 0.15) %>%
  filter(grepl("Neuron Projection|Axonal Growth Cone|G Protein-Coupled Receptor Dimeric Complex|Synaptic Vesicle Membrane|Postsynaptic Density Membrane", Term, ignore.case = TRUE)) %>%
  arrange(`Adjusted.P.value`)

# Mostrar los términos filtrados y ordenados
print(top_terms_cc)


# Ver los términos filtrados y ordenados
print(relevant_terms_cc)

# Filtrar por p-value ajustado menor a 0.05 y por términos relevantes para ADAR
top_terms_mf <- go_mf %>%
  filter(`Adjusted.P.value` < 0.05) %>%
  filter(grepl("Cytokine Activity|Chemokine Activity|Chemokine Receptor Binding|CXCR Chemokine Receptor Binding", Term, ignore.case = TRUE)) %>%
  arrange(`Adjusted.P.value`)

# Mostrar los términos filtrados y ordenados
print(top_terms_mf)




library(patchwork)

# Crear los gráficos de barras con la paleta de colores lavanda a púrpura y etiquetas A, B, C
plot_bp <- ggplot(top_terms_bp, aes(x = reorder(Term, -`P.value`), y = -log10(`P.value`), fill = -log10(`P.value`))) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lavender", high = "purple") +
  coord_flip() +
  labs(x = "Procesos Biológicos", y = "-log10(P value)") +
  theme_pubr() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1)
  ) +
  scale_y_continuous(expand = c(0, 0))

plot_cc <- ggplot(top_terms_cc, aes(x = reorder(Term, -`P.value`), y = -log10(`P.value`), fill = -log10(`P.value`))) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lavender", high = "purple") +
  coord_flip() +
  labs(x = "Componentes Celulares", y = "-log10(P value)") +
  theme_pubr() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1)
  ) +
  scale_y_continuous(expand = c(0, 0))

plot_mf <- ggplot(top_terms_mf, aes(x = reorder(Term, -`P.value`), y = -log10(`P.value`), fill = -log10(`P.value`))) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lavender", high = "purple") +
  coord_flip() +
  labs(x = "Funciones Moleculares", y = "-log10(P Value)") +
  theme_pubr() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.y = element_text(size = 10, hjust = 1)
  ) +
  scale_y_continuous(expand = c(0, 0))

# Combinar los gráficos en una cuadrícula con etiquetas A, B, C y añadir un título general
combined_plot <- wrap_plots(plot_bp, plot_cc, plot_mf, ncol = 1) +
  plot_annotation(
    title = "Análisis de Enriquecimiento de GO\nGenes Subexpresados (N=421)",
    tag_levels = 'A',  # Agrega etiquetas A, B, C
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

# Mostrar el gráfico final
print(combined_plot)



##############################################################

library(ggplot2)
library(dplyr)
library(ggpubr)

# Leer los archivos
go_bp <- read.table("/home/administrador/Escritorio/ontologias/over/GO_Biological_Process_2023_table (1).txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
go_cc <- read.table("/home/administrador/Escritorio/ontologias/over/GO_Cellular_Component_2023_table.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
go_mf <- read.table("/home/administrador/Escritorio/ontologias/over/GO_Molecular_Function_2023_table.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Filtrar por p-value ajustado menor a 0.07 y por términos relevantes para ADAR
top_terms_bp <- go_bp %>%
  filter(`Adjusted.P.value` < 0.07) %>%
  filter(grepl("Adenosine To Inosine Editing|Regulation Of T Cell Differentiation In Thymus|DNA Methylation|Response To Stress|Regulation Of Transcription From RNA Polymerase II Promoter In Response To Stress|Positive Regulation Of Interleukin-6 Production", Term, ignore.case = TRUE)) %>%
  arrange(`Adjusted.P.value`)

# Mostrar los términos filtrados y ordenados
print(top_terms_bp)


top_terms_cc <- go_cc %>%
  arrange(`Adjusted.P.value`)

# Filtrar por p-value ajustado menor a 0.20 y por términos relevantes para ADAR
top_terms_mf <- go_mf %>%
  filter(`Adjusted.P.value` < 0.20) %>%
  filter(grepl("Adenosine Deaminase Activity|tRNA-specific Adenosine Deaminase Activity|Double-Stranded RNA Binding|Transcription Regulatory Region Nucleic Acid Binding|Chromatin DNA Binding", Term, ignore.case = TRUE)) %>%
  arrange(`Adjusted.P.value`)

# Mostrar los términos filtrados y ordenados
print(top_terms_mf)

library(patchwork)

# Crear los gráficos de barras sin títulos individuales
plot_bp <- ggplot(top_terms_bp, aes(x = reorder(Term, -`P.value`), y = -log10(`P.value`), fill = -log10(`P.value`))) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
  coord_flip() +
  labs(x = "Procesos Biológicos", y = "-log10(P value)") +
  theme_pubr() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none",
    plot.title = element_blank(),  # Eliminamos el título
    axis.text.y = element_text(size = 10, hjust = 1),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  scale_y_continuous(expand = c(0, 0))

plot_cc <- ggplot(top_terms_cc, aes(x = reorder(Term, -`P.value`), y = -log10(`P.value`), fill = -log10(`P.value`))) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
  coord_flip() +
  labs(x = "Componentes Celulares", y = "-log10(P value)") +
  theme_pubr() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none",
    plot.title = element_blank(),  # Eliminamos el título
    axis.text.y = element_text(size = 10, hjust = 1),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  scale_y_continuous(expand = c(0, 0))

plot_mf <- ggplot(top_terms_mf, aes(x = reorder(Term, -`P.value`), y = -log10(`P.value`), fill = -log10(`P.value`))) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
  coord_flip() +
  labs(x = "Funciones Moleculares", y = "-log10(P value)") +
  theme_pubr() +
  theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none",
    plot.title = element_blank(),  # Eliminamos el título
    axis.text.y = element_text(size = 10, hjust = 1),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  scale_y_continuous(expand = c(0, 0))

# Combinar los gráficos sin títulos individuales y ajustar los márgenes
combined_plot <- plot_bp / plot_cc / plot_mf +
  plot_annotation(
    title = "Análisis de Enriquecimiento de GO\nGenes Sobreexpresados (N=17)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(t = 20, b = 20)))
  ) & 
  theme(plot.margin = margin(20, 20, 20, 20))

# Mostrar el gráfico final
print(combined_plot)


