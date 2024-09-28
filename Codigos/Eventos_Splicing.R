# Cargar librerías necesarias
library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(grid)

# Definir la ruta base
base_path <- "/home/administrador/Escritorio/rmats_post/"

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
#write.csv(datos_jcec_comb, "/home/administrador/Escritorio/rnaseq/SPLICING.csv")

library(ggplot2)
library(gridExtra)

# Gráfico de Área para la Distribución de p-values con colores y etiquetas
plot1 <- ggplot(datos_jcec_comb, aes(x = PValue, fill = Tipo)) +
  geom_area(stat = "bin", bins = 100, alpha = 0.7, position = 'identity') +
  theme_minimal() +
  labs(title = "Distribución de p-values en Eventos de Splicing",
       x = "p-value",
       y = "Cantidad de Eventos",
       fill = "Tipo de Conteo",
       caption = "Todos los eventos de splicing") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ajustar el eje Y para comenzar en 0
  scale_x_continuous(expand = c(0, 0)) +  # Ajustar el eje X para categorías
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

# Gráfico de Área para la Distribución de FDR con colores y etiquetas
plot2 <- ggplot(datos_jcec_comb, aes(x = IncLevelDifference, fill = Tipo)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7, position = 'identity') +
  theme_minimal() +
  labs(title = "Distribución de PSI en Eventos de Splicing",
       x = "PSI (Percent Spliced In)",
       y = "Cantidad de Eventos",
       fill = "Tipo de Conteo",
       caption = "Todos los eventos de splicing con PSI diferente de 0") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-1, NA)) +  # Ajustar el eje X para comenzar en 0
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

# Gráfico de Histograma para la Distribución de FDR
plot3 <- ggplot(datos_jcec_comb, aes(x = FDR, fill = Tipo)) +
  geom_area(stat = "bin", bins = 100, alpha = 0.7, position = 'identity') +
  theme_minimal() +
  scale_x_continuous(trans='log10') +
  labs(title = "Distribución de FDR en Eventos de Splicing",
       x = "FDR (escala logarítmica)",
       y = "Cantidad de Eventos",
       fill = "Tipo de Conteo",
       caption = "Todos los eventos de splicing") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  # Ajustar el eje Y para comenzar en 0
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

# Gráfico de Barras para el Conteo de Tipos de Eventos de Splicing
plot4 <- ggplot(datos_jcec_comb, aes(x = Tipo, fill = Tipo)) +
  geom_bar(binwidth = 0.1, alpha = 0.7, position = 'identity') +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Conteo de Tipos de Eventos de Splicing (163.029)",
       x = "Tipo de Evento",
       y = "Conteo",
       fill = "Tipo de Evento",
       caption = "Fuente: datos JCEC combinados") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 120000)) +  # Ajustar el eje Y para comenzar en 0 y extender hasta 120000
  scale_x_discrete(expand = c(0, 0)) +  # Ajustar el eje X para categorías
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

# Título general para la cuadrícula
titulo_grid <- textGrob("Análisis de eventos de splicing", gp = gpar(fontsize = 16, fontface = "bold"))

# Combinar los gráficos en una cuadrícula
grid.arrange(
  titulo_grid,
  arrangeGrob(plot1, plot2, plot3, plot4, ncol = 2, padding = unit(1, "line")),
  ncol = 1,
  heights = unit.c(unit(1, "lines"), unit(1, "npc") - unit(1, "lines"))
)



datos_filtrados_Pvalue <- datos_jcec_comb %>% filter(PValue <= 0.05)
datos_filtrados_Pvalue_PSI <- datos_filtrados_Pvalue %>% filter(abs(IncLevelDifference) >= 0.1)
datos_filtrados_Pvalue_PSI_FDR <- datos_filtrados_Pvalue_PSI %>% filter(FDR <= 0.20)


library(UpSetR)

head(datos_filtrados_Pvalue_PSI_FDR)
# Filtrar los datos relevantes y crear una lista de eventos por gen
eventos_por_gen <- datos_filtrados_Pvalue_PSI_FDR %>%
  dplyr::select(geneSymbol, Tipo) %>%
  distinct() %>%
  group_by(geneSymbol) %>%
  summarise(eventos = list(Tipo))

total_genes <- datos_filtrados_Pvalue_PSI_FDR %>%
  dplyr::select(geneSymbol) %>%
  distinct() %>%
  nrow()

# Mostrar el total de genes
total_genes

ggplot(eventos_por_gen, aes(x = eventos)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +  # Agrega las etiquetas sobre las barras
  geom_point(stat = "count") +
  ggupset::scale_x_upset() +
  labs(
    title = "UpSet Plot de Eventos de Splicing",
    x = "Tipos de Eventos",
    y = "Número de Genes"
  ) +
  theme_minimal()



#write.csv(datos_filtrados_Pvalue_PSI_FDR, "/home/administrador/Escritorio/rnaseq/SIG_SPLICING.csv")

head(datos_filtrados_Pvalue_PSI_FDR)
library(ggplot2)
library(gridExtra)

# Gráfico de Área para la Distribución de p-values
plot1 <- ggplot(datos_filtrados_Pvalue, aes(x = PValue, fill = Tipo)) +
  geom_area(stat = "bin", bins = 100, alpha = 0.7, position = 'identity') +
  theme_minimal() +
  labs(title = paste("Distribución de p-values (Eventos: ", nrow(datos_filtrados_Pvalue), ")", sep = ""),
       x = "p-value",
       y = "Cantidad de Eventos",
       fill = "Tipo de Conteo",
       caption = "Todos los eventos de splicing") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))

# Gráfico de Área para la Distribución de PSI
plot2 <- ggplot(datos_filtrados_Pvalue_PSI, aes(x = IncLevelDifference, fill = Tipo)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7, position = 'identity') +
  theme_minimal() +
  labs(title = paste("Distribución de PSI (Eventos: ", nrow(datos_filtrados_Pvalue_PSI), ")", sep = ""),
       x = "PSI (Percent Spliced In)",
       y = "Cantidad de Eventos",
       fill = "Tipo de Conteo",
       caption = "Todos los eventos de splicing con PSI diferente de 0") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-1, NA))

# Gráfico de Histograma para la Distribución de FDR
plot3 <- ggplot(datos_filtrados_Pvalue_PSI_FDR, aes(x = FDR, fill = Tipo)) +
  geom_area(stat = "bin", bins = 100, alpha = 0.7, position = 'identity') +
  theme_minimal() +
  scale_x_continuous(trans='log10') +
  labs(title = paste("Distribución de FDR (Eventos: ", nrow(datos_filtrados_Pvalue_PSI_FDR), ")", sep = ""),
       x = "FDR (escala logarítmica)",
       y = "Cantidad de Eventos",
       fill = "Tipo de Conteo",
       caption = "Todos los eventos de splicing") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

# Gráfico de Barras para el Conteo de Tipos de Eventos de Splicing
datos_filtrados_Pvalue_PSI_FDR <- datos_filtrados_Pvalue_PSI_FDR %>% 
  filter(!is.na(Tipo))
scale_y_continuous(expand = c(0, 0))

plot4 <- ggplot(datos_filtrados_Pvalue_PSI_FDR, aes(x = Tipo, fill = Tipo)) +
  geom_bar(binwidth = 0.1, alpha = 0.7, position = 'identity') +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  theme_minimal() +
  labs(title = "Conteo de Tipos",
       x = "Tipo de Evento",
       y = "Conteo",
       fill = "Tipo de Evento",
       caption = "Fuente: datos JCEC combinados") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2300)) + 
  scale_x_discrete(expand = c(0, 0)) +  # Ajustar el eje X para categorías
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"))
plot4
# Título general para la cuadrícula
titulo_grid <- textGrob("Análisis de eventos de splicing post-filtrado", gp = gpar(fontsize = 16, fontface = "bold"))

# Combinar los gráficos en una cuadrícula
grid.arrange(
  titulo_grid,
  arrangeGrob(plot1, plot2, plot3, plot4, ncol = 2),
  ncol = 1,
  heights = unit.c(unit(1, "lines"), unit(1, "npc") - unit(1, "lines"))
)


# Instalar plotly desde Github para obtener gráficos de embudo si es necesario
# devtools::install_github("ropensci/plotly")
library(plotly)

# Crear un data frame con los conteos
conteos <- data.frame(
  Filtro = c("Total", "PValue <= 0.05", "PSI >= 0.1", "FDR <= 0.1"),
  Cantidad = c(nrow(datos_jcec_comb),
               nrow(datos_filtrados_Pvalue),
               nrow(datos_filtrados_Pvalue_PSI),
               nrow(datos_filtrados_Pvalue_PSI_FDR))
)



library(plotly)

# Crear un data frame con los conteos
conteos <- data.frame(
  Filtro = c("Splicing Brutos", "PSI >= 0.1", "FDR <= 0.2"),
  Cantidad = c(nrow(datos_jcec_comb),
               nrow(datos_filtrados_Pvalue_PSI),
               nrow(datos_filtrados_Pvalue_PSI_FDR))
)

# Escalar los valores para mejorar la visibilidad
conteos$ScaledCantidad <- sqrt(conteos$Cantidad) * 10

# Crear el gráfico de embudo de área con valores escalados y etiquetas con la cantidad de conteos
fig <- plot_ly(
  type = "funnelarea",
  values = conteos$ScaledCantidad,
  text = paste(conteos$Filtro, ": ", conteos$Cantidad),
  marker = list(colors = c("#a0e2b36D", "#82c2946D", "#dcb0eaFF"),
                line = list(color = c("wheat", "wheat", "wheat"), width = c(1, 1, 1))),
  textinfo = "text",
  textfont = list(family = "Old Standard TT, serif", size = 20, color = "black"),
  opacity = 0.9
)

fig <- fig %>% layout(
  title = list(text = "Conteo de eventos filtrados", font = list(size = 17)),
  funnelmode = "stack",
  yaxis = list(showgrid = FALSE, zeroline = FALSE),
  xaxis = list(showgrid = FALSE, zeroline = FALSE)
)

fig

# Cargar los paquetes
library(webshot)
library(webshot2)

# Guardar el gráfico como archivo HTML
htmlwidgets::saveWidget(as_widget(fig), "funnel_plot.html")




#####################################################################################################################################





