#caracte

# Cargar el archivo VCF
vcf_data <- read.table("/home/administrador/Escritorio/MASKER/combined.vcf", sep="\t", header=FALSE, stringsAsFactors=FALSE, comment.char = "#")
colnames(vcf_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "A1", "A2", "A3", "G1", "G2", "G3")

# Ver las primeras filas para confirmar la carga correcta
head(vcf_data)

is_exclusive_to_A <- function(A1, A2, A3, G1, G2, G3) {
  # Verifica si al menos una de las muestras A tiene un genotipo no homocigoto para el alelo de referencia (es decir, 0/1 o 1/1)
  in_A <- any(grepl("1/1|0/1", c(A1, A2, A3)))
  
  # Verifica si todas las muestras G son estrictamente homocigotas para el alelo de referencia (es decir, 0/0)
  not_in_G <- all(grepl("0/0|\\./\\.", c(G1, G2, G3)))
  
  # Retorna TRUE si la variante es exclusiva de las muestras A
  return(in_A & not_in_G)
}

# Aplicar la función a cada fila del dataframe
exclusive_variants <- vcf_data[apply(vcf_data[, c("A1", "A2", "A3", "G1", "G2", "G3")], 1, function(x) is_exclusive_to_A(x[1], x[2], x[3], x[4], x[5], x[6])), ]

# Ver los resultados
head(exclusive_variants)


write.table(exclusive_variants, "/home/administrador/Escritorio/MASKER/Solo_A.vcf", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)




# Función para identificar si una variante es compartida entre A y G
is_shared <- function(A1, A2, A3, G1, G2, G3) {
  # Verifica si al menos una muestra en ambos grupos A y G tiene un genotipo no homocigoto para el alelo de referencia
  in_A <- any(grepl("1/1|0/1", c(A1, A2, A3)))
  in_G <- any(grepl("1/1|0/1", c(G1, G2, G3)))
  
  # Retorna TRUE si la variante está presente en ambos grupos
  return(in_A & in_G)
}

# Aplicar la función para identificar variantes compartidas
shared_variants <- vcf_data[apply(vcf_data[, c("A1", "A2", "A3", "G1", "G2", "G3")], 1, function(x) is_shared(x[1], x[2], x[3], x[4], x[5], x[6])), ]

# Ver cuántas variantes son compartidas
num_shared_variants <- nrow(shared_variants)
print(paste("Número de variantes compartidas:", num_shared_variants))
write.table(shared_variants, "/home/administrador/Escritorio/MASKER/Compartidas_A_G.vcf", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)



# Función para identificar si una variante es exclusiva del grupo G
is_exclusive_to_G <- function(A1, A2, A3, G1, G2, G3) {
  # Verifica si al menos una de las muestras G tiene un genotipo no homocigoto para el alelo de referencia
  in_G <- any(grepl("1/1|0/1", c(G1, G2, G3)))
  
  # Verifica si todas las muestras A son homocigotas para el alelo de referencia (es decir, 0/0) o no tienen datos (./.)
  not_in_A <- all(grepl("0/0|\\./\\.", c(A1, A2, A3)))
  
  # Retorna TRUE si la variante es exclusiva del grupo G
  return(in_G & not_in_A)
}

# Aplicar la función para identificar variantes exclusivas de G
exclusive_variants_G <- vcf_data[apply(vcf_data[, c("A1", "A2", "A3", "G1", "G2", "G3")], 1, function(x) is_exclusive_to_G(x[1], x[2], x[3], x[4], x[5], x[6])), ]

# Ver cuántas variantes son exclusivas de G
num_exclusive_variants_G <- nrow(exclusive_variants_G)
print(paste("Número de variantes exclusivas del grupo G:", num_exclusive_variants_G))
write.table(exclusive_variants_G, "/home/administrador/Escritorio/MASKER/Solo_G.vcf", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)





# Contar variantes totales
total_variants <- nrow(vcf_data)

# Contar variantes compartidas
num_shared_variants <- nrow(shared_variants)

# Contar variantes exclusivas de G
num_exclusive_variants_G <- nrow(exclusive_variants_G)

# Contar variantes exclusivas de A
num_exclusive_variants_A <- total_variants - (num_shared_variants + num_exclusive_variants_G)


# Datos para el diagrama de torta
variant_counts <- c(num_shared_variants, num_exclusive_variants_A, num_exclusive_variants_G)
total_variants <- sum(variant_counts)
percentages <- round(variant_counts / total_variants * 100, 1)

# Crear etiquetas con cantidad y porcentaje
labels <- paste(c("Compartidas", "Exclusivas de A", "Exclusivas de G"), 
                "\n", variant_counts, " variantes\n", percentages, "%", sep="")

# Crear el diagrama de torta con etiquetas y mejorar la presentación
pie(variant_counts, 
    labels = labels, 
    main = "Distribución de Variantes en Grupos A y G", 
    col = c("skyblue", "orange", "lightgreen"), 
    cex = 0.8,        # Tamaño de las etiquetas
    radius = 1,       # Tamaño del gráfico
    border = "white", # Borde blanco para mejor contraste
    clockwise = TRUE) # Girar en sentido horario
