library(dplyr)

# Función para crear la matriz de ediciones
create_matrix <- function(site){ 
  ref <- c()    
  alt <- c()
  
  for(i in 1:length(site)){
    ref_tmp <- as.numeric(strsplit(site[[i]], split=",")[[1]][1])
    alt_tmp <- as.numeric(strsplit(site[[i]], split=",")[[1]][2])
    ref <- append(ref, ref_tmp)
    alt <- append(alt, alt_tmp)
  }
  mat_editions <- matrix(c(alt, ref), 6, 2)  # Matriz de 6 muestras y 2 columnas (ALT y REF)
  mat_editions[is.na(mat_editions)] <- 0  # Reemplaza NA por 0
  mat_ed_group <- rbind(colSums(mat_editions[1:3,]), colSums(mat_editions[4:6,]))  # Agrupa en dos grupos
  return(mat_ed_group)
}    

# MAIN
# Lee el archivo de entrada

# Inicializa un dataframe para almacenar los resultados
fisher_out <- data.frame()

# Aplica el test de Fisher a cada fila del dataframe
for(i in 1:nrow(data)){
  mat_ed <- create_matrix(data[i, 10:15])  # Columna 10 a 15 son HEK_A1-3 y HEK_G1-3
  tmp_fisher <- fisher.test(mat_ed)  # Realiza el test de Fisher
  fisher_out <- rbind(fisher_out, data.frame(data[i,1:9], pvalue=tmp_fisher$p.value, odds=tmp_fisher$estimate[[1]], data[i,10:15]))
}

# Visualiza los primeros resultados
head(fisher_out)

# Guarda los resultados en un archivo de salida
write.table(fisher_out, "/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/Merge/formatted_fil_cleaned.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


for(i in 1:5){  # Solo para las primeras 5 filas
  print(create_matrix(data_filter[i, 10:15]))
}







# Función para crear la matriz de ediciones con manejo de NA
create_matrix <- function(site){ 
  ref <- c()    
  alt <- c()
  
  for(i in 1:length(site)){
    # Separa ALT y REF, manejando posibles NA
    values <- strsplit(site[[i]], split=",")[[1]]
    
    if(length(values) == 2) {
      alt_tmp <- as.numeric(values[1])
      ref_tmp <- as.numeric(values[2])
    } else {
      alt_tmp <- 1
      ref_tmp <- 1
    }
    
    # Si se detecta un NA, se reemplaza por 0
    if(is.na(alt_tmp)) alt_tmp <- 
    if(is.na(ref_tmp)) ref_tmp <- 1
    
    ref <- append(ref, ref_tmp)
    alt <- append(alt, alt_tmp)
  }
  
  # Matriz de 6 muestras y 2 columnas (ALT y REF)
  mat_editions <- matrix(c(alt, ref), 6, 2)  
  mat_editions[is.na(mat_editions)] <- 1  # Reemplaza cualquier NA por 0
  
  # Agrupa las ediciones por grupos: HEK_A1-3 y HEK_G1-3
  mat_ed_group <- rbind(colSums(mat_editions[1:3,]), colSums(mat_editions[4:6,]))
  return(mat_ed_group)
}



# Inicializa un dataframe para almacenar los resultados
fisher_out <- data.frame()

# Aplica el test de Fisher a cada fila del dataframe
for(i in 1:nrow(data)){
  mat_ed <- create_matrix(data[i, 10:15])  # Columna 10 a 15 son HEK_A1-3 y HEK_G1-3
  tmp_fisher <- fisher.test(mat_ed)  # Realiza el test de Fisher
  fisher_out <- rbind(fisher_out, data.frame(data[i,1:9], pvalue=tmp_fisher$p.value, odds=tmp_fisher$estimate[[1]], data[i,10:15]))
}

# Visualiza los primeros resultados
head(fisher_out)


data <- read.table("/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/Merge/formatted_fil_cleaned.txt", 
                    sep="\t", header=TRUE, stringsAsFactors=FALSE, 
                    check.names=FALSE, colClasses = "character")


write.table(fisher_out, "/home/administrador/Escritorio/Variantes/Filtered/ADAR/listo/Merge/fishers.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Filtrar variantes con p-value significativo
significant_variants <- fisher_out[fisher_out$pvalue < 0.05, ]
# Clasificar según el Odds Ratio
high_odds <- significant_variants[significant_variants$odds > 1, ]  # OR > 1
low_odds <- significant_variants[significant_variants$odds < 1, ]   # OR < 1
neutral_odds <- significant_variants[significant_variants$odds == 1, ] # OR == 1


head(fisher_out)
# Mostrar una fila específica en R
data[1, ]  # Por ejemplo, la primera fila


# Construir la tabla de contingencia
# Construir la tabla de contingencia
tabla <- matrix(c(1, 0, 2228, 1955), nrow = 2)

# Realizar el test de Fisher
resultado <- fisher.test(tabla)
print(resultado)

# Realizar el test de Fisher
resultado <- fisher.test(tabla)
print(resultado)

# Construir la tabla de contingencia
tabla <- matrix(c(18, 20, 41, 19), nrow = 2)

# Realizar el test de Fisher
resultado <- fisher.test(tabla)
print(resultado)



