
# Proyecto de Análisis de Eventos Genéticos

Este repositorio contiene los archivos de datos generados para el análisis de varios eventos genéticos, incluyendo el uso diferencial de exones (DEU), uso diferencial de transcritos (DTU), análisis de expresión diferencial de genes (DEA) y eventos de splicing. A continuación, se describen los archivos principales y su contenido.

## Archivos

### 1. `Master_Table.tsv`
Este archivo contiene la tabla maestra combinada que incluye los resultados de diferentes análisis como DEA, DTU, DEU, splicing diferencial, y la cantidad de ediciones por gen. Es el archivo principal que reúne toda la información relevante de los diferentes enfoques analíticos. Las columnas principales son:

- **Columnas**:  
  `gene_symbol`, `ensembl_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `n_transcripts`, `n_differential_transcripts`, `n_overexpressed_transcripts`, `n_underexpressed_transcripts`, `n_exons`, `n_differential_exons`, `n_overexpressed_exons`, `n_underexpressed_exons`, `n_events`, `n_differential_events`, `n_overexpressed_events`, `n_underexpressed_events`, `Ediciones`

### 2. `SIG_DEA.tsv`
Este archivo contiene los resultados significativos del análisis de expresión diferencial de genes (DEA). Cada fila representa un gen, con detalles sobre:

- **Columnas**:  
  `ensembl_id`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `Name_genes`

### 3. `SIG_DTU.tsv`
Este archivo incluye los resultados significativos del uso diferencial de transcritos (DTU). El archivo reporta los transcritos que muestran cambios significativos entre las condiciones comparadas. Las columnas incluidas son:

- **Columnas**:  
  `transcript_id`, `gene_symbol`, `logFC`, `logCPM`, `LR`, `PValue`, `FDR`

### 4. `Sig_DEU.tsv`
Este archivo presenta los resultados del uso diferencial de exones (DEU). Cada fila corresponde a un exón que ha mostrado una variación significativa entre las condiciones, junto con las métricas estadísticas asociadas. Las columnas principales son:

- **Columnas**:  
  `Geneid`, `gene_symbol`, `Chr`, `Start`, `End`, `Strand`, `Length`, `logFC`, `P.Value`, `FDR`

### 5. `Sig_Splicing.tsv`
Este archivo contiene los eventos de splicing significativo identificados en el análisis. Cada evento es caracterizado por el tipo de splicing (como salto de exón o retención de intrón) y su relevancia estadística. Las columnas principales incluyen:

- **Columnas**:  
  `GeneID`, `geneSymbol`, `chr`, `strand`, `longExonStart_0base`, `longExonEnd`, `shortES`, `shortEE`, `flankingES`, `flankingEE`, `IJC_SAMPLE_1`, `SJC_SAMPLE_1`, `IJC_SAMPLE_2`, `SJC_SAMPLE_2`, `IncFormLen`, `SkipFormLen`, `PValue`, `FDR`, `IncLevel1`, `IncLevel2`, `IncLevelDifference`, `Tipo`

## Uso de los Archivos

1. **Master_Table.tsv**: Se utiliza como la tabla unificada para realizar comparaciones integradas entre DEA, DTU, DEU y splicing, incluyendo la cantidad de ediciones por gen. Útil para análisis globales de los resultados.
   
2. **SIG_DEA.tsv**: Se puede usar para enfocar análisis en genes que muestran diferencias significativas en su expresión entre las condiciones experimentales.

3. **SIG_DTU.tsv**: Permite estudiar los transcritos que tienen un uso diferencial significativo y pueden estar asociados con cambios en la funcionalidad del gen.

4. **Sig_DEU.tsv**: Útil para explorar cómo los exones individuales están siendo utilizados diferencialmente entre condiciones, lo cual puede influir en la estructura final de las proteínas.

5. **Sig_Splicing.tsv**: Ayuda a comprender los cambios en los patrones de splicing que pueden estar relacionados con diferentes condiciones celulares, lo que puede afectar las isoformas proteicas resultantes.

## Contacto
Para más información o preguntas sobre el proyecto, por favor contacta al autor del proyecto.
