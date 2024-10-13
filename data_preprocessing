preprocessing_script.R

### 01. Carga de paquetes/librerías
```{r include=FALSE}
library(BiocManager)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(dplyr)
library(NOISeq)
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(ggbiplot)
```

# 02. Este código es para obtener información sobre los genes del genoma humano relacionados al cáncer de mama mediante la base de datos de ensembl, en resumen, este código utiliza el paquete biomaRt de R para obtener información sobre los genes relacionados con el cáncer de mama a partir de la base de datos Ensembl del genoma humano. La información obtenida incluye el ID del gen, su posición en el cromosoma, su longitud, su símbolo HGNC, su contenido de GC y su tipo de gen.
```{r echo=FALSE}
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www") # Usa Ensembl del paquete Biomart para conectar la bade deatos de Ensembl

features <- c("ensembl_gene_id", "chromosome_name", # Se especifica las características que se desean obtener de Ensembl para cada gen.
              "start_position", "end_position", "hgnc_symbol",	
              "percentage_gene_gc_content", "gene_biotype")

chrs <- c(1:22, "X", "Y") # Se especifican los cromosomas de donde se obtiene info

annot <- getBM(attributes = features, # Uso de getBM para obtener información sobre los genes que se encuentran en los cromosomas 'chrs'
      filters = "chromosome_name",
      values = chrs, 
      mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type") # Cambio de nombres
annot$Length <- abs(annot$End - annot$Start) # Colocación de una columna llamada Length
```

# 03. Este código se usa para descargar los datos de expresión génica (en forma de recuentos de expresión génica) de un estudio cáncer de mama del proyecto TCGA (The Cancer Genome Atlas) a través de la plataforma GDC (Genomic Data Commons), es decir, utiliza la plataforma GDC para descargar datos de expresión génica de TCGA para cáncer de mama y preparar los datos descargados para su análisis posterior. 

```{r include=FALSE}
query <- GDCquery(project = "TCGA-BRCA", # Se especifican parámetros de búsqueda para datos de expresión de TCGA-BRCA
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts", )

samplesDown <- getResults(query,cols=c("cases")) # getResults obtiene los casos ('cases') de las muestras seleccionadas en query

# Como no se pudo descargar toda la data se realiza la descarga en 2 partes (2 objetos: half y other_half, respectivamente para los 1231 casos)
half <- samplesDown[1:553] # Primeros 553 muestras (casos)
other_half <- samplesDown[554:1231] # Los casos restantes

# Se debe indicar que ahora se especifican los parámetros para cada data descargada (half y other_half)
queryDown_half <- GDCquery(project = "TCGA-BRCA", # queryDown_half especifica los mismos paramtros de busqueda que query pero ahora solo con half
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts",  
                      barcode = half)

queryDown_other_half <- GDCquery(project = "TCGA-BRCA", # lo mismo que el anterior pero con other_half
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts",  
                      barcode = other_half)

GDCdownload(query = queryDown_half) # Descarga los datos de expresión génica queryDown_half (con todos los parámetros especificados)
GDCdownload(query = queryDown_other_half) #  Descarga los datos de expresión génica queryDown_other_half

dataPrep1_1 <- GDCprepare(query = queryDown_half) # Prepara los datos en el environment para queryDown_half con el nombre de 'dataPrep1_1' 
dataPrep1_2 <- GDCprepare(query = queryDown_other_half) # # Prepara los datos en el environment para queryDown_other_half con el nombre de 'dataPrep1_2'
```

# 04. El siguiente código se utiliza para seleccionar los datos de expresión génica correspondientes a tumores primarios y tejidos normales en el estudio del cáncer de mama de TCGA. Los datos seleccionados se guardan en cuatro variables diferentes para cada uno de los dos grupos.

```{r}
gc() # Libera memoria
# De dataPrep_1_1 filtra datos para segun 'sample_type' o sea, Primary Tumor o Normal (para diferenciar segun sample_type si es Primary Tumor o Solid Tissue Normal)
BC1.T <- dataPrep1_1[ , dataPrep1_1$sample_type == "Primary Tumor"]
BC1.N <- dataPrep1_1[ , dataPrep1_1$sample_type == "Solid Tissue Normal"]

# De dataPrep_1_2 filtra datos para segun el sample type sea, Primary Tumor o Normal, igual que al anterior
BC2.T <- dataPrep1_2[ , dataPrep1_2$sample_type == "Primary Tumor"]
BC2.N <- dataPrep1_2[ , dataPrep1_2$sample_type == "Solid Tissue Normal"]
```
# 05. SOLO SI ES NECESARIO, (sin revisión de pureza 0,0,0,0,0) - El siguiente código se utiliza para estimar la pureza tumoral de las muestras de tumor primario del cáncer de mama utilizando la función TCGAtumor_purity() y luego guardar los códigos de barras de las muestras con una pureza estimada mayor o igual a 0.6 en la variable Purity.BRCA.1. 
```{r}
Purity.BRCA.1<-TCGAtumor_purity(colnames(BC1.T), 0,0,0,0,0)$pure_barcodes # Evalúa los tumores (BC1.T) que cumplan con criterios de pureza y pure_codes indica cuáles son esas muestras.
### 424 al filtrar pureza
```

```{r}
Purity.BRCA.2<-TCGAtumor_purity(colnames(BC2.T), 0,0,0,0,0)$pure_barcodes # # Evalúa los tumores (BC2.T) que cumplan con criterios de pureza y pure codes indica cuáles son esas muestras.
### 516 al filtrar pureza
```

# 06. Calcula las muestras de tumores primarios en el conjunto half que no pertenecen a los subtipos moleculares filtrados por la función TCGA_MolecularSubtype(). La función TCGA_MolecularSubtype() se utiliza para clasificar los tumores en diferentes subtipos moleculares en función de la expresión génica.

```{r}
# diff.1 y diff.2 son vectores que contienen los códigos de barras de las muestras que cumplen con el criterio de pureza especificado pero no cumplen con los criterios de subtipo molecular especificados.
diff.1 <- setdiff(Purity.BRCA.1, # Da como resultado un vector que contiene elementos que están en Purity.BRCA.1 pero no en el resultado de la función TCGA_MolecularSubtype(half)$filtered.
                  TCGA_MolecularSubtype(half)$filtered) 

diff.2 <- setdiff(Purity.BRCA.2, # Da como resultado un vector que contiene elementos que están en Purity.BRCA.2 pero no en el resultado de la función TCGA_MolecularSubtype(half)$filtered.
                  TCGA_MolecularSubtype(other_half)$filtered) 
```

# 07. Se crea una matriz llamada 'rnas' mediante la función cbind() que une las matrices de expresión génica assay() correspondientes a las muestras tumorales y normales para ambos conjuntos de datos, que se seleccionan mediante el uso de [,diff.1] y [,diff.2].
```{r}
rnas <- cbind(assay(BC1.T)[,diff.1], assay(BC2.T)[,diff.2], assay(BC1.N), assay(BC2.N))
head(rnas)
dim(rnas)
# Específicamente de BC1.T se selecciona y luego se extrae la información de expresión génica de las columnas correspondientes a barcodes que estén en diff.1 (todos los colnames presentes en diff.1), lo mismo para BC2.T
gc()
```

# 08. Se crea un objeto llamado mol_subtypes que contiene información sobre los subtipos moleculares de los tumores y las muestras normales.
```{r}
mol_subtypes.1 <-TCGA_MolecularSubtype(colnames((BC1.T)[,diff.1]))$subtypes$subtype # Identifican subtipos moleculares de tumores filtrados del conjunto de datos BC1.T
mol_subtypes.2 <-TCGA_MolecularSubtype(colnames((BC2.T)[,diff.2]))$subtypes$subtype # Identifican subtipos moleculares de tumores filtrados del conjunto de datos BC2.T
normal = data.frame(subtype = rep('normal', 113)) # crea un marco de datos con una columna "subtype" que contiene la etiqueta "normal" 
mol_subtypes.1 <- data.frame(subtype = mol_subtypes.1) # Se convierten en dataframes para ser unidos posteriormente
mol_subtypes.2 <- data.frame(subtype = mol_subtypes.2)

mol_subtypes = rbind(mol_subtypes.1, mol_subtypes.2, normal) # Unión (rbind) de subtipos moleculares (1 y 2) y normal en un dataframe que llamado mol_subtypes
table(mol_subtypes) 
mol_subtypes <- make.names(mol_subtypes) # Hace que los nombres de los subtipos moleculares sean legales y se formateen 
head(mol_subtypes)
gc()
```

# 09. En este bloque de código se prepara un dataframe llamado "factors" que contiene información sobre las muestras, su grupo y subtipo molecular.
```{r}
# Crea un dataframe llamado factorBC con dos columnas llamadas Group y Sample. La columna Group contiene la cadena "PT" (que significa "Primary Tumor") y la columna Sample contiene los nombres de las columnas de la matriz assay(BC1.T) y assay(BC2.T).
factorBC <- data.frame(Group = "PT", Sample = c(colnames(assay(BC1.T)[,diff.1]), colnames(assay(BC2.T)[,diff.2])))
# Crea un dataframe llamado factorsNormalBC que contiene dos columnas: "Group" y "Sample". La columna "Group" tiene el valor "Normal" y la columna "Sample" contiene los nombres de las muestras de tejidos normales
factorsNormalBC <- data.frame(Group = "Normal", Sample = c(colnames(BC1.N),colnames(BC2.N))) 
factors <- rbind(factorBC, factorsNormalBC) # Combinación de los dataframes creados
factors = cbind(factors, mol_subtypes) # Se añade la columna mol_subtypes
rownames(factors) <- factors$Sample 
Ready_factors <- as.data.frame(factors$Group) # Ready_factors que contiene únicamente la información de la columna "Group" de factors
```

# 10. Este código realiza un filtrado de datos de expresión génica en una matriz "rnas" mediante el método "quantile". El objetivo es eliminar genes que presentan una expresión muy baja en una cantidad significativa de muestras.
```{r}
dataFilt <- TCGAanalyze_Filtering(tabDF = rnas, # Aplica filtrado "quantile" sobre la matriz de expresión génica 'rnas'
                                  method = "quantile",
                                  qnt.cut = 0.25)
dim(dataFilt)
threshold <- round(dim(rnas)[2]/2) # Calcula la mitad del número total de columnas en rnas, y luego redondea al número entero más cercano. 
threshold
ridx <- rowSums(dataFilt == 0) <= threshold # Se indica qué filas de dataFilt tienen valores deben mantenerse en función de treshold
table(ridx)
dataFilt <- dataFilt[ridx, ] # Selección de filas que cumplen con ridx
ridx <- rowMeans(dataFilt) >= 50 # Filtro para quedarse con valor de expresión más altos según las medias de los genes
table(ridx)
dataFilt <- dataFilt[ridx, ] # Selección de filas que cumplen con ridx
print(dim(dataFilt))
rnas <- rnas[rownames(rnas) %in% rownames(dataFilt), ] # filtra la matriz rnas para mantener solo las filas correspondientes a las muestras que también se encuentran en la matriz dataFilt.
dim(rnas) # 'rnas' contiene los valores de expresión filtrados por las condiciones descritas líneas arriba
```
# 11. Este código se utiliza para hacer coincidir los identificadores de genes en dos conjuntos de datos diferentes (rnas y annot) y crear subconjuntos de datos filtrados (rnas1 y annot1) con los que se puede trabajar posteriormente. Hace check para los posibles duplicados
```{r}
rownames(rnas) = gsub("\\..*","", rownames(rnas)) # Elimina extensión de archivos de valores de filas de 'rnas' y les da formato
inter <- intersect(rownames(rnas), annot$ensembl_gene_id) # Intersección entre rownames(rnas) y ensembl_gene_id de annot (annot$ensembl_gene_id)
length(inter)
rnas1 <- rnas[rownames(rnas) %in% inter,] # Selecciona las filas de 'rnas' que se encuentran en inter.
dim(rnas1)
annot1 <- annot[annot$ensembl_gene_id  %in% inter,] # Selecciona las filas de rnas de los genes presentes en annot$ensembl_gene_id de inter
dim(annot1)
annot1 <- annot1[!duplicated(annot1$ensembl_gene_id),] # Mantiene las filas únicas en "annot1" basadas en la columna "ensembl_gene_id"
dim(annot1)
annot1[annot1 == ""] <- NA # Asigna valores "NA" a todas las celdas en "annot1" donde el valor es una cadena vacía ""
gc()
```

# 12. Este código normaliza y corrige los datos de expresión génica utilizando una combinación de métodos de normalización de dentro y entre carriles (dentro y entras las muestras), y finalmente aplica el algoritmo NOISeq para detectar genes diferencialmente expresados (Para PT y Normal)
```{r}
ln.data <- withinLaneNormalization(rnas1, annot1$Length, which = "full") # Se normaliza por la longitud del gen
gcn.data <- withinLaneNormalization(ln.data , annot1$GC, which = "full") # Normalización por contenido de GC
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") # Normalización entre carriles (muestras) para complementar 'gcn.data'
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0) # Normalización TMM 
norm.counts
gc()
noiseqData <- NOISeq::readData(norm.counts, factors = Ready_factors) # Se especifican los factores a comparar para los DEGs
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE) # Normalizacion de datos por ARSyNSeq
rnas2 <- exprs(mydata2corr1) # rnas2 contiene los datos normalizados y corregidos de expresión génica de los datos originales de RNA-seq.
```

# 12.00 Creación de Ready_factorssubtype
```{r}
# Obtener subtipos moleculares para las muestras de tumor en 'factors'
subtypes_tumor <- ifelse(factors$Group == "PT", TCGA_MolecularSubtype(factors$Sample[factors$Group == "PT"])$subtypes$subtype, NA)

# Crear una nueva columna 'subtype' en 'factors'
factors$subtype <- NA

# Asignar los subtipos a las muestras de tumor
factors$subtype[factors$Group == "PT"] <- subtypes_tumor

# Etiquetar específicamente las muestras normales
factors$subtype[factors$Group == "Normal"] <- "Normal"

# Crear 'factors1' con la información correcta de subtipos
factors1 <- factors

# Visualizar 'factors1'
head(factors1)

# Filtrar muestras con subtipos faltantes
factors_filtered <- factors[!is.na(factors$subtype), ]
rnas1_filtered <- rnas1[rownames(rnas1) %in% factors_filtered$Sample, ]

# Crear un nuevo objeto 'Ready_factorssubtype' con los subtipos filtrados
Ready_factorssubtype_filtered <- as.data.frame(factors_filtered$subtype)
```

# 12.01 Preparación de Ready_factorssubtype - para que coja los subtipos en lugar de groups (PT y Normal)
```{r}
Ready_factorssubtype <- as.data.frame(factors_filtered$subtype)
subtypes <- levels(factors_filtered$subtype)
# noiseqData <- NOISeq::readData(norm.counts, factors = Ready_factorssubtype)

# Filtrar 'rnas1' para que contenga solo las muestras presentes en 'factors_filtered'
rnas1_filtered <- rnas1[, colnames(rnas1) %in% factors_filtered$Sample]
```

# 12.1 Corregido para normalizar en NOISeq con SUBTIPOS en lugar de group
```{r}
ln.data <- withinLaneNormalization(rnas1_filtered, annot1$Length, which = "full") # Se normaliza por la longitud del gen
gcn.data <- withinLaneNormalization(ln.data , annot1$GC, which = "full") # Normalización por contenido de GC
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") # Normalización entre carriles (muestras) para complementar 'gcn.data'
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0) # Normalización TMM 
norm.counts
gc()
noiseqData <- NOISeq::readData(norm.counts, factors = Ready_factorssubtype) # Se especifican los factores a comparar para los DEGs
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE) # Normalizacion de datos por ARSyNSeq
rnas2 <- exprs(mydata2corr1) # rnas2 contiene los datos normalizados y corregidos de expresión génica de los datos originales de RNA-seq.
gc()
```

### Quality control 
# 13. Control de calidad - El código utiliza la librería ggbiplot para realizar un análisis de componentes principales (PCA) antes y después de la normalización de los datos de expresión génica, incluido (BRCA.Normal)
```{r}
gc()
library(ggbiplot) 
before.pca <- prcomp(t(rnas1_filtered),center = TRUE,scale. = TRUE) # Análisis de componentes principales a rnas1 (antes de ARSyNSeq) transpuesta
summary(before.pca) 
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=factors_$Group) # Visualización de análisis pca a before.pca
 
after.pca <- prcomp(t(rnas2),center = TRUE,scale. = TRUE) # Análisis de componentes principales a rnas2 transpuesta (después de ARSyNSeq)
summary(after.pca) 
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=factors$Group) # Vidualización de análisis pca a after.pca 
```

# 13.1 Normalizado a Subtipos moleculares (Todos los subtipos, incluido BRCA.Normal)
```{r}
# PCA para subtipos moleculares (Con todos los subtipos moleculares)
before.pca <- prcomp(t(rnas1_filtered),center = TRUE,scale. = TRUE) # Análisis de componentes principales a rnas1 (antes de ARSyNSeq) transpuesta
summary(before.pca) 
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=factors_filtered$subtype) # Visualización de análisis pca a before.pca
 
after.pca <- prcomp(t(rnas2),center = TRUE,scale. = TRUE) # Análisis de componentes principales a rnas2 transpuesta (después de ARSyNSeq)
summary(after.pca) 
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=factors_filtered$subtype) # Vidualización de análisis pca a after.pca 
```

# 13.2 PCA para subtipos moleculares de cancer de mama
```{r}
# PCA para subtipos moleculares 2
########################################################## BEFORE
# Definir los subtipos seleccionados
selected_subtypes <- c("BRCA.Basal", "BRCA.Her2", "BRCA.LumA", "BRCA.LumB", "Normal")

# Crear una matriz de datos de expresión para los subtipos seleccionados
selected_data <- rnas1_filtered[, factors_filtered$subtype %in% selected_subtypes]

# Realizar análisis de PCA
before.pca <- prcomp(t(selected_data), center = TRUE, scale. = TRUE)

# Crear un vector de grupos para los subtipos seleccionados
subtype_groups <- factor(factors_filtered$subtype[factors_filtered$subtype %in% selected_subtypes])

# Visualizar el gráfico de PCA
pca_before <- ggbiplot(before.pca, var.axes = FALSE, ellipse = TRUE, groups = subtype_groups)

# Guardar el gráfico de PCA
ggsave('pca_before.jpg', pca_before, dpi = 300)

############################################################ AFTER
# Definir los subtipos seleccionados
selected_subtypes <- c("BRCA.Basal", "BRCA.Her2", "BRCA.LumA", "BRCA.LumB", "Normal")

# Crear una matriz de datos de expresión para los subtipos seleccionados en 'rnas2'
selected_data <- rnas2[, factors_filtered$subtype %in% selected_subtypes]

# Realizar análisis de PCA en los datos seleccionados
after.pca <- prcomp(t(selected_data), center = TRUE, scale. = TRUE)

# Crear un vector de grupos para los subtipos seleccionados
subtype_groups <- factor(factors_filtered$subtype[factors_filtered$subtype %in% selected_subtypes])

# Visualizar el gráfico de PCA
pca_after <- ggbiplot(after.pca, var.axes = FALSE, ellipse = TRUE, groups = subtype_groups)

# Guardar el gráfico de PCA
ggsave("pca_after.jpg", pca_after, dpi = 300)
```






































  
