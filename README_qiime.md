# *Qiime2 homework*
## Víctor Manuel López Molina
## Fecha: 20/05/2022
## Ejercicio

En este ejercicio de la asignatura, dados una serie de ficheros y carpetas proporcionados por el profesor, vamos a redactar un informe en el que se muestren los comandos, con sus correspondientes parámetros, ejecutados por el profesor para obtener dichos archivos. Adicionalmente, junto con este informe, se rellenará el cuestionario disponible en Moodle. 

## Informe

### 1. Importamos los datos en qiime2

En primer lugar, debemos recordar que, para esta práctica, se han tomado datos de un experimento en el que se analizó la influencia del frío sobre la microbiota intestinal en ratones, y la secuenciación realizada fue rDNA 16S. Dicho esto, vamos a importar los datos proporcionados en qiime, por medio del fichero llamado "samplemanifest", el cual nos da información acerca de la localización de cada archivo y el identificador de la muestra correspondiente, por medio del siguiente comando:
```
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
                   --input-path samplemanifest \
                   --output-path paired-end-demux.qza \
                   --input-format PairedEndFastqManifestPhred33
```
Este comando generará el primer artefacto y, tras este paso, vamos a generar un resumen del proceso de importación. Este nos permitirá conocer cuántas secuencias fueron obtenidas por muestra, así como la distribución de las calidades de la secuencia por posición. Para ello, ejecutamos esto:
```
qiime demux summarize --i-data paired-end-demux.qza --o-visualization paired-end-demux.qzv
```
Este fichero .qzv puede ser visualizado por medio del comando:
```
qiime tools view paired-end-demux.qzv
```
Por medio de este fichero podemos, efectivamente, conocer el número de muestras del experimento, cuántas secuencias existen por muestra, el número total de secuencias, la distribución de las calidades por posición, etc.

### 2. Determinamos los ASVs por medio de DADA2

A continuación, una vez ya hemos importado los datos en qiime2, vamos a determinar los ASVs por medio de DADA2. Este programa lleva a cabo un análisis de calidad, una eliminación de las lecturas que sean de baja calidad y de las quimeras, etc. Para ello, ejecutamos el siguiente comando en qiime:
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 19 \
  --p-trunc-len-f 240 \
  --p-trim-left-r 20 \
  --p-trunc-len-r 155 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza \
  --p-n-threads 7 \
  --p-n-reads-learn 1921748
```
Donde:
* --i-demultiplexed-seqs paired-end-demux.qza: es el nombre del artefacto obtenido anteriormente, que tiene las secuencias asociadas a las muestras junto con otra información adicional.
* --p-trim-left-f 19: número de nucleótidos eliminados en el extremo 5' en las lecturas *forward*.
* --p-trunc-len-f 240: posición a partir de la cual eliminamos las lecturas *forward*.
* --p-trim-left-r 20: número de nucleótidos eliminados en el extremo 5' en las lecturas *reverse*.
* --p-trunc-len-r 155: posición a partir de la cual eliminamos las lecturas *reverse*.
* --o-representative-sequences rep-seqs.qza: nombre del fichero que tendrá las secuencias de cada ASV.
* --o-table table.qza: nombre del fichero correspondiente a la tabla de abundancias de cada ASV por muestra.
* --o-denoising-stats stats.qza: nombre del fichero que contiene información sobre la evolución del proceso.
* --p-n-threads 7: número de procesos llevados a cabo en paralelo.
* --p-n-reads-learn 1921748: número de secuencias usadas para determinar la tasa de error (en este caso, 1921748 es el 25% del número total de secuencias, 7686994, siendo un porcentaje adecuado).

Tal y como puede apreciarse, los valores de estos parámetros son correspondientes a los usados por el profesor cuando ejecutó este comando; esto es algo que podemos comprobar al abrir, por ejemplo, rep-seqs.qza en view.qiime2.org y pinchar en la pestaña "Provenance":

![image](https://user-images.githubusercontent.com/98259577/166155264-623d0bf1-c40b-43d7-a72b-021e58824f25.png)

Este será el procedimiento que sigamos para determinar qué parámetros usó el profesor al ejecutar aquellos comandos que dependan de parámetros determinados.

### 3. Importamos los datos de DADA2 a Qiime2

Obtenidos los ficheros rep-seqs.qza, table.qza y stats.qza tras la ejecución de DADA2, vamos a importar estos datos a qiime2. Para ello, ejecutaremos los siguientes comandos en qiime2 a partir de los tres ficheros obtenidos tras la ejecución del protocolo del fichero dada.rmd aportado por el profesor, los cuales son rep-seqs.fna, seqtab-nochim y stats.txt:
```
qiime tools import \
      --input-path rep-seqs.fna \
      --type 'FeatureData[Sequence]' \
      --output-path rep-seqs.qza
```
Con el comando anterior, hemos transformado los archivos relativos a las secuencias y las abundancias en artefactos para poder usarlos en qiime2. A continuación, vamos a modificar ligeramente el fichero seqtab-nochim.txt e indexar el fichero resultante (biom-table.txt), por medio de los siguientes comandos:
```
echo -n "#OTU Table" | cat - seqtab-nochim.txt > biom-table.txt
biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5
```
Tras obtener table.biom que contiene la información organizada para qiime2, vamos a importar la tabla de ASVs en qiime2:
```
qiime tools import \
      --input-path table.biom \
      --type 'FeatureTable[Frequency]' \
      --input-format BIOMV210Format \
      --output-path table.qza
```
Y finalmente obtener los ficheros de visualización de table.qza y rep-seqs.qza:
```
qiime feature-table summarize \
      --i-table table.qza \
      --o-visualization table.qzv \
      --m-sample-metadata-file metadata

qiime feature-table tabulate-seqs \
      --i-data rep-seqs.qza \
      --o-visualization rep-seqs.qzv
```
Si echamos un vistazo a table.qzv, vemos que este archivo nos aporta información relativa al número de muestras, el número de features, cómo se distribuyen en función de las muestras, etc. Por su parte, rep-seqs.qzv nos permite conocer los IDs de los features, sus longitudes, sus secuencias, etc. 

### 4. Determinación de las distancias filogenéticas por medio de MAFFT y FastTree

A continuación, con el fin de tener una idea acerca de la similitud/disimilitud existente entre los ASVs, vamos a alinear las secuencias representativas (rep-seqs.qza) usando MAFFT para el alineamiento y FastTree para construir el árbol filogenético. Para ello, ejecutamos el siguiente comando:
```
qiime phylogeny align-to-tree-mafft-fasttree \
                --i-sequences rep-seqs.qza \
                --o-alignment aligned-rep-seqs.qza \
                --o-masked-alignment masked-aligned-rep-seqs.qza \
                --o-tree unrooted-tree.qza \
                --o-rooted-tree rooted-tree.qza
```
Donde:
* --i-sequences rep-seqs.qza: archivo con las secuencias de los ASVs.
* --o-alignment aligned-rep-seqs.qza: fichero con las secuencias alineadas.
* --o-masked-alignment masked-aligned-rep-seqs.qza: fichero de las secuencias enmascaradas.
* --o-tree unrooted-tree.qza: fichero correspondiente al árbol sin enraizar.
* --o-rooted-tree rooted-tree.qza: fichero correspondiente al árbol enraizado.

### 5. Asignación taxonómica

Ya tenemos las lecturas agrupadas en ASVs, así como los árboles enraizado y no enraizado. A continuación, nos interesa conocer la taxonomía de cada ASV identificado, y eso es lo que vamos a hacer a continuación. 

En primer lugar, es necesario tomar una base de datos de secuencias y sus correspondientes asignaciones taxonómicas. En este caso, usaremos la base de datos Silva y el plugin rescript. Para ello, vamos a ejecutar el siguiente comando en qiime2:
```
qiime rescript get-silva-data \
      --p-version '138' \
      --p-target 'SSURef_NR99' \
      --p-include-species-labels \
      --o-silva-sequences silva-138-ssu-nr99-seqs.qza \
      --o-silva-taxonomy silva-138-ssu-nr99-tax.qza
```
Ahora, vamos a eliminar las secuencias que tengan bases ambiguas:
```
qiime rescript cull-seqs \
      --i-sequences silva-138-ssu-nr99-seqs.qza \
      --o-clean-sequences silva-138-ssu-nr99-seqs-cleaned.qza
```
Tras esto, vamos a eliminar también las secuencias pequeñas o inútiles:
```
qiime rescript filter-seqs-length-by-taxon \
      --i-sequences silva-138-ssu-nr99-seqs-cleaned.qza \
      --i-taxonomy silva-138-ssu-nr99-tax.qza \
      --p-labels Archaea Bacteria Eukaryota \
      --p-min-lens 900 1200 1400 \
      --o-filtered-seqs silva-138-ssu-nr99-seqs-filt.qza \
      --o-discarded-seqs silva-138-ssu-nr99-seqs-discard.qza
```
Finalmente, eliminamos también aquellas secuencias que sean idénticas pero con distinta taxonomía, con el fin de evitar asignaciones ambiguas:
```
qiime rescript dereplicate \
      --i-sequences silva-138-ssu-nr99-seqs-filt.qza  \
      --i-taxa silva-138-ssu-nr99-tax.qza \
      --p-rank-handles 'silva' \
      --p-mode 'uniq' \
      --o-dereplicated-sequences silva-138-ssu-nr99-seqs-derep-uniq.qza \
      --o-dereplicated-taxa silva-138-ssu-nr99-tax-derep-uniq.qza
```
Los ficheros silva-138-ssu-nr99-seqs-derep-uniq.qza y silva-138-ssu-nr99-tax-derep-uniq.qza son equivalentes a los ficheros 85_otus.qza y ref-taxonomy.qza aportados por el profesor.

___

Tras estos pasos, vamos a importar la base de datos Greegenes 85% (las secuencias están agrupadas en base a una homología del 85%), con el fin de acelerar el proceso de asignación taxonómica, si bien lo más apropiado habría sido usar la base de datos 99% SILVA.

Para ello, importamos las secuencias en el artefacto correspondiente, así como los datos relativos a la asignación taxonómica:
```
qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path 85_otus.fasta \
      --output-path 85_otus.qza
      
qiime tools import \
     --type 'FeatureData[Taxonomy]' \
     --input-format HeaderlessTSVTaxonomyFormat \
     --input-path 85_otu_taxonomy.txt \
     --output-path ref-taxonomy.qza
```
Con el fin de evitar errores en la asignación, vamos a entrenar el clasificador para que contenga solo la parte de cada secuencia que hayamos amplificado. Si consultamos el fichero v4 primers.txt, podremos ver que los primers usados son 515F (5′-GTGYCAGCMGCCGCGGTAA-3′) y 806R (5′-GGACTACNVGGGTWTCTAAT-3′), que amplifican la región V4. Por tanto, extraeremos del artefacto correspondiente la región de interés por medio del siguiente comando:
```
qiime feature-classifier extract-reads \
      --i-sequences 85_otus.qza \
      --p-f-primer GTGYCAGCMGCCGCGGTAA \
      --p-r-primer GGACTACHVGGGTWTCTAAT \
      --p-min-length 100 \
      --p-max-length 400 \
      --o-reads ref-seqs.qza
```
Como vemos, el comando anterior ha tomado el artefacto correspondiente a las secuencias, así como los primers indicados en el párrafo anterior, proporcionando el fichero ref-seqs.qza con las secuencias de 85_otus.qza que pueden ser amplificadas usandos dichos primers permitiendo un rango de longitud entre 100 y 400 nucleótidos, algo que hemos determinado al abrir el fichero ref-seqs.qza en view.qiime2.org y pinchar en "Provenance".

Tras esto, ya podemos generar el clasificador final, en el que asociamos las secuencias con su asignación taxonómica correspondiente, por medio del siguiente comando de qiime2:
```
qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier.qza
```
Una vez ya hemos generado nuestro clasificador, vamos a llevar a cabo la asignación taxonómica de las secuencias representativas de cada ASV, por medio del siguiente comando, que generará un fichero llamado "classification.qza" con dicha información en una carpeta llamada "taxa":
```
qiime feature-classifier classify-sklearn --i-reads rep-seqs.qza \
                                          --i-classifier classifier.qza \
                                          --p-n-jobs 2 \
                                          --output-dir taxa
```
Obtenemos así una carpeta "taxa" en la que tenemos las asignaciones taxonómicas de los ASVs y las representaciones gráficas. Con --p-n-jobs 2, indicamos al sistema que analice 2 ASVs al mismo tiempo.

Tras esto, puesto que estamos interesados en conocer la abundancia de cada ASV en cada muestra y representar las mismas, vamos a generar el fichero taxa_barplot.qzv que nos permitirá visualizar esto, y que ha sido proporcionado por el profesor:
```
qiime taxa barplot --i-table table.qza \
                   --i-taxonomy taxa/classification.qza \
                   --m-metadata-file metadata \
                   --o-visualization taxa/taxa_barplot.qzv
```
Adicionalmente, si bien la gráfica de taxa_barplot.qzv muestra la proporción de ASVs para todas las muestras, también podemos ver esta información para una de las categorías del fichero de metadatos. Concretamente, el profesor agrupó en base a la columna "Day_Temp", por medio del siguiente comando:
```
qiime feature-table group --i-table table.qza \
                          --p-axis sample \
                          --p-mode mean-ceiling \
                          --m-metadata-file metadata \
                          --m-metadata-column Day_Temp \
                          --o-grouped-table table_sample.qza
```
Los resultados pueden ser representados en base a esta variable, por medio del siguiente comando en el que se ha usado otro fichero metadata llamado metadata_sample:
```
qiime taxa barplot --i-table table_sample.qza \
                   --i-taxonomy taxa/classification.qza \
                   --m-metadata-file metadata_sample \
                   --o-visualization taxa/taxa_sample_barplot.qzv
```
Ademas, la tabla anterior la podemos visualizar también en view.qiime2.org por medio de:
```
qiime feature-table summarize \
      --i-table table_sample.qza \
      --o-visualization table_sample.qzv \
      --m-sample-metadata-file metadata_sample
```
### 6. Estudio de la diversidad

