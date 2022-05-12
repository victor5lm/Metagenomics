# *QIIME2 homework*
## Víctor Manuel López Molina
## Fecha: 20/05/2022
## Ejercicio

En este ejercicio de la asignatura, dados una serie de ficheros y carpetas proporcionados por el profesor, vamos a redactar un informe en el que se muestren los comandos, con sus correspondientes parámetros, ejecutados por el profesor para obtener dichos archivos. Adicionalmente, junto con este informe, se rellenará el cuestionario disponible en Moodle. 

## INFORME

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
* --p-n-reads-learn 1921748: número de secuencias usadas para determinar la tasa de error (en este caso, 1921748 es el 25% del número total de secuencias, 7686994, siendo éste, por tanto, un porcentaje adecuado).

Tal y como puede apreciarse, los valores de estos parámetros son correspondientes a los usados por el profesor cuando ejecutó este comando; esto es algo que podemos comprobar al abrir, por ejemplo, rep-seqs.qza en view.qiime2.org y pinchar en la pestaña "Provenance":

![image](https://user-images.githubusercontent.com/98259577/166155264-623d0bf1-c40b-43d7-a72b-021e58824f25.png)

Este será el procedimiento que sigamos para determinar qué parámetros usó el profesor al ejecutar aquellos comandos que dependan de parámetros determinados.

Por otro lado, y antes de continuar, conviene prestar atención a los valores indicados por el profesor para dichos parámetros, y analizar por qué han sido considerados éstos. Por un lado, con los parámetros --p-trim-left-f y --p-trim-left-r, hemos indicado a qiime2 que elimine las primeras 19 y 20 bases, respectivamente, de los extremos 5' de las lecturas *forward* y *reverse*. Esto se debe a que, si consultamos la pestaña "Interactive Quality Plot" de paired-end-demux.qzv, podremos ver que estas 19-20 primeras bases de dichos extremos de las lecturas presentan algo menos de calidad que el resto de las posiciones; de ahí que sean eliminadas. 

Por otra parte, con los parámetros --p-trunc-len-f y --p-trunc-len-r, hemos indicado a qiime2 que tome, como posiciones a partir de las cuales cortar las lecturas *forward* y *reverse*, 240 y 155, respectivamente. Esto se debe a que, consultando dicha pestaña de paired-end-demux.qzv, la calidad de las lecturas comienza a decrecer considerablemente a partir de dichas posiciones; de ahí que se tomen estos valores para dichos parámetros, que son responsables de truncar las secuencias *forward* en posición 240 y las *reverse* en posición 155. Todo esto se puede comprobar a continuación:

![image](https://user-images.githubusercontent.com/98259577/166160724-874594af-97f0-49a2-ac37-d9c9c42c0c1d.png)

![image](https://user-images.githubusercontent.com/98259577/166160737-a770eb21-53f2-43c1-8fce-fcff07146b75.png)

### 3. Importamos los datos de DADA2 a Qiime2

Obtenidos los ficheros rep-seqs.qza, table.qza y stats.qza tras la ejecución de DADA2, vamos a importar estos datos a qiime2. Para ello, podemos bien ejecutar los siguientes comandos para visualizar directamente dichos ficheros en view.qiime2.org:
```
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```

Sin embargo, también podríamos ejecutar los siguientes comandos en qiime2 a partir de los tres ficheros obtenidos tras la ejecución del protocolo del fichero dada.rmd aportado por el profesor, los cuales son rep-seqs.fna, seqtab-nochim y stats.txt:
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
Tras obtener table.biom, que contiene la información organizada para qiime2, vamos a importar la tabla de ASVs en qiime2:
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
Si echamos un vistazo a table.qzv, vemos que este archivo nos aporta información relativa al número de muestras, el número de features, cómo se distribuyen en función de las muestras, etc. Por su parte, rep-seqs.qzv nos permite conocer los IDs de los features, sus longitudes, sus secuencias, etc. Respecto a stats.qzv, este fichero aporta información sobre cómo ha ido variando el número de secuencias asociado a cada muestra a lo largo de los distintos pasos de pre-procesado.

A continuación, vamos a construir el árbol filogenético relativo a las secuencias representativas. Cabe puntualizar, por tanto, que no se ha llevado a cabo el paso de eliminación de los *singletons* y de los ASVs de baja frecuencia ya que, por ejemplo, para el fichero table.qza, lo que se observa como último recuadro en la pestaña "Provenance" en view.qiiem2.org es lo siguiente:

![image](https://user-images.githubusercontent.com/98259577/168052259-052a7e4c-cb0f-4af8-9821-fffdfc410611.png)

Como vemos, este fichero ha surgido tras ejecutar dada2 denoise-paired, y no se han realizado pasos posteriores, por lo que no se ha ejecutado qiime feature-table filter-features con este archivo como input. Lo mismo ocurre para rep-seqs.qza, ya que no consta que se haya ejecutado el comando qiime feature-table filter-seqs. Respecto al fichero table_sample.qza, este sí presenta un paso adicional, además de dada2 denoise-paired, en la pestaña "Provenance" de view.qiime2.org, pero este paso se corresponde con el comando qiime feature-table group, algo que se detallará posteriormente en este informe. Por consiguiente, no se han eliminado los *singletons* ni los ASVs de baja frecuencia de los ficheros relativos a la tabla y secuencias representativas originales.

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
Tal y como puede apreciarse, los valores de estos parámetros son correspondientes a los usados por el profesor cuando ejecutó este comando; esto es algo que podemos comprobar al abrir ref-seqs.qza en view.qiime2.org y pinchar en la pestaña "Provenance":

![image](https://user-images.githubusercontent.com/98259577/167139574-bfb8e095-25cc-4a16-8d23-c9c2c30ca736.png)

Como vemos, el comando anterior ha tomado el artefacto correspondiente a las secuencias, así como los primers indicados en el párrafo superior al mismo, proporcionando el fichero ref-seqs.qza con las secuencias de 85_otus.qza que pueden ser amplificadas usandos dichos primers permitiendo un rango de longitud entre 100 y 400 nucleótidos.

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

Finalmente, el profesor también proporcionó archivos correspondientes a estudios de alfa y beta diversidad. Por tanto, vamos a reconstruir estos pasos también en este informe.

En primer lugar, con el fin de determinar si existen diferencias significativas de diversidad entre las muestras del experimento, vamos a prestar atención a dos conceptos: la alfa diversidad (variación de las poblaciones dentro de una muestra concreta) y la beta diversidad (variación de las poblaciones entre distintas muestras).

Puesto que el estudio de la alfa-diversidad es útil para determinar si la profundidad del mismo es suficiente para cubrir todos los taxones de la muestra, una forma de averiguar esto de forma visual es mediante un test de rarefacción, en el que muestras aleatorias son tomadas de un experimento al modificar la profundidad, con el fin de determinar el número de ASVs que aparecen. Las gráficas de las curvas de rarefacción crecen hasta cierto valor de profundidad, de forma que por encima de dicho valor no se accede a información adicional acerca de la biodiversidad. Para obtener estas gráficas, el profesor ejecutó el siguiente comando el qiime2:
```
qiime diversity alpha-rarefaction --i-table table.qza \
                                  --p-max-depth 288000 \
                                  --p-steps 100 \
                                  --i-phylogeny rooted-tree.qza \
                                  --m-metadata-file metadata \
                                  --o-visualization rarefaction_curves.qzv
```
Los valores de estos parámetros son correspondientes a los usados por el profesor cuando ejecutó este comando; esto es algo que podemos comprobar al abrir rarefaction_curves.qzv en view.qiime2.org y pinchar en la pestaña "Provenance":

![image](https://user-images.githubusercontent.com/98259577/167140229-b079d548-6db4-4850-82e2-fe500a449daf.png)

El profesor, en el parámetro --p-max-depth, introdujo el valor 288000, ya que, si consultamos la pestaña "Interactive Sample Detail" de table.qzv, vemos que el número de reads que la muestra más abundante posee es aproximadamente 288000, por lo que este será el máximo valor de profundidad a ser evaluado. Por otro lado, indicando --p-steps 100, el comando ha contado el número de ASVs cada 2880 secuencias.

Tras ejecutar este comando y observar las curvas de rarefacción, también podemos obtener diversos índices de medición de la alfa-diversidad (como Shannon) y de la beta-diversidad (como Unifrac, Jaccard, etc), por medio del siguiente comando:
```
qiime diversity core-metrics-phylogenetic --i-table table.qza \
                                          --i-phylogeny rooted-tree.qza \
                                          --p-sampling-depth 85000 \
                                          --m-metadata-file metadata \
                                          --p-n-jobs-or-threads 2 \
                                          --output-dir diversity
```
Los valores de estos parámetros son correspondientes a los usados por el profesor cuando ejecutó este comando; esto es algo que podemos comprobar al abrir cualquiera de los ficheros de la carpeta "diversity" en view.qiime2.org y pinchar en la pestaña "Provenance":

![image](https://user-images.githubusercontent.com/98259577/167140535-e195248a-f3f6-4709-9f66-e76d911f7ca0.png)

Para la ejecución de este comando, el profesor ha usado --p-sampling-depth 85000. Este parámetro es fundamental para este paso ya que la mayoría de las medidas de diversidad son sensibles a las diferencias existentes en la profundidad de secuenciación de las muestras, por lo que, por medio de este parámetro, todas las muestras acaban teniendo el mismo número de counts, en este caso, 85000 (es una forma de normalización). Si el número total de counts de una muestra es menor a dicho valor, la muestra no se tiene en cuenta para los análisis de diversidad, por lo que lo ideal es tomar un valor que sea lo más alto posible pero que excluya el menor número de muestras posible. Por tanto, si echamos un vistazo de nuevo a la pestaña "Interactive Sample Detail" del fichero table.qzv, veremos que la muestra que tiene menor número de counts tiene 85612 counts; de ahí que el profesor haya escogido 85000 para el parámetro --p-sampling-depth, para así evitar la pérdida de muestras a la hora de determinar la alfa-diversidad por medio de diversas medidas.

Los índices obtenidos, guardados en la carpeta "diversity", pueden ser, como hemos comprobado, visualizados y consultados en view.qiime2.org.

Estos índices de medición de las alfa y beta diversidades también los podemos obtener para la tabla agrupada según la variable "Day_Temp", por medio del siguiente comando:
```
qiime diversity core-metrics-phylogenetic --i-table table_sample.qza \
                                          --i-phylogeny rooted-tree.qza \
                                          --p-sampling-depth 101046 \
                                          --m-metadata-file metadata \
                                          --p-n-jobs-or-threads 2 \
                                          --output-dir diversity_sample
```
Los valores de estos parámetros son correspondientes a los usados por el profesor cuando ejecutó este comando; esto es algo que podemos comprobar al abrir cualquiera de los ficheros de la carpeta "diversity_sample" en view.qiime2.org y pinchar en la pestaña "Provenance":

![image](https://user-images.githubusercontent.com/98259577/167140813-72cd2e15-f362-4cd8-8469-28823eb8088b.png)

Para la ejecución de este comando, el profesor ha usado --p-sampling-depth 101046. Este parámetro es fundamental para este paso ya que la mayoría de las medidas de diversidad son sensibles a las diferencias existentes en la profundidad de secuenciación de las muestras, por lo que, por medio de este parámetro, todos los grupos acaban teniendo el mismo número de counts, en este caso, 101046 (es una forma de normalización). Si el número total de counts de un grupo es menor a dicho valor, el grupo no se tiene en cuenta para los análisis de diversidad, por lo que lo ideal es tomar un valor que sea lo más alto posible pero que excluya el menor número de, en este caso, grupos posible. Por tanto, si echamos un vistazo a la pestaña "Interactive Sample Detail" del fichero table_sample.qzv, veremos que el grupo que tiene menor número de counts tiene 101046 counts; de ahí que el profesor haya escogido 101046 para el parámetro --p-sampling-depth, para así evitar la pérdida de grupos y, por consiguiente, de información, a la hora de determinar las diversidades por medio de diversas medidas.

Los índices obtenidos, guardados en la carpeta diversity_sample, pueden ser, como hemos podido comprobar, visualizados y consultados en view.qiime2.org.

Finalmente, también podemos analizar la composición de las muestras usando PERMANOVA por medio del siguiente comando. Éste evaluará si las distancias entre las muestras dentro de un mismo grupo son más similares entre sí que entre éstas y otras muestras de otros grupos:
```
qiime diversity beta-group-significance --i-distance-matrix diversity/weighted_unifrac_distance_matrix.qza \
                                        --m-metadata-file metadata \
                                        --m-metadata-column Day_Temp \
                                        --o-visualization diversity/weighted_unifrac_condition_significance.qzv \
                                        --p-method permanova \
                                        --p-pairwise
```
Por medio del parámetro --p-pairwise, llevamos a cabo tests *pairwise* que permiten determinar qué pares de grupos específicos difieren de otros, si los hay. Este comando se ha llevado a cabo específicamente para la columna "Day_Temp" de los metadatos para agilizar el proceso. Por otro lado, este comando de qiime2 presenta otro parámetro denominado --p-permutations, correspondiente al número de permutaciones para comuputar el p-valor. Al abrir el fichero weighted_unifrac_condition_significance.qzv en view.qiime2.org y pinchar en la pestaña "Provenance", podemos comprobar que el valor de este parámetro es 999, que es el valor por defecto, de ahí que no venga incluido en el comando expuesto en la celda superior a este párrafo.

### 7. Análisis en R con DESeq2

Finalmente, para terminar la práctica, vamos a echar un vistazo al fichero "DESeq2.html" proporcionado por el profesor y tratar de interpretar la gráfica obtenida (el código no se comentará al ser ya aportado por el profesor). Este archivo se corresponde con el análisis de los resultados, concretamente, en este caso, con el fin de determinar diferencias, dentro de un nivel taxónomico concreto, entre las abundancias de diversos taxones. Específicamente, como se indica en este fichero, se ha ejecutado este programa para conocer qué ASVs o taxones son más abundantes en función del día de toma de la muestra (día 0 desde el inicio del experimento o día 31), siendo dichas muestras de heces y de ratones expuestos a frío. Tal y como muestra el profesor al final del archivo, tras ejecutar todos los comandos expuestos se obtiene la siguiente gráfica:

![image](https://user-images.githubusercontent.com/98259577/167214694-016c9b30-7fe9-4f57-b97e-4908f58e4c5d.png)

La forma de interpretar esta representación gráfica depende de lo que hayamos indicado en este comando:
```
rescold <- results(colddds, cooksCutoff = FALSE, contrast=c("Collection_Day", "Day_31", "Day_0"))
```
Al haber indicado como numerador la variable "Day_31" y como denominador la variable "Day_0", las familias taxonómicas que se sitúan por encima de 1 (línea naranja) son más abundantes en muestras de heces, recogidas en frío, tomadas en el día 31 del experimento, mientras que las familias taxonómicas situadas por debajo de -1 (línea rosa) son más abundantes en muestras de heces, recogidas en frío, tomadas en el mismo día del inicio del experimento (día 0). Estos aspectos parecen indicar que, efectivamente, **la exposición de los ratones al frío afecta a la microbiota a medida que pasan los días** ya que, por ejemplo, para la familia *Lachnospiraceae*, del phylum *Firmicutes*, ésta es mucho más abundante en las heces recogidas en frío tras 31 días de experimento en comparación con las heces recogidas en frío el mismo día del experimento. Lo contrario ocurre, por ejemplo, con la familia *S24-7* del phylum *Bacteroidetes*, esto es, ésta es más abundante en heces recogidas en frío el mismo día del experimento en comparación con las heces recogidas en frío tras 31 días de experimento (lo cual equivaldría a decir que dichas heces no difieren respecto a las tomadas a temperatura ambiente ya que, al haber sido recogidas el mismo día del inicio del experimento, no pueden haberse visto afectadas por el frío en tan poco tiempo).

### 8. Conclusión final

Por medio de esta práctica, hemos podido llevar a cabo el análisis de datos metagenómicos por medio de diversas herramientas, sobre todo QIIME2, obteniendo información altamente relevante, como el número de muestras, el número de counts por muestra, el número de ASVs por muestra, las alfa y beta diversidades, las taxonomías de los ASVs identificados, etc. Además, con DESeq2, hemos podido comprobar cómo la temperatura es capaz de modular la microbiota intestinal, siendo ésta por tanto una más de las múltiples variables externas que influyen en la composición de la misma. Esto demuestra, en última instancia, la enorme utilidad de herramientas computacionales como éstas a la hora de obtener información acerca de la composición microbiana de muestras metagenómicas y su variación dadas unas condiciones externas.
