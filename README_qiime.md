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

