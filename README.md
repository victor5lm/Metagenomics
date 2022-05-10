# *De novo assembly homework*
## Víctor Manuel López Molina
## Fecha: 20/05/2022
## Ejercicio

Repeat all the steps with a viral metagenome from a human saliva sample (Virome.zip in Unit 3 of Moodle). Compare different de novo assemblies options (try _--meta--) or different kmer values. You must perform at least 3 different assemblies. Write a brief summary describing the bioinformatic pipeline you have followed (trimming, decontamination, improve in quality, number of reads remove in each step, etc.). Compare different de novo assemblies with QUAST and choose the best based on the obtained metrics (smaller number of contigs, higher N50, smaller L50, longest total assembly length, etc.).

>NOTE: In the quality filtering step, modify the MINLEN argument considering the original read length. Consider that reads with a minimum of 50% of the average original size are ok for subsequent analyses.

>NOTE: Importantly, you do not have a reference genome for a metagenome.

### 1. Pre-procesado
#### 1.1. Descargamos el *dataset* de interés y lo descomprimimos

```
mkdir -p ~/Documents/tarea_de_novo
mv ~/Downloads/virome.zip ~/Documents/tarea_de_novo
cd ~/Documents/tarea_de_novo
unzip virome.zip
```
Una vez descomprimido el fichero *virome.zip*, vemos que el contenido del mismo es:
```
ls
```
>virome_R1.fastq virome_R2.fastq

#### 1.2 Comprobamos cuál es la integridad de los archivos
A continuación, vamos a comprobar los *hashes* de los dos ficheros expuestos anteriormente:
```
md5sum virome_R1.fastq
```
>cd891dfc865f01c8f3923edd55dacde5  virome_R1.fastq

```
md5sum virome_R2.fastq
```
>28f0d0ace9fb45af8b2595cbc78fa2bd virome_R2.fastq

Conocer la estructura de estos *hashes* es relevante pues éstos no cambian a no ser que el contenido de los ficheros haya sido modificado. Para determinar si existe algo inusual en el fichero, basta con fijarnos en los últimos 5-6 dígitos del *hash* proporcionado por md5.

#### 1.3 Contamos el número de lecturas

En este paso, vamos a determinar si ambos archivos poseen el mismo número de lecturas, por medio de los siguientes comandos:

```
wc -l virome_R1.fastq | awk '{print $1/4}'
```
>100000

```
wc -l virome_R2.fastq | awk '{print $1/4}'
```
>100000

Efectivamente, ambos ficheros poseen el mismo número de *reads*.

#### 1.4 Comprobamos la calidad de las secuencias con la herramienta fastqc
```
mkdir virome_quality
fastqc virome_R1.fastq -o virome_quality/
fastqc virome_R2.fastq -o virome_quality/
```
Si abrimos los ficheros .html obtenidos, obtenemos las siguientes gráficas de calidad, además de numerosa información adicional (%G+C, distribución de la longitud de las lecturas, etc):

**R1 quality plot**

![image](https://user-images.githubusercontent.com/98259577/165406826-d21bf882-c2b1-4d3f-b714-bc3d462cb807.png)

**R2 quality plot**

![image](https://user-images.githubusercontent.com/98259577/165406753-3aa52c07-0258-4f80-867b-fff206035a15.png)

#### 1.5 Eliminación de las lecturas cortas y los extremos de baja calidad (*Trimmomatic*)
Siguiendo el consejo del profesor, vamos a determinar cuál es la media de longitud de las lecturas, para así seleccionar como valor del parámetro MINLEN el 50% de dicho valor.
```
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' virome_R1.fastq
```
>298.939
```
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' virome_R2.fastq
```
>299.249

Como la media de la longitud de las lecturas es 299, el 50% de ésta es 149, por lo que tomaremos MINLEN=149 en *Trimmomatic*. Dejaremos el parámetro SLIDINGWINDOW dentro del rango 4:20. Además, tal y como indica el enunciado del ejercicio, aparte de ejecutar este programa con dichos valores de ambos parámetros, vamos a ejecutarlo también con otros valores. Por ejemplo, vamos a hacerlo, en segundo lugar, con SLIDINGWINDOW:4:20 y MINLEN:70, y luego con SLIDINGWINDOW:4:15 y MINLEN:36.

```
trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qf_paired.fq virome_R1_qf_unpaired.fq virome_R2_qf_paired.fq virome_R2_qf_unpaired.fq SLIDINGWINDOW:4:20 MINLEN:149
```
>Multiple cores found: Using 4 threads
>
>Input Read Pairs: 100000 Both Surviving: 56547 (56.55%) Forward Only Surviving: 26011 (26.01%) Reverse Only Surviving: 1968 (1.97%) Dropped: 15474 (15.47%)
>
>TrimmomaticPE: Completed successfully

```
trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qf1_paired.fq virome_R1_qf1_unpaired.fq virome_R2_qf1_paired.fq virome_R2_qf1_unpaired.fq SLIDINGWINDOW:4:20 MINLEN:70
```
>Multiple cores found: Using 4 threads
>
>Input Read Pairs: 100000 Both Surviving: 80850 (80.85%) Forward Only Surviving: 14158 (14.16%) Reverse Only Surviving: 1026 (1.03%) Dropped: 3966 (3.97%)
>
>TrimmomaticPE: Completed successfully

```
trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qf2_paired.fq virome_R1_qf2_unpaired.fq virome_R2_qf2_paired.fq virome_R2_qf2_unpaired.fq SLIDINGWINDOW:4:15 MINLEN:36
```
>Multiple cores found: Using 4 threads
>
>Input Read Pairs: 100000 Both Surviving: 93796 (93.80%) Forward Only Surviving: 5870 (5.87%) Reverse Only Surviving: 143 (0.14%) Dropped: 191 (0.19%)
>
>TrimmomaticPE: Completed successfully

A continuación, vamos a ejecutar de nuevo *fastqc* para analizar la calidad de las lecturas, para los tres casos, tras la eliminación de los extremos de baja calidad y las lecturas cortas:
```
mkdir virome_QF_quality
fastqc virome_R1_qf_paired.fq -o virome_QF_quality
fastqc virome_R2_qf_paired.fq -o virome_QF_quality
fastqc virome_R1_qf1_paired.fq -o virome_QF_quality
fastqc virome_R2_qf1_paired.fq -o virome_QF_quality
fastqc virome_R1_qf2_paired.fq -o virome_QF_quality
fastqc virome_R2_qf2_paired.fq -o virome_QF_quality
```
Tras la ejecución de estos comandos, obtenemos las siguientes gráficas de calidad:

**virome_R1_qf_paired.fq**

![image](https://user-images.githubusercontent.com/98259577/165759339-d00f26cb-5c5e-4e6a-93f0-c364a21c59e2.png)

**virome_R2_qf_paired.fq**

![image](https://user-images.githubusercontent.com/98259577/165759513-bf9e0e4a-f59e-4664-acd3-f083e84c0a38.png)

**virome_R1_qf1_paired.fq**

![image](https://user-images.githubusercontent.com/98259577/165761454-8ab67688-f96e-4ca1-9c80-88e42b724004.png)

**virome_R2_qf1_paired.fq**

![image](https://user-images.githubusercontent.com/98259577/165761550-0d5d20a1-2b90-4273-898c-436b1dbccccf.png)

**virome_R1_qf2_paired.fq**

![image](https://user-images.githubusercontent.com/98259577/165761658-bec14acb-7de0-4386-ae88-59ec819b6084.png)

**virome_R2_qf2_paired.fq**

![image](https://user-images.githubusercontent.com/98259577/165761720-ef4e39e4-c46e-4bf7-b1a5-79d3a34e7849.png)

Tal y como se puede apreciar, a medida que disminuye el valor del parámetro MINLEN (longitud mínima a tener las lecturas para no ser desechadas), mayor es el número de lecturas que se conservan tras el filtrado de calidad.

#### 1.6. Eliminación de las lecturas que alinean con el genoma humano o el de phiX174 con *Bowtie2*

A continuación, para eliminar estas lecturas correspondientes a contaminaciones, vamos a usar Bowtie2. Para ello, vamos a construir primero un índice que contenga las secuencias de referencia a usar para el alineamiento, en este caso procedentes del ser humano y de phiX174. En primer lugar, ejectuamos estos comandos:
```
awk '{if(/^>/) printf("\n%s\n",$0); else printf("%s",$0)} END {printf("\n")}' < GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna > HumanGenome_ONELINE.fasta
awk '{if(/^>/) printf("\n%s\n",$0); else printf("%s",$0)} END {printf("\n")}' < phiX174.fasta > phiX174_genome_ONELINE.fasta
cat HumanGenome_ONELINE.fasta phiX174_genome_ONELINE.fasta > Human-phiX174.fasta
```
Tras esto, generamos el índice para Bowtie2:
```
bowtie2-build Human-phiX174.fasta human-phix174
```
Finalmente, ejecutamos Bowtie2, usando este índice recién generado, para eliminar las contaminaciones en todos los ficheros .fq generados anteriormente tras la ejecución de *Trimmomatic* con distintos valores de SLIDINGWINDOW y MINLEN:
```
bowtie2 -x human-phix174 -1 virome_R1_qf_paired.fq -2 virome_R2_qf_paired.fq --un-conc virome_clean.fq -S tmp.sam
bowtie2 -x human-phix174 -1 virome_R1_qf1_paired.fq -2 virome_R2_qf1_paired.fq --un-conc virome_clean_1.fq -S tmp.sam
bowtie2 -x human-phix174 -1 virome_R1_qf2_paired.fq -2 virome_R2_qf2_paired.fq --un-conc virome_clean_2.fq -S tmp.sam
```
A continuación, contamos el número de lecturas restantes para cada tipo de fichero obtenido en función de dichos parámetros en *Trimmomatic*:
```
wc -l *clean.* | awk '{print $1/4}'
```
>56474 
>56474 
>112948

```
wc -l *clean_1* | awk '{print $1/4}'
```
>80761
>80761
>161522

```
wc -l *clean_2* | awk '{print $1/4}'
```
>93700
>93700
>187400

Tal y como podemos apreciar, la proporción de contaminaciones existentes entre las lecturas es muy pequeña, ya que el número de lecturas restante, tanto en las *forward* como en las *reverse*, tras la eliminación de las mismas, es muy similar al existente tras el filtrado de calidad realizado anteriormente. Esto puede observarse para los tres filtrados de calidad realizados con *Trimmomatic*, así como los correspondientes alineamientos con los genomas humano y de phiX174.

### 2. Ensamblaje *de novo* con *Spades* y comparación de los distintos ensamblajes con *QUAST*

A continuación, vamos a utilizar el programa *Spades* para el ensamblaje *de novo*. Para ello, ejecutaremos estos comandos, usando asimismo los parámetros -t (para indicar el número de *threads* a usar por el ensamblador), -m (para indicar la memoria que el programa puede usar) y -k (para indicar la longitud de los kmeros a usar durante el ensamblado), así como el *flag* --meta. Llevaremos a cabo, para los ficheros virome_clean.1.fq y virome_clean.2.fq, 3 ensamblajes distintos:

**Ficheros virome_clean.1.fq y virome_clean.2.fq**
```
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_meta_77 -k77
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_meta_99 -k99
spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o virome_spades_meta_127 -k127
```
Tras esto, vamos a pasar todos los contigs/scaffolds a la misma carpeta, procedentes de los 3 ensamblajes, para posteriormente cargarlos en QUAST.
```
mkdir quast
ln -rs ./virome_spades_meta_77/contigs.fasta ./quast/contigs_meta_77.fasta
ln -rs ./virome_spades_meta_77/scaffolds.fasta ./quast/scaffolds_meta_77.fasta
ln -rs ./virome_spades_meta_99/contigs.fasta ./quast/contigs_meta_99.fasta
ln -rs ./virome_spades_meta_99/scaffolds.fasta ./quast/scaffolds_meta_99.fasta
ln -rs ./virome_spades_meta_127/contigs.fasta ./quast/contigs_meta_127.fasta
ln -rs ./virome_spades_meta_127/scaffolds.fasta ./quast/scaffolds_meta_127.fasta
```
Cargamos estos ficheros en QUAST, obteniendo el siguiente report:

**LINK:** http://cab.cc.spbu.ru/quast/reports/29_Apr_2022_16:46:36_462922_1st-trimmomatic-minlen-149-/report.html

![image](https://user-images.githubusercontent.com/98259577/166062721-a6297d1b-5eb3-4203-81f5-5fb5e3f49f14.png)

![image](https://user-images.githubusercontent.com/98259577/166062786-a194788d-af56-423b-b217-8cb47f7149e5.png)

Analizando la tabla aportada por el report de QUAST, podemos determinar cuál de los 3 ensamblajes es mejor para estos ficheros, obtenidos tras ejecutar *Trimmomatic* con los parámetros SLIDINGWINDOW:4:20 y MINLEN:149. Si nos fijamos en el número de contigs obtenidos para k=77, k=99 y k=127 (columnas 2, 3 y 4), podemos apreciar que el ensamblaje con el menor número de contigs es el obtenido con k=127. Además, el mayor contig obtenido también es con k=127. Asimismo, N50 y N75 son mayores para k=127, así como L50 y L75 son menores para dicho valor de k (N75 y L75 son relativos al 75% del ensamblado). Todos estos hechos indican que, para este caso, el mejor ensamblaje es con un tamaño de kmero igual a 127. Esto se fundamenta en que, a menor número de contigs y mayor tamaño del más grande de ellos, más sencillo resulta ensamblar el genoma del virus que estamos analizando en cuestión. Por otro lado, hemos de recordar la definición de N50: longitud para la cual los contigs de dicha longitud o mayor cubren el 50% del ensamblado. Naturalmente, cuanto mayor sea N50, mejor será el ensamblado, porque se necesitarán menos contigs para reconstruir el 50% del genoma. Por su parte, siendo L50 el número de contigs mínimo necesario para cubrir dicho 50%, cuanto menor sea éste, mejor. Es en base a esto por lo que el ensamblado con k=127 resulta ser más apropiado. Esto también puede apreciarse en las siguientes gráficas de QUAST:

![image](https://user-images.githubusercontent.com/98259577/166112500-713f1565-e486-4a71-b070-ebafc71c7dab.png)

![image](https://user-images.githubusercontent.com/98259577/166112514-0833f767-3fd1-496f-af85-9291a04051f4.png)

Analizando estas gráficas, podremos ver que, para la gráfica Nx, N50 es mayor para k=127. Para la otra gráfica, podemos apreciar que, para k=127, el 1º contig es el mayor y el número de contigs es menor para este valor de k.

A continuación, vamos a hacer analizar los resultados para el resto de archivos obtenidos para las otras dos ejecuciones de *Trimmomatic* con distintos valores de los parámetros SLIDINGWINDOW y MINLEN:

**Ficheros virome_clean_1.1.fq y virome_clean_1.2.fq**
```
spades.py --meta -m 3 -t 2 -1 virome_clean_1.1.fq -2 virome_clean_1.2.fq -o virome_spades_meta_77_1 -k77
spades.py --meta -m 3 -t 2 -1 virome_clean_1.1.fq -2 virome_clean_1.2.fq -o virome_spades_meta_99_1 -k99
spades.py --meta -m 3 -t 2 -1 virome_clean_1.1.fq -2 virome_clean_1.2.fq -o virome_spades_meta_127_1 -k127
```
Tras esto, vamos a pasar todos los contigs/scaffolds a la misma carpeta, procedentes de los 3 ensamblajes, para posteriormente cargarlos en QUAST.
```
mkdir quast_1
ln -rs ./virome_spades_meta_77_1/contigs.fasta ./quast_1/contigs_meta_77_1.fasta
ln -rs ./virome_spades_meta_77_1/scaffolds.fasta ./quast_1/scaffolds_meta_77_1.fasta
ln -rs ./virome_spades_meta_99_1/contigs.fasta ./quast_1/contigs_meta_99_1.fasta
ln -rs ./virome_spades_meta_99_1/scaffolds.fasta ./quast_1/scaffolds_meta_99_1.fasta
ln -rs ./virome_spades_meta_127_1/contigs.fasta ./quast_1/contigs_meta_127_1.fasta
ln -rs ./virome_spades_meta_127_1/scaffolds.fasta ./quast_1/scaffolds_meta_127_1.fasta
```
Cargamos estos ficheros en QUAST, obteniendo el siguiente report:

**LINK:** http://cab.cc.spbu.ru/quast/reports/29_Apr_2022_16:47:30_283421_2nd-trimmomatic-minlen-70-/report.html

![image](https://user-images.githubusercontent.com/98259577/166063696-8dc2f0aa-c3f1-4111-9074-3bea65292a01.png)

![image](https://user-images.githubusercontent.com/98259577/166064254-390ae184-6568-4f1a-a807-3f0dce3cb77f.png)

Analizando la tabla aportada por el report de QUAST, podemos determinar cuál de los 3 ensamblajes es mejor para estos ficheros, obtenidos tras ejecutar *Trimmomatic* con los parámetros SLIDINGWINDOW:4:20 y MINLEN:70. Si nos fijamos en el número de contigs obtenidos para k=77, k=99 y k=127 (columnas 2, 3 y 4), podemos apreciar que el ensamblaje con el menor número de contigs es el obtenido con k=127. Además, el mayor contig obtenido también es con k=127. Asimismo, N50 y N75 son mayores para k=127, así como L50 y L75 son menores para dicho valor de k. Todos estos hechos indican que, para este caso, el mejor ensamblaje es con un tamaño de kmero igual a 127. Esto se fundamenta en que, a menor número de contigs y mayor tamaño del más grande de ellos, más sencillo resulta ensamblar el genoma del virus que estamos analizando en cuestión. Por otro lado, hemos de recordar la definición de N50: longitud para la cual los contigs de dicha longitud o mayor cubren el 50% del ensamblado. Naturalmente, cuanto mayor sea N50, mejor será el ensamblado, porque se necesitarán menos contigs para reconstruir el 50% del genoma. Por su parte, siendo L50 el número de contigs mínimo necesario para cubrir dicho 50%, cuanto menor sea éste, mejor. Es en base a esto por lo que el ensamblado con k=127 resulta ser más apropiado. Esto también puede apreciarse en las siguientes gráficas de QUAST:

![image](https://user-images.githubusercontent.com/98259577/166113323-b6a54b15-1f40-4a8b-be7b-7859b4353a33.png)

![image](https://user-images.githubusercontent.com/98259577/166113337-1c220330-018d-4b87-be3d-6f5a01461f0a.png)

Analizando estas gráficas, podremos ver que, para la gráfica Nx, N50 es mayor para k=127. Para la otra gráfica, podemos apreciar que, para k=127, el 1º contig es el mayor y el número de contigs es menor para este valor de k.

Finalmente, vamos a hacer analizar los resultados para los archivos obtenidos para la última ejecución de *Trimmomatic*:

**Ficheros virome_clean_2.1.fq y virome_clean_2.2.fq**
```
spades.py --meta -m 3 -t 2 -1 virome_clean_2.1.fq -2 virome_clean_2.2.fq -o virome_spades_meta_77_2 -k77
spades.py --meta -m 3 -t 2 -1 virome_clean_2.1.fq -2 virome_clean_2.2.fq -o virome_spades_meta_99_2 -k99
spades.py --meta -m 3 -t 2 -1 virome_clean_2.1.fq -2 virome_clean_2.2.fq -o virome_spades_meta_127_2 -k127
```
Tras esto, vamos a pasar todos los contigs/scaffolds a la misma carpeta, procedentes de los 3 ensamblajes, para posteriormente cargarlos en QUAST.
```
mkdir quast_2
ln -rs ./virome_spades_meta_77_2/contigs.fasta ./quast_2/contigs_meta_77_2.fasta
ln -rs ./virome_spades_meta_77_2/scaffolds.fasta ./quast_2/scaffolds_meta_77_2.fasta
ln -rs ./virome_spades_meta_99_2/contigs.fasta ./quast_2/contigs_meta_99_2.fasta
ln -rs ./virome_spades_meta_99_2/scaffolds.fasta ./quast_2/scaffolds_meta_99_2.fasta
ln -rs ./virome_spades_meta_127_2/contigs.fasta ./quast_2/contigs_meta_127_2.fasta
ln -rs ./virome_spades_meta_127_2/scaffolds.fasta ./quast_2/scaffolds_meta_127_2.fasta
```
Cargamos estos ficheros en QUAST, obteniendo el siguiente report:

**LINK:** http://cab.cc.spbu.ru/quast/reports/29_Apr_2022_16:48:26_624564_3rd-trimmomatic-minlen-36-/report.html

![image](https://user-images.githubusercontent.com/98259577/166067785-ddb7789c-2868-4698-87c6-db1901d62024.png)

![image](https://user-images.githubusercontent.com/98259577/166067850-b467fddb-8a77-4ec5-867e-21046382ffe4.png)

Analizando la tabla aportada por el report de QUAST, podemos determinar cuál de los 3 ensamblajes es mejor para estos ficheros, obtenidos tras ejecutar *Trimmomatic* con los parámetros SLIDINGWINDOW:4:20 y MINLEN:70. Si nos fijamos en el número de contigs obtenidos para k=77, k=99 y k=127 (columnas 2, 3 y 4), podemos apreciar que el ensamblaje con el menor número de contigs es el obtenido con k=127. Además, el mayor contig obtenido también es con k=127. Asimismo, N50 y N75 son mayores para k=127, así como L50 y L75 son menores para dicho valor de k. Todos estos hechos indican que, para este caso, el mejor ensamblaje es con un tamaño de kmero igual a 127. Esto se fundamenta en que, a menor número de contigs y mayor tamaño del más grande de ellos, más sencillo resulta ensamblar el genoma del virus que estamos analizando en cuestión. Por otro lado, hemos de recordar la definición de N50: longitud para la cual los contigs de dicha longitud o mayor cubren el 50% del ensamblado. Naturalmente, cuanto mayor sea N50, mejor será el ensamblado, porque se necesitarán menos contigs para reconstruir el 50% del genoma. Por su parte, siendo L50 el número de contigs mínimo necesario para cubrir dicho 50%, cuanto menor sea éste, mejor. Es en base a esto por lo que el ensamblado con k=127 resulta ser más apropiado. Esto también puede apreciarse en las siguientes gráficas de QUAST:

![image](https://user-images.githubusercontent.com/98259577/166113448-8f995bf8-589d-447e-8a67-6d6ea6b7102d.png)

![image](https://user-images.githubusercontent.com/98259577/166113454-7083f334-2cd3-4002-a942-ced2026d1140.png)

Analizando estas gráficas, podremos ver que, para la gráfica Nx, N50 es mayor para k=127. Para la otra gráfica, podemos apreciar que, para k=127, el 1º contig es el mayor y el número de contigs es menor para este valor de k.

### 3. Conclusión final

A modo de resumen, hemos partido de datos metagenómicos de un virus con el fin de realizar un ensamblado *de novo*. En concreto, las lecturas obtenidas tras la secuenciación se han pre-procesado de tres formas distintas en base a los distintos parámetros del programa *Trimmomatic*. Para cada uno de los resultados obtenidos, hemos realizado el ensamblado del genoma, a su vez, de tres formas diferentes, por medio del programa SPAdes, en base a la longitud del kmero escogida. 

En los 3 casos, hemos podido comprobar que, con **k=127**, el ensamblado es mucho **mejor** que con otros valores de k, comparando parámetros como N50, L50, N75, L75, el tamaño del mayor contig, el número de contigs, etc. Si comparamos, asimismo, estos tres mejores ensamblados entre sí, podemos comprobar que, habiendo realizado el control de calidad con *Trimmomatic* con los parámetros SLIDINGWINDOW:4:20 y MINLEN:149, es decir, indicando para MINLEN el 50% del tamaño medio de las lecturas originales, **este ensamblado es el** **mejor**.

Concretamente, podemos apreciar que, si bien el N50 no es el mayor de todos para este ensamblado, L50 y L75 sí que son los mejores (es decir, los más bajos en comparación). Aunque buscamos que N50 sea lo más alto posible, al fin y al cabo también buscamos que L50 sea lo más bajo posible, ya que es el mínimo número de contigs que cubren la mitad del ensamblado. Por eso, aunque N50 no sea el más alto para este ensamblado en comparación con los otros dos, L50 y L75 (es el mínimo número de contigs que cubren el 75% del ensamblado) sí que son los más bajos, y en definitiva esto resulta ventajoso al reconstruir el genoma de este virus.

**En conclusión, el mejor ensamblado ha demostrado ser el obtenido tras ejecutar *Trimmomatic* con SLIDINGWINDOW:4:20 y MINLEN:149, además de SPAdes con el *flag* --meta y k=127.**

Como hemos podido comprobar, el empleo de todas estas herramientas de pre-procesado de las lecturas, ensamblado (con SPAdes) y comparación de distintos ensamblados (con QUAST) han resultado ser altamente útiles para la reconstrucción de un metagenoma, en este caso de un virus, al poder guiarnos en torno a qué parámetros son más óptimos a la hora de realizar dichos ensamblados con SPAdes.
