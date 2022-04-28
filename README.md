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
>10000

```
wc -l virome_R2.fastq | awk '{print $1/4}'
```
>10000

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
trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qf1_paired.fq virome_R1_qf1_unpaired.fq virome_R2_qf1_paired.fq virome_R2_qf1_unpaired.fq SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 virome_R1.fastq virome_R2.fastq virome_R1_qf2_paired.fq virome_R1_qf2_unpaired.fq virome_R2_qf2_paired.fq virome_R2_qf2_unpaired.fq SLIDINGWINDOW:4:15 MINLEN:36
```
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

### 2. Ensamblaje *de novo* con *Spades*

spades.py -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o
virome_spades_default

spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o
virome_spades_meta_33 -k33

spades.py --meta -m 3 -t 2 -1 virome_clean.1.fq -2 virome_clean.2.fq -o
virome_spades_meta_21 -k21
