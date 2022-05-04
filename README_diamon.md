# Asignación taxonómica con DIAMOND+MEGAN
## Víctor Manuel López Molina
## Fecha: 20/05/2022
## Ejercicio

Follow the workflow of this tutorial for the taxonomic binning of the virome reads and contigs used as homework in the unit_3.

>Reads from unit_3 homework are paired_end reads. To perform this task you can do it in several ways; joining all decontaminated and QF reads and then run DIAMOND or run DIAMOND using both files and then join DIAMOND outputs, etc..

Write a brief summary describing the bioinformatic pipeline you have followed (trimming, decontamination, improve in quality, number of reads remove in each step, etc.) and the most relevant results (two figures maximun) with the taxonomy assessment of the virome at family level.

**Optional**: run the same analysis but using the assembled contigs and try to compare the results using MEGAN (File/Compare and load both rma6 files created while making the individual analysis).

### 1. Pre-procesado
En primer lugar, partimos de los ficheros virome_R1.fastq y virome_R2.fastq que habíamos obtenido en la práctica anterior a partir del fichero virome.zip. Al tratarse de lecturas paired-end, y al deber trabajar con un único archivo, vamos a fusionar el contenido del ambos archivos en uno solo (llamado virome.fq), por medio del siguiente comando:
```
cat virome_R1.fastq virome_R2.fastq > virome.fq
```
Tras esto, vamos a llevar a cabo el filtrado de las lecturas en base a su calidad por medio de Trimmomatic, utilizando los parámetros que demostraron ser más adecuados en la práctica anterior:
```
trimmomatic SE -phred33 virome.fq virome_qf.fq SLIDINGWINDOW:4:20 MINLEN:149
```
>Input Reads: 200000
>
>Surviving: 141073 (70.54%) 
>
>Dropped: 58927 (29.46%)

Finalmente, tras el filtrado de calidad, vamos a eliminar posibles contaminaciones del fichero virome.fq por medio de bowtie2:
```
for f in ./Index/*bt2; do ln -s $f .; done
bowtie2 -x human-phix174 -q virome_qf.fq --un virome_qf_clean.fq -S tmp.sam
```
>141073 reads; of these:
>
>141073 (100.00%) were unpaired; of these:
>
>140903 (99.88%) aligned 0 times
>
>162 (0.11%) aligned exactly 1 time
>
>8 (0.01%) aligned >1 times
>
>0.12% overall alignment rate
---
### 2. Alineamiento de las lecturas de alta calidad y descontaminadas con una base de datos de proteínas virales

#### 2.1: Preparación de la base de datos de proteínas virales

En primer lugar, vamos a descargar los archivos relativos a la base de datos a partir del siguiente enlace: https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/

Dentro de los archivos listados, descargamos los siguientes:
* viral.1.protein.faa.gz
* viral.2.protein.faa.gz
* viral.3.protein.faa.gz
* viral.4.protein.faa.gz

Tras haber descargado todos los archivos, los descomprimimos y los unimos en un único archivo .faa:
```
gunzip viral.*.protein.faa.gz
cat viral.*.protein.faa > viral.protein.faa
grep -c ">" *faa
```
#### 2.2: Creamos la base de datos a usar por DIAMOND

```
diamond makedb --in viral.protein.faa -d viralproteins
```
#### 2.3: Alineamiento por medio de blastx
```
diamond blastx -d viralproteins.dmnd -q virome_qf_clean.fq -o virome_qf_clean_vs_viralprotein.m8
```
>Total time = 16.779s
>
>Reported 347842 pairwise alignments, 347842 HSPs.  
>
>34649 queries aligned.

Tras el alineamiento, podemos echar un vistazo al resultado obtenido:
```
head virome_qf_clean_vs_viralprotein.m8
```
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| M02255:131:000000000-AJC6R:1:1101:8786:11380_1:N:0:AGTCAA	| YP_009216661.1 | 50.8 | 63 | 30 | 1 | 212 | 27 | 104 | 166 | 3.7e-09 | 61.2 |
| M02255:131:000000000-AJC6R:1:1101:8786:11380_1:N:0:AGTCAA | YP_006383483.1 | 50.0 | 48 | 23 | 1 | 167 | 27 | 120 | 167 | 1.2e-04 | 46.2 |
| M02255:131:000000000-AJC6R:1:1101:8786:11380_1:N:0:AGTCAA	| YP_009837323.1 | 47.9 | 48 | 24 | 1 | 167 | 27 | 120 | 167 | 3.6e-04 | 44.7 |
| M02255:131:000000000-AJC6R:1:1101:8786:11380_1:N:0:AGTCAA | YP_009840443.1 | 47.9 | 48 | 24 | 1 | 167 | 27 | 120 | 167 | 3.6e-04 | 44.7 |
| M02255:128:000000000-AG7E5:1:2116:15291:16348_1:N:0:AGTCAA | NP_795683.1 | 78.7 | 94 | 20 | 0 | 284 | 3 | 81 | 174 | 7.3e-37 | 153.7 |
| M02255:128:000000000-AG7E5:1:2116:15291:16348_1:N:0:AGTCAA | YP_001469206.1 | 78.7 | 94 | 20 | 0 | 284 | 3 | 81 | 174 | 7.3e-37 | 153.7 |
| M02255:128:000000000-AG7E5:1:2116:15291:16348_1:N:0:AGTCAA | NP_803283.1 | 73.4 | 94 | 25 | 0 | 284 | 3 | 97 | 190 | 2.6e-34 | 145.2 |
| M02255:128:000000000-AG7E5:1:2116:15291:16348_1:N:0:AGTCAA | YP_001285355.1 | 73.4 | 94 | 25 | 0 | 284 | 3 | 97 | 190 | 2.6e-34 | 145.2 |
| M02255:128:000000000-AG7E5:1:2116:15291:16348_1:N:0:AGTCAA | YP_239566.1 | 73.4 | 94 | 25 | 0 | 284 | 3 | 97 | 190 | 2.6e-34 | 145.2 |
| M02255:128:000000000-AG7E5:1:2116:15291:16348_1:N:0:AGTCAA | YP_239642.1 | 73.4 | 94 | 25 | 0 | 284 | 3 | 97 | 190 | 2.6e-34 | 145.2 |

A continuación, vamos a analizar la taxonomía de estas secuencias por medio de MEGAN6.

---
### 3. Análisis de la taxonomía de las lecturas con MEGAN6

MEGAN6 analiza la taxonomía de las lecturas alineadas con la base de datos de proteínas virales tomada de NCBI y genera un árbol taxonómico.

#### 3.1: Importamos los archivos de DIAMOND en MEGAN6

Importamos los ficheros pertinentes en MEGAN6, concretamente de esta forma, tras haber pinchado en File->Import from BLAST:

![image](https://user-images.githubusercontent.com/98259577/166725000-374b2d61-7ac0-4117-982b-996fbc2f6244.png)

Tal y como puede apreciarse, hemos seleccionado el formato BlastTab y el modo BlastX. Posteriormente, hacemos click en "Next" y subimos también el archivo prot_acc2tax-Jul2019X1.abin.zip (descomprimido previamente), necesario para aportar información taxonómica a MEGAN6, pulsando en el símbolo de la carpeta que podemos ver al lado de "Load Accession mapping file", el cual se puede apreciar en la siguiente imagen:

![image](https://user-images.githubusercontent.com/98259577/166725783-b781e520-fedc-456d-ae08-348bedf1a7f4.png)

Tras pulsar en "Apply", obtenemos el siguiente resultado:

![image](https://user-images.githubusercontent.com/98259577/166731863-9374c8de-5885-4376-bf72-4a1335f4a2bd.png)

A continuación, lo que realmente nos interesa es conocer la composición taxonómica de nuestras lecturas (ya pre-procesadas para quedarnos solamente con las de buena calidad) a nivel de familia, algo que podemos conseguir por medio de los siguientes pasos:
1. En "Tree"->"Rank", indicamos como nivel taxonómico "Family".
2. En "Tree", hacemos click en "Show number of summarized".
3. En "Options", hacemos click en "Change LCA Parameters", tras lo cual indicamos como minimal score=60, max e-value=10e-10 y min complexity=0.5, con el fin de reducir falsos positivos.
Tras esto, obtenemos la siguiente representación gráfica en MEGAN6:

**Resultados a nivel de familia**

![image](https://user-images.githubusercontent.com/98259577/166743649-7eb4944e-0f35-480d-bc13-37153666648f.png)

---
### 4. Conclusión final


