# Asignación taxonómica con DIAMOND+MEGAN
## Víctor Manuel López Molina
## Fecha: 20/05/2022
## Ejercicio

Follow the workflow of this tutorial for the taxonomic binning of the virome reads and contigs used as homework in the unit_3.

>Reads from unit_3 homework are paired_end reads. To perform this task you can do it in several ways; joining all decontaminated and QF reads and then run DIAMOND or run DIAMOND using both files and then join DIAMOND outputs, etc..

Write a brief summary describing the bioinformatic pipeline you have followed (trimming, decontamination, improve in quality, number of reads remove in each step, etc.) and the most relevant results (two figures maximun) with the taxonomy assessment of the virome at family level.

**Optional**: run the same analysis but using the assembled contigs and try to compare the results using MEGAN (File/Compare and load both rma6 files created while making the individual analysis).

### 1. Pre-procesado
En primer lugar, partimos de los ficheros virome_R1.fastq y virome_R2.fastq que habíamos obtenido en la práctica anterior a partir del fichero virome.zip. Al tratarse de lecturas paired-end, y al deber trabajar con un único archivo, vamos a fusionar el contenido del ambos archivos en uno solo, por medio del siguiente comando:
```
cat virome_R1.fastq virome_R2.fastq > virome.fq
```

PREPROCESSING
cat virome_R1.fastq virome_R2.fastq > virome.fq
trimmomatic SE -phred33 virome.fq virome_qf.fq SLIDINGWINDOW:4:20 MINLEN:149
for f in ./Index/*bt2; do ln -s $f .; done
bowtie2 -x human-phix174 -q virome_qf.fq --un virome_qf_clean.fq -S tmp.sam

PREPARATION OF THE DATABASE
gunzip viral.*.protein.faa.gz
cat viral.*.protein.faa > viral.protein.faa
grep -c ">" *faa

REFERENCE DATABASE
diamond makedb --in viral.protein.faa -d viralproteins

ALIGNMENT TASK
diamond blastx -d viralproteins.dmnd -q virome_qf_clean.fq -o virome_qf_clean_vs_viralprotein.m8

TAKE A LOOK AT THE FILE
head virome_qf_clean_vs_viralprotein.m8

---

HASTA AQUÍ SON LOS COMANDOS, AHORA ABRES MEGAN
