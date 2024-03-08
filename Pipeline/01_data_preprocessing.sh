#!/bin/bash

##QC step1: filter low quality sequecne and remove adapters by fastp
while read id
do
    fastp -i ${id}_1.fastq.gz -I ${id}_2.fastq.gz -o ${id}_1.clean.fastq.gz -O ${id}_2.clean.fastq.gz -j ${id}.fastp.josn -h ${id}.fastp.html
done<id.txt

##QC step2: filter host sequence by bowtie2
###index of reference genome
bowtie2-build reference.fna reference.db
while read id
do
    bowtie2 --very-sensitive-local -x ~/genome/reference.db -p 2 -1 ${id}_1.fastq.clean.gz -2 ${id}_2.fastq.clean.gz  --un-conc-gz ${id}_clean2.gz >/dev/null 
done<id.txt

## QC step3: Generate Reports
while read id
do
    fastp -i ${id}_clean2.1.gz -I ${id}clena2.2.gz -o ${id}_1.clean2.fastq.gz -O ${id}_2.clean2.fastq.gz -j ${id}.fastp.final.json -h ${id}.fastp.final.html
done<id.txt
