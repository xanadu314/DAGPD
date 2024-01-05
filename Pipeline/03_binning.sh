# for bacterias
## binning

while read id
do 
    bowtie2-build ${id}.fa ${id}.bowtie2.index
done<id.txt

while read id
do 
    bowtie2 -x ${id}.bowtie2.index ${id}_1.clean2.fastq.gz ${id}2_2.clean2.fastq.gz | samtools sort -m 8000000000 -O bam -@ 2 -o - > ${id}.bam
done<id.txt

while read id
do 
    jgi_summarize_bam_contig_depths --outputDepth ./metabat2_result/${id}_depth.txt ${id}.bam
done<id.txt

while read id
do 
    metabat2 -t 8 -i ${id}.fa -a ./metabat2_result/${id}_depth.txt -o ./bins/${id}
done<id.txt


### copy all the bins to drep/ with public MAGs

## dereplicate
dereplicate drep_out  -p 32 --ignoreGenomeQuality --S_algorithm fastANI -sa 0.99 -nc 0.30 -g drep/* 
