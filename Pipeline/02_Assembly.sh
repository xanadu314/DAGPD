## Assemble cleanData to contigs
while read id
do 
    megahit --min-contig-len 3000 --tmp-dir ./ --memory 0.3 --num-cpu-threads 40 -1 ${id}_1.clean2.fastq.gz -2 ${id}2_2.clean2.fastq.gz -o ./megahit_result/${id}
done<id.txt

### rename results name
while read id1 id2
do
	mv ${id1} ${id2}.fa
done<id12.txt

### rename contigs name
while read id
do sed -i 's/>/>${id}/g' ${id}.fa
done<id.txt
