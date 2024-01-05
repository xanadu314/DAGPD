# vOTU abundance
## only calculate abundance in chicken, pig, and ruminant
salmon index -t chick.phage.fasta -i chick.phage.salmon.index
salmon index -t pig.phage.fasta -i pig.phage.salmon.index
salmon index -t cow.phage.fasta -i cow.phage.salmon.index

while read id1 id2
do salmon quant --validateMapping --minScoreFraction 0.9 --meta -p 12 -l IU -i chick.phage.index -1 ${id1}_1.fastq.clean.gz -2 ${id1}_2.fastq.clean.gz  -o ${id2}/${id1}.phage.quant
done<id12.txt

while read id1 id2
do salmon quant --validateMapping --minScoreFraction 0.9 --meta -p 12 -l IU -i pig.phage.index -1 ${id1}_1.fastq.clean.gz -2 ${id1}_2.fastq.clean.gz  -o ${id2}/${id1}.phage.quant
done<id12.txt

while read id1 id2
do salmon quant --validateMapping --minScoreFraction 0.9 --meta -p 12 -l IU -i cow.phage.index -1 ${id1}_1.fastq.clean.gz -2 ${id1}_2.fastq.clean.gz  -o ${id2}/${id1}.phage.quant
done<id12.txt
