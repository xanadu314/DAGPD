# for virus
## filter contigs more than 9kb
while read id
do
    seqkit seq -m 9000 {id}.fa -o ${id}.9k.fasta
done<id.txt

## virsorter2
while read id
do virsorter run --prep-for-dramv -w ${id}.out -i ${id}.9k.fasta -j 4 all \
--include-groups "dsDNAphage,ssDNA" \
--provirus-off --max-orf-per-seq 20 all
done

## PPR_Meta
while read id
do ./PPR_Meta ${id}.9k.fa ${id}.ppr.fasta
done<id.txt

### merge PPR and vs2 results and rename as all.ppr.vs2.fasta

## cd-hit
cd-hit-est  -i all.ppr.vs2.9k.fasta -o cdhit2.ppr.vs2.9k.fasta -c 0.95 -M 640000 -T 48 -n 11 -d 0 -aS 0.95 -g 1 -sc 1 -sf 1

## what the phage
nextflow run What_the_Phage-master/phage.nf -profile local,singularity --cores 16 --max_cores 48 --sk --vf --vs --dv --mp --pp --vn --sm --identify --fasta ./cdhit2.ppr.vs2.fasta 

### filter the what the phage results and only remain the contigs annotated by at least 2 tools.

## checkv
checkv end_to_end ./cdhit2.ppr.vs2.morethan2.9k.fasta  cdhit2.ppr.vs2.morethan2.9k.fasta.out -t 32 
