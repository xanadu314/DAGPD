# crAss-like phage annotation
## filter contigs with at least 3 conserved proteins
blastp -query proteins.faa -db  ~/metavirus/10crasslike/02genefamliy/portal_1556.fa  -outfmt 6 -evalue 1e-10 -num_threads 32  -max_target_seqs 1 -out proteins.crAss_portal.blastout.txt 
blastp -query proteins.faa -db  ~/metavirus/10crasslike/02genefamliy/terL_1556.fa  -outfmt 6 -evalue 1e-10 -num_threads 32  -max_target_seqs 1 -out proteins.crAss_terL.blastout.txt 
blastp -query proteins.faa -db  ~/metavirus/10crasslike/02genefamliy/MCP_1556.fa  -outfmt 6 -evalue 1e-10 -num_threads 32  -max_target_seqs 1 -out proteins.crAss_MCP.blastout.txt 

hmmsearch --cpu 30 -E 1e-10 --tblout portal_1556.hmm proteins.faa > portal_1556.hmm.txt
hmmsearch --cpu 30 -E 1e-10 --tblout terl_1556.hmm proteins.faa > terl_1556.hmm.txt
hmmsearch --cpu 30 -E 1e-10 --tblout MCP_1556.hmm proteins.faa > MCP_1556.hmm.txt

## specific annotation
while read id
do
    prokka --outdir ${id} --prefix ${id}  --kingdom Viruses --cpus 24 --gcode 11 ${id}.fasta --locustag ${id}@ --force --hmms  ~/miniconda3/envs/prokka/db/hmm/PHROGs_VOG.hmm 
done<id.txt

## phylogenetic analysis

### TerL protein
mafft --maxiterate 1000 --thread  12 --localpair  DAGPD_crAss_Terl.faa > DAGPD_crAss_Terl.aln
seqkit seq -m 100 DAGPD_crAss_Terl.aln > DAGPD_crAss_Terl.100aa.aln
iqtree2 -s DAGPD_crAss_Terl.100aa.aln -B 1000 --bnni -T 48

## genome
roary -i 60 -e --mafft -p 24 *.gff -f roary_crAss_group_a --group_limit 1000000
mv core_gene_alignment.aln roary_crAss_group_a.aln
iqtree2 -s roary_crAss_group_a.aln -B 1000 --bnni -T 48
