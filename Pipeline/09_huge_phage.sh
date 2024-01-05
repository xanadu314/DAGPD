# huge phage annotation

## filter genome size more than 200kb
seqkit seq -m 200000 cdhit2.ppr.vs2.morethan2.checkv.9k.fasta >  cdhit2.ppr.vs2.morethan2.checkv.200k.fasta

## circularization estimate by vs2
virsorter run  -w huge_phage -i cdhit2.ppr.vs2.morethan2.checkv.200k.fasta -j 4 all


## CRISPRFinder
perl ~/download/crisprfinder/CRISPRCasFinder/CRISPRCasFinder.pl -so ~/download/crisprfinder/CRISPRCasFinder/sel392v2.so  -drpt ~/download/crisprfinder/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts ~/crisprfinder/CRISPRCasFinder/supplementary_files/Repeat_List.csv -def G -out CRISPR_out -in CRISPRCasFinder/cdhit2.ppr.vs2.morethan2.checkv.200k.fasta

## tRNA 
tRNAscan-SE -B -o tRNA.out -m tRNA.stats cdhit2.ppr.vs2.morethan2.checkv.200k.fasta
