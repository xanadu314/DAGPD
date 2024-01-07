# taxonomy prediction
python run_Speed_up.py --contigs cdhit2.ppr.vs2.morethan2.checkv.9k.fasta #phaGCN2

# host prediction

## iphop
iphop predict --fa_file cdhit2.ppr.vs2.morethan2.checkv.9k.fasta --db_dir iphop_db/Test_db_rw/  --out_dir iphop_res

## minimap2
cat ./MAG/*.fa > bac.ref.fa
minimap2 -c bac.ref.fa cdhit2.ppr.vs2.morethan2.checkv.9k.fasta > phage2bac.minimap2.paf 

## vcontact2
### chick.pig.cow.protein.fa were rename form checkv results.

vcontact2 --raw-proteins chick.pig.cow.protein.fa --proteins-fp chick.pig.cow_proteome_contig.csv --db ProkaryoticViralRefSeq85-ICTV -t 2 --output-dir chick_pig_cow #DAGPD

vcontact2 --raw-proteins chick.pig.cow.human.protein.fa --proteins-fp chick.pig.cow.human_proteome_contig.csv --db ProkaryoticViralRefSeq85-ICTV -t 2 --output-dir chick_pig_cow_human #GPD&DAGPD
