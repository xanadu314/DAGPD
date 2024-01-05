# Gene Catalog

## COG annotation by EGGNOG online

## KEGG
hmmsearch --cpu 30 -E 0.001  --tblout kofam_profiles.hmm chick.pig.cow.protein.fa > kofam_out.txt 

## ARGs
rgi main -i chick.pig.cow.protein.fa -t protein -o output.card.out.txt --clean -n 32

## Cazy
run_dbcan chick.pig.cow.protein.fa protein --out_dir dbcan_out
