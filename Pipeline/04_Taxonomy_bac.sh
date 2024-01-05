# quality control by checkM

## checkm
checkm lineage_wf --threads 32 --tmpdir ./ --extension fa dereplicated_genomes/ checkm/ > checkM.sh.log 2>&1

### remain high-quality bins

## gtdbtk for high-quality MAGs
gtdbtk classify_wf --cpus 36 --pplacer_cpus 2 --genome_dir ./dereplicated_genomes_newid/ -x fa --out_dir gtdb_result/
