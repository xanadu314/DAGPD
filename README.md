# DAGPD
Domestic Animals Gut Phage Database

## Overview of pipline
The scripts of metagenome and virome analysis are placed in [Pipeline](https://github.com/xanadu314/DAGPD/tree/main/Pipeline) directory.
There are three main modules in the pipeline:

- metagenomes assembly 
- generation the MAGs from metagenomes
- phages idenfication from metagenomes

### metagenomes assembly 

#### 01 data_preprocessing
Metagenomes assembly: reads trimming and host (animals and feed) reads removal.

#### 02 Assembly
Assemble short reads into contigs

### generation the MAGs from metagenomes

#### 03 binning
Process of grouping reads or contigs and assigning them to individual genome (bin)

#### 04 Taxonomy of bacteria
Bacteria taxonomic classification were determined based on [GTDBtk](http://gtdb.ecogenomic.org/)

### phages idenfication from metagenomes

#### 05 virus_annoatation
Identify virus (phages) sequences from contigs by a series tools.

#### 06 Taxonomy and Hosts of phage
Taxonomy and host predicition of phages

#### 07 Abundance
Calculate the relative abundance of phages in metagenomes.

#### 08 Gene annotation
The KEGG Orthology and pathway, CAZymes family, EggNOG Orthology,antibiotic Resistance genes annotation by [KEGG](https://www.genome.jp/tools/kofamkoala/), [dbCAN](https://bcb.unl.edu/dbCAN/download/), [EGGNOG](http://eggnog-mapper.embl.de), [CARD](https://card.mcmaster.ca/analyze/rgi) database.

#### 09 Huge phage annotation
Huge phage annotation and tRNA annotation

#### 10 crAss phage annotation
CrAss-like phages annotation by blastp and hmmer and Evolution analysis.
