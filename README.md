# YaHS: yet another Hi-C scaffolding tool

## Overview
YaHS is scaffolding tool using Hi-C data. It relies on a new algothrim for contig joining detection which considers the topological distribution of Hi-C signals aiming to distingush real interaction signals from mapping nosies. YaHS has been tested in a wide range of genome assemblies. Compared to other Hi-C scaffolding tools, it usually generates more contiguous scaffolds - especially with a higher N90 and L90 statistics. It is also super fast - takes less than 5 minutes to reconstruct the human genome from an assembly of 5,483 contigs with ~45X Hi-C data. See the [poster](https://drive.google.com/file/d/1jPhSi1k4ROGb1OSfKDurFIUvn-x1G6Kg/view?usp=sharing) presented in the [Bioversity Genemics 2021 conference](https://www.darwintreeoflife.org/news_item/biodiversity-genomics-2021-sequencing-genomes-across-the-planet/) for more information.

## Installation
You need to have a C compiler, GNU make and zlib development files installed. Download the source code from this repo or with `git clone https://github.com/c-zhou/yahs.git`. Then type `make` in the source code directory to compile.

## Run YaHS
YaHS has two required inputs: a FASTA format file with contig sequences which need to be indexed (with [samtools faidx](http://www.htslib.org/doc/samtools-faidx.html) for example) and a BAM/BED/BIN file with the alignment results of Hi-C reads to the contigs. A recommended way to generate the alignment file is to use the [Arima Genomics' mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline). It is also recommened to mark PCR/optical duplicates. Several tools are available out there for marking duplicates such as `bammarkduplicates2` from [biobambam2](https://bio.tools/biobambam) and `MarkDuplicates` from [Picard](https://broadinstitute.github.io/picard/). The resulted BAM file need to be sorted by read name before feeding to YaHS. This could be done with [samtools sort](http://www.htslib.org/doc/samtools-sort.html) with `-n` option. The BED format file can be generated from the BAM file (with [bedtools bamtobed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) for example), but do NOT forget to filter out the PCR/optical duplicates. The BED format is accepted mainly to keep consistent with other Hi-C scaffolding tools. There is no need to convert the BAM format to BED format unless you want to compare YaHS to other tools. The BIN format is a binary format specific to YaHS. If the input file is BAM (with `.bam` extension) or BED (with `.bed` extension) format, the first step of YaHS is to convert them to BIN format (with `.bin` extension). This is to save running time as multiple rounds of file IO are needed during the scaffolding process. If you have run YaHS and need to rerun it, the BIN file in the output directory could be reused to save some time - although might be just a few minutes.

Here is an example to run YaHS,

    yahs contigs.fa hic-to-contigs.bam

The outputs include several [AGP format](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) files and a FASTA format file. The `*_inital_break_[0-9]{2}.agp` AGP files are for initial assembly error corrections. The `*_r[0-9]{2}.agp` and related `*_r[0-9]{2}_break.agp` AGP files are for scaffolding results in each round. The `*_scaffolds_final.agp` and `*_scaffolds_final.fa` files are for the final scaffolding results.

There are some optional parameters.

With `-o` option, you can specify the prefix of the output files. It is `./yash.out` by default. If a directory structure is included, the directory needs to be existed.

With `-a` option, you can specify a AGP format file to ask YaHS to do scaffolding with the scaffolds in the AGP file as the start point.

With `-r` option, you can specify a range of resultions. It is `50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000` by default and the upper limit is automatically adjusted with the genome size.

With `--no-contig-ec` option, you can skip the initial assembly error correction step. With `-a` option, this will be set automatically.

With `--no-scaffold-ec` option, YaHS will skip the scaffolding error check in each round. There will no `*_r[0-9]{2}_break.agp` AGP output files.

## Limitations
YaHS is still under development and only tested with genome assemblies limited to a few species. You are welcomed to use it and report failures. Any suggestions would be appreciated.
