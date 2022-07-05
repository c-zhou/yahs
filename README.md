# YaHS: yet another Hi-C scaffolding tool [![DOI](https://zenodo.org/badge/411044095.svg)](https://zenodo.org/badge/latestdoi/411044095)

## Overview
YaHS is a scaffolding tool using Hi-C data. It relies on a new algothrim for contig joining detection which considers the topological distribution of Hi-C signals aiming to distingush real interaction signals from mapping nosies. YaHS has been tested in a wide range of genome assemblies. Compared to other Hi-C scaffolding tools, it usually generates more contiguous scaffolds - especially with a higher N90 and L90 statistics. It is also super fast - takes less than 5 minutes to reconstruct the human genome from an assembly of 5,483 contigs with ~45X Hi-C data. See the [poster](https://drive.google.com/file/d/1jPhSi1k4ROGb1OSfKDurFIUvn-x1G6Kg/view?usp=sharing) presented in the [Bioversity Genemics 2021 conference](https://www.darwintreeoflife.org/news_item/biodiversity-genomics-2021-sequencing-genomes-across-the-planet/) for more information.

## Installation
You need to have a C compiler, GNU make and zlib development files installed. Download the source code from this repo or with `git clone https://github.com/c-zhou/yahs.git`. Then type `make` in the source code directory to compile.

## Run YaHS
YaHS has two required inputs: a FASTA format file with contig sequences which need to be indexed (with [samtools faidx](http://www.htslib.org/doc/samtools-faidx.html) for example) and a BAM/BED/BIN file with the alignment results of Hi-C reads to the contigs. A recommended way to generate the alignment file is to use the [Arima Genomics' mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline). The resulted BAM file is recommened to mark PCR/optical duplicates before feeding to YaHS. Several tools are available out there for marking duplicates such as `bammarkduplicates2` from [biobambam2](https://bio.tools/biobambam) and `MarkDuplicates` from [Picard](https://broadinstitute.github.io/picard/). The BED format is accepted mainly to keep consistent with other Hi-C scaffolding tools such as [SALSA2](https://github.com/marbl/SALSA). Each line of the BED file should contain at least four columns, i.e., contig name the read mapped to, the start position of the alignment, the end position of the alignment and the read name. The first and last read from a read pair is optionally marked by '/1' and '/2' suffix to the read name. All the information after the fourth column are ignored. Each read pair should be placed in two consecutive lines. The BED format file can be generated from the BAM file with [bedtools bamtobed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) for example. There is no need to convert the BAM format to BED format unless you want to compare YaHS to other tools. The BIN format is a binary format specific to YaHS. If the input file is BAM (with `.bam` extension) or BED (with `.bed` extension) format, the first step of YaHS is to convert them to BIN format (with `.bin` extension). This is to save running time as multiple rounds of file IO are needed during the scaffolding process. If you have run YaHS and need to rerun it, the BIN file in the output directory could be reused to save some time - although might be just a few minutes.

> **_NOTE 1:_** The input BAM could either sorted by read names ([samtools sort](http://www.htslib.org/doc/samtools-sort.html) with `-n` option) or not. The behaviours of the program are slightly different, which might lead to slightly different scaffolding results. For a BAM input sorted by read names, with each mapped read pair, a Hi-C link is counted between the middle positions of the read alignments; while for a BAM input sorted by coordinates or unsorted, Hi-C links are counted between the start positions of the read alignments. Also, for a BAM input not sorted by read names, the mapping quality filtering is suppressed (`-q` option).

> **_NOTE 2:_** The BAM file used to genereate BED file need to be filtered out unmapped reads, supplementary/secondary alignment records, and PCR/optical duplicates, and sorted by read names (otherwise the resulted BED file need to be sorted by the read name column).

Here is an example to run YaHS,

    yahs contigs.fa hic-to-contigs.bam

The outputs include several [AGP format](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) files and a FASTA format file. The `*_inital_break_[0-9]{2}.agp` AGP files are for initial assembly error corrections. The `*_r[0-9]{2}.agp` and related `*_r[0-9]{2}_break.agp` AGP files are for scaffolding results in each round. The `*_scaffolds_final.agp` and `*_scaffolds_final.fa` files are for the final scaffolding results.

There are some optional parameters.

With `-o` option, you can specify the prefix of the output files. It is `./yash.out` by default. If a directory structure is included, the directory needs to be existed.

With `-a` option, you can specify a AGP format file to ask YaHS to do scaffolding with the scaffolds in the AGP file as the start point.

With `-r` option, you can specify a range of resultions. It is `50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000` by default and the upper limit is automatically adjusted with the genome size.

With `-e` option, you can specify the restriction enzyme(s) used by the Hi-C experiment. For example, `GATC` for the DpnII restriction enzyme used by the Dovetail Hi-C Kit; `GATC,GANT` and `CGATC,GANTC,CTNAG,TTAA` for Arima genomics 2-enzyme and 4-enzyme protocol, respectively. Sometimes, the specification of enzymes may not change the scaffolding result very much if not make it worse, especically when the base quality of the assembly is not very good, e.g., assembies constructed from noisy long reads.

With `-l` option, you can specify the minimum contig length included for scaffolding.

With `-q` option, you can set the minimum read mapping quality (for BAM input only).

With `--no-contig-ec` option, you can skip the initial assembly error correction step. With `-a` option, this will be set automatically.

With `--no-scaffold-ec` option, YaHS will skip the scaffolding error check in each round. There will be no `*_r[0-9]{2}_break.agp` AGP output files.

## Generate HiC contact maps
YaHS offers some auxiliary tools to help generating HiC contact maps for visualisation. A demo is provided in the bash script `scripts/run_yahs.sh`. To generate and visualise a HiC contact map, the following tools are required.

* samtools: https://github.com/samtools/samtools
* juicer_tools: https://github.com/aidenlab/juicer/wiki/Download
* Juicebox: https://github.com/aidenlab/Juicebox

The first step is to convert the HiC alignment file (BAM/BED/BIN) to a file required by `juicer_tools` using the tool `juicer pre` provided by YaHS. To save time, BIN file is recommended which has already been generated in the scaffolding step. Here is an example bash command:

    (juicer pre hic-to-contigs.bin scaffolds_final.agp contigs.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

The tool `juicer pre` takes three positional parameters: the alignments of HiC reads to contigs, the scaffold AGP file and the contig FASTA index file. With `-o` option, it will write the results to a file. Here, the outputs are directed to `stdout` as we need a sorted (by scaffold names) file for `juicer_tools`.

For sorting, we use 8 threads, 32Gb memory and the current directory for temporaries. You might need to adjust these settings according to your device.

The next step is to generate HiC contact matrix using `juicer_tools`. Here is an example bash command:

    (java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic)

The `juicer_tools`'s `pre` command takes three positional parameters: the sorted alignment file generated in the first step, the output file name and the file for scaffold sizes. The file for scaffold sizes should contain two columns - scaffold name and scaffold size, which can be taken from the first two columns of the FASTA index file.

Finally, the output file `out.hic` could be used for visualisation with Juicebox. More information about `juicer_tools` and Juicebox can be found [here]( https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start).

## Manual curation with Juicebox (JBAT)
You can generate a HiC contact mapfile that can be loaded by Juicebox (JBAT) for manual editing with `juicer pre` by adding `-a` parameter. For example, 
    
    juicer pre -a -o out_JBAT hic-to-contigs.bin scaffolds_final.agp contigs.fa.fai >out_JBAT.log 2>&1
    
where `hic-to-contigs.bin` (can be replaced with the original BED/BAM file with some sacrifice in running time) and `scaffolds_final.agp` are the outputs of YaHS. This results in five output files,

    out_JBAT.txt             - BED format file for hic links
    out_JBAT.liftover.agp    - coordinate file between new and old contigs
    out_JBAT.assembly        - assembly annotation file for Juicebox
    out_JBAT.assembly.agp    - AGP file contains same information as the assembly annotation file. Not a real AGP file as there are no gaps.
    out_JBAT.log             - the output log file
   
You will then need to run `juicer_tools pre` with `out_JBAT.txt` file.

    (java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

There will be a line like `[I::main_pre] JUICER_PRE CMD: java -Xmx36G -jar ${juicer_tools} pre out_JBAT.txt out_JBAT.hic <(echo "assembly 183277074")` in the `out_JBAT.log` file telling you how to run juicer_tools. The `<(echo "assembly 183277074")` part can be replaced by a chromosome size file containing only one line - assembly 183277074.

The `out_JBAT.hic` can be loaded into Juicebox along with `out_JBAT.assembly` for manual editing.
> **_NOTE 3:_**  if your total assembly size is larger than 2Gb, there will be a scale factor applied to the HiC contact map file to make it loadable by Juicebox. This scale factor can be found in the log file (out_JBAT.log) and always a power of two, e. g. `[I::main_pre] scale factor: 4`. You need to set this parameter in Juicebox through `Assembly > Set Scale`. Otherwise the HiC contact map and assembly file will not match.
    
Once completed editing, there should be a file named something like `out_JBAT.review.assembly` generated by Juicebox, which can be fed into `juicer post` command to generate AGP and FASTA files for the final genome assembly. You also need the `out_JBAT.liftover.agp` coordinate file previously generated with `juicer pre` command.

    juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp contigs.fa

This will end up with two files `out_JBAT.FINAL.agp` and `out_JBAT.FINAL.fa`. Together with `hic-to-contigs.bin` or the original BED/BAM file, you can regenerate a HiC contact map for the final assembly as described in the previous section.

You can find more information about mannual editing with Juicebox here [Issue 4](https://github.com/c-zhou/yahs/issues/4).

## Other tools
* ***juicer*** is a tool used to quickly generate HiC alignment file required for HiC contact map generation with tools like [Juicebox](https://github.com/aidenlab/Juicebox), [PretextMap](https://github.com/wtsi-hpag/PretextMap) and [Higlass](https://github.com/higlass/higlass) (`juicer pre`). It can be also used to generate AGP and FASTA files after manual editing with Juicebox JBAT (`juicer post`).
* ***agp_to_fasta*** creates a FASTA file from a AGP file. It takes two positional parameters: the AGP file and the contig FASTA file. By default, the output will be directed to `stdout`. You can write to a file with `-o` option. It allows changing the FASTA line width with `-l` option, which by default is 60. If the AGP file contains sequence components of unknown orientations ('?', '0' or 'na' identifiers, see [AGP format](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/)), you will need `-u` option, with which components with unknown orientation are treated as if they had '+' orientation.

## Limitations
* In rare cases, YaHS has been seen making telomere-to-telomere false joins.
* YaHS can only handle up to 45,000 contigs. Consider excluding short contigs from scaffolding (with `-l` option) if the contig number exceeds this limit.
* The memory consumption might be very high when the genome size is large and at the same time the contig number is large. In this case, consider starting from a lower resolution (larger `-r` values).

## Citation

Chenxi Zhou, Shane A. McCarthy, Richard Durbin. YaHS: yet another Hi-C scaffolding tool. 2022. bioRxiv doi: https://doi.org/10.1101/2022.06.09.495093
