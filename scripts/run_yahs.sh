#!/bin/bash

out=yahs_test
outdir="./outs/"
contigs="./data/test.fa.gz" # need to be indexed, i.e., test.fa.gz.fai is available
hicaln="./data/test.bed" # could be .bed, .bam or .bin file

## need juicer_tools and samtools if want to do hic plot
## juicer_tools: https://github.com/aidenlab/juicer/wiki/Download
## samtools: https://github.com/samtools/samtools
## please adjust the path to juicer_tools and samtools
## here we use 12 CPUs and 32Gb memory for juicer_tools pre - adjust it according to your device
## see more information for juicer tools https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start
## output file will be ${outdir}/${out}.hic
## the output hic file could be viewed with JuiceBox https://github.com/aidenlab/Juicebox
noplot=0
#juicer_tools="java -Xmx32G -jar ~/bin/juicer_tools_1.22.01.jar pre --threads 12"
## v1.9.9 seems much faster than v1.22.01
juicer_tools="java -Xmx32G -jar ~/bin/juicer_tools.1.9.9_jcuda.0.8.jar pre"
samtools="~/bin/samtools"

#### run yahs scaffolding
../yahs -o ${outdir}/${out} ${contigs} ${hicaln} >${outdir}/${out}.log 2>&1

if [ ${noplot} -ne 0 ]; then exit 0; fi

#### this is to generate input file for juicer_tools - non-assembly mode
## here we use 8 CPUs and 32Gb memory for sorting - adjust it according to your device
(../juicer_pre ${outdir}/${out}.bin ${outdir}/${out}_scaffolds_final.agp ${contigs}.fai 2>${outdir}/tmp_juicer_pre.log | LC_ALL=C sort -k2,2d -k6,6d -T ${outdir} --parallel=8 -S32G | awk 'NF' > ${outdir}/alignments_sorted.txt.part) && (mv ${outdir}/alignments_sorted.txt.part ${outdir}/alignments_sorted.txt)
## do juicer hic map
## prepare chromosome size file from samtools index file
# ${samtools} faidx ${outdir}/${out}_scaffolds_final.fa
# cut -f1-2 ${outdir}/${out}_scaffolds_final.fa.fai >${outdir}/${out}_scaffolds_final.chrom.sizes
## another way to prepare chromosome size file
## this is an easier way especially when we have >2G scaffolds which need scaling 
cat ${outdir}/tmp_juicer_pre.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${outdir}/${out}_scaffolds_final.chrom.sizes
(${juicer_tools} ${outdir}/alignments_sorted.txt ${outdir}/${out}.hic.part ${outdir}/${out}_scaffolds_final.chrom.sizes) && (mv ${outdir}/${out}.hic.part ${outdir}/${out}.hic)

#### this is to generate input file for juicer_tools - assembly (JBAT) mode (-a)
../juicer_pre -a -o ${outdir}/${out}_JBAT ${outdir}/${out}.bin ${outdir}/${out}_scaffolds_final.agp ${contigs}.fai 2>${outdir}/tmp_juicer_pre_JBAT.log
cat ${outdir}/tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${outdir}/${out}_JBAT.chrom.sizes
(${juicer_tools} ${outdir}/${out}_JBAT.txt ${outdir}/${out}_JBAT.hic.part ${outdir}/${out}_JBAT.chrom.sizes) && (mv ${outdir}/${out}_JBAT.hic.part ${outdir}/${out}_JBAT.hic)

