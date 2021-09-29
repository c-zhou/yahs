#!/bin/bash

out=yash.test
outdir="./outs/"
contigs="./data/test.fa" # need to be indexed, i.e., test.fa.fai is available
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
juicer_tools="java -jar -Xmx32G ~/bin/juicer_tools_1.22.01.jar pre --threads 12"
samtools="~/bin/samtools"

## run yahs scaffolding
../yahs -o ${outdir}/${out} ${contigs} ${hicaln} >${outdir}/${out}.log 2>&1

if [ ${noplot} -ne 0 ]; then exit 0; fi

## this is to generate input file for juicer_tools
## here we use 8 CPUs and 32Gb memory for sorting - adjust it according to your device
(../juicer_pre ${outdir}/${out}.bin ${outdir}/${out}_scaffolds_final.agp ${contigs}.fai | sort -k2,2d -k6,6d -T ${outdir} --parallel=8 -S32G | awk 'NF' > ${outdir}/alignments_sorted.txt.part) && (mv ${outdir}/alignments_sorted.txt.part ${outdir}/alignments_sorted.txt)

## do juicer hic map
${samtools} faidx ${outdir}/${out}_scaffolds_final.fa
(${juicer_tools} ${outdir}/alignments_sorted.txt ${outdir}/${out}.hic.part ${outdir}/${out}_scaffolds_final.fa.fai) && (mv ${outdir}/${out}.hic.part ${outdir}/${out}.hic) && (rm ${outdir}/alignments_sorted.txt)

