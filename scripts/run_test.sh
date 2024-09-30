#!/bin/bash


ROOT_DIR=$(dirname "$0")
TEST_OUTDIR=$(mktemp -d ./TEST.XXXXXX)

################################ START OF YAHS SCAFFOLDING TEST ################################
#### download the test data
test="LYZE01"
(echo "Downloading test data..."
wget -P ${TEST_OUTDIR} -q https://zenodo.org/record/7079219/files/${test}.contigs.fasta.gz && \
wget -P ${TEST_OUTDIR} -q https://zenodo.org/record/7079219/files/${test}.contigs.fasta.gz.fai && \
wget -P ${TEST_OUTDIR} -q https://zenodo.org/record/7079219/files/${test}.hic.bed.gz && \
echo "Downloading test data DONE.") || { echo "Downloading test data FAILED."; rm -rf ${TEST_OUTDIR}; exit 1; }

out="test_out"
contigs="${test}.contigs.fasta.gz" # need to be indexed, i.e., ${test}.contigs.fasta.gz.fai is presented
hicaln="${test}.hic.bed.gz" # could be .bed, .bam or .bin file

#### run yahs scaffolding
(${ROOT_DIR}/../yahs -o ${TEST_OUTDIR}/${out} ${TEST_OUTDIR}/${contigs} ${TEST_OUTDIR}/${hicaln} >/dev/null 2>&1 && echo "YaHS scaffolding DONE.") || { echo "YaHS scaffolding FAILED."; rm -rf ${TEST_OUTDIR}; exit 1; }

(${ROOT_DIR}/../agp_to_fasta ${TEST_OUTDIR}/test_out_scaffolds_final.agp ${TEST_OUTDIR}/${contigs} -o ${TEST_OUTDIR}/${out}.fa >/dev/null 2>&1 && echo "AGP to FASTA DONE.") || { echo "AGP to FASTA FAILED."; rm -rf ${TEST_OUTDIR}; exit 1; }

(${ROOT_DIR}/../juicer pre ${TEST_OUTDIR}/test_out.bin ${TEST_OUTDIR}/test_out_scaffolds_final.agp ${TEST_OUTDIR}/LYZE01.contigs.fasta.gz.fai >${TEST_OUTDIR}/test_out.aln.txt 2>/dev/null && echo "Juicer pre DONE.") || { echo "Juicer pre FAILED."; rm -rf ${TEST_OUTDIR}; exit 1; }

(${ROOT_DIR}/../juicer pre ${TEST_OUTDIR}/test_out.bin ${TEST_OUTDIR}/test_out_scaffolds_final.agp ${TEST_OUTDIR}/LYZE01.contigs.fasta.gz.fai -a -o ${TEST_OUTDIR}/test_out.JBAT >/dev/null 2>&1 && echo "Juicer pre -a DONE.") || { echo "Juicer pre -a FAILED."; rm -rf ${TEST_OUTDIR}; exit 1; }

(${ROOT_DIR}/../juicer post ${TEST_OUTDIR}/test_out.JBAT.assembly ${TEST_OUTDIR}/test_out.JBAT.liftover.agp ${TEST_OUTDIR}/LYZE01.contigs.fasta.gz -o ${TEST_OUTDIR}/test_out.JBAT >/dev/null 2>&1 && echo "Juicer post DONE.") || { echo "Juicer post FAILED."; rm -rf ${TEST_OUTDIR}; exit 1; }

rm -rf ${TEST_OUTDIR}

echo "Successful."

################################# END OF YAHS SCAFFOLDING TEST #################################

