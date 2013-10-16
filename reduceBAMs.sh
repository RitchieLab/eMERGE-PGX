#! /bin/bash

# EDIT THE FOLLOWING VARIABLES:

PBS_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/pbs_output/reduce_bam/96_control_mtsinai"
OUTPUT_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/96_control/mtsinai/reduced_bam"

# Reference that everything was aligned to
REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/ucsc.hg19.fasta"
#REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/human_g1k_v37_decoy.fasta"

INPUT_DIR="/gpfs/group1/m/mdr23/FTP_root/eMERGE/Uploads/MtSinai/PGRNseq2013/BAMs"

###==========================================================================================================
# DO NOT EDIT ANYTHING BELOW THIS LINE (unless you're brave)!!

if test ! -d "$OUTPUT_DIR"; then
	mkdir -p $OUTPUT_DIR
fi

if test ! -d "$PBS_DIR"; then
	mkdir -p $PBS_DIR
fi

for d in $INPUT_DIR; do
	N_BAMS=$(ls -1 $d/*.bam | wc -l)
	qsub -v BAM_DIR=$d,OUT_DIR=$OUTPUT_DIR,REFERENCE=$REFERENCE -t 1-$N_BAMS -w $PBS_DIR reduceBAMDir.pbs
done
