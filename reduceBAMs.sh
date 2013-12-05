#! /bin/bash

# EDIT THE FOLLOWING VARIABLES:

SITE="20131105_vanderbilt"

INPUT_DIR="/gpfs/group1/m/mdr23/datasets/eMERGE-PGRN/20131105_vanderbilt/e_roden_pgxvu_seqcustom_808/rawdataset_to_PI_CC/BAM_BAI/"

# Reference that everything was aligned to
REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/human_g1k_v37_decoy.fasta"
#REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/ucsc.hg19.fasta"

#============================================================================================================

PBS_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/pbs_output/reduce_bam/$SITE"
OUTPUT_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/input/$SITE/reduced_bam"

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
	qsub -v BAM_DIR=$d,OUT_DIR=$OUTPUT_DIR,REFERENCE=$REFERENCE -t 1-$N_BAMS -w $PBS_DIR /gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/reduceBAMDir.pbs
done
