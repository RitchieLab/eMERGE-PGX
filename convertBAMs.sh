#! /bin/bash

# EDIT THE FOLLOWING VARIABLES:

SITE="mayo"

INPUT_DIR="/gpfs/group1/m/mdr23/datasets/eMERGE-PGX-control/mayo/last_3"

#============================================================================================================

PBS_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/pbs_output/convert_bam/$SITE"
OUTPUT_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/input/$SITE/converted_bam"

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
	qsub -v BAM_DIR=$d,OUT_DIR=$OUTPUT_DIR -t 1-$N_BAMS -w $PBS_DIR /gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/convertBAMDir.pbs
done
