#! /bin/bash

# EDIT THE FOLLOWING VARIABLES:

# Reference that everything was aligned to
#REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/ucsc.hg19.fasta"
REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/human_g1k_v37_decoy.fasta"

# The base directory for the analysis.  (bgzipped) VCF files will go in "$BASE_DIR/data/PGX_batch_N.vcf.gz"
# These can be modified with the OUTPUT_DIR and PREFIX variables
#BASE_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/96_control/mtsinai/"
BASE_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/96_control/uw/"
#BASE_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/input/mayo_09262013"

# Set this environment variable is using reduced BAMs (leave unset if using raw BAMs)
#REDUCED=1

# This is a file of BAM files, one per line.
BAM_LIST="$BASE_DIR/raw_bamlist"
#BAM_LIST="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/96_control/mtsinai/bamlist"
#BAM_LIST="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/96_control/cidr/raw_bamlist"

# Number of batches to run.  If you'd prefer to specify the number of samples / batch, set to 0
N_BATCHES=1

# Number of samples / batch.  WARNING: If you use this parameter, the last batch may have significantly fewer samples
# NOTE: this is only used when N_BATCHES < 1
N_SAMPLES=500

#OUTPUT_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/output/VCF/raw"
OUTPUT_DIR="$BASE_DIR/VCF"
PBS_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/pbs_output/call_variants/$(echo $BASE_DIR | sed -e 's|/*$||' -e 's|.*/||')"

#PREFIX="$(echo $BASE_DIR | sed -e 's|/*$||' -e 's|.*/||')"
PREFIX="raw"

# Directory to put grouped BAM file lists
SHUF_DIR="$HOME/scratch/PGX/${PREFIX}/"

# Use this to set specific PBS ARRAYID limits
#ARRAY_ARGS="291"

###==========================================================================================================
# DO NOT EDIT ANYTHING BELOW THIS LINE (unless you're brave)!!

if test ! -d "$PBS_DIR"; then
	mkdir -p $PBS_DIR
fi

if test ! -d "$SHUF_DIR"; then
	mkdir -p $SHUF_DIR
fi

N_BAM=$(wc -l $BAM_LIST| cut -d' ' -f1)

BAM_FN_ARRAY=( "$BAM_LIST" ) 

if test "$N_BATCHES" -ne 1; then
	if test "$N_BATCHES" -lt 1; then
		slice=$N_SAMPLES
	else
		slice=$((N_BAM / N_BATCHES + ( N_BAM % N_BATCHES > 0 ) ))
	fi
		
	echo "Submitting in batches of $slice"
	
	if test ! -d "${SHUF_DIR}"; then 
		mkdir "${SHUF_DIR}"; 
	else
		rm -rf "${SHUF_DIR}/bamlist.*"	
	fi
	
	cat $BAM_FN | shuf | split -a $((N_BATCHES / 10 + 1)) -l $slice -d - "${SHUF_DIR}/bamlist.${PREFIX}"
	
	
	OIFS="$IFS"
	IFS=$'\n'
	BAM_FN_ARRAY=( $(ls -1 "${SHUF_DIR}/bamlist.*") )
	IFS="$OIFS"
fi

N_BATCH="${#BAM_FN_ARRAY[@]}"
if test "$N_BATCH" -ne 1; then
	PREFIX="${PREFIX}_batch"
else
	USE_PREFIX="$PREFIX"	
fi

for i in "${!BAM_FN_ARRAY[@]}"; do 
	if test "$N_BATCH" -ne 1; then
		USE_PREFIX="${PREFIX}_$i"
	fi
	
	# Calculate the time needed
	N_SAMPLES=$(wc -l ${BAM_FN_ARRAY[$i]} | cut -d' ' -f1)
	N_sec=$(awk "END {printf \"%d\", $N_SAMPLES*(585*log($N_SAMPLES)/log(10) + 2350) }" </dev/null)
	N_min=$(( N_sec / 60 + 180 ))

	# If using reduced BAMs, we can expect ~12x speedup
	if test ! -z "$REDUCED"; then
		N_min=$(( N_min / 12 ))
	fi

	TIME_STR=$(printf "%02d:%02d:00" $((N_min / 60)) $((N_min % 60)))
	
	qsub -N call_variants -l walltime=${TIME_STR} -v PREFIX=${OUTPUT_DIR}/${USE_PREFIX},BAM_LIST=${BAM_FN_ARRAY[$i]},REFERENCE=$REFERENCE -w $PBS_DIR runPipeline.pbs
done
