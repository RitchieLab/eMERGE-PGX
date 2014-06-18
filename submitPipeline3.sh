#! /bin/bash

# EDIT THE FOLLOWING VARIABLES:

# set to 1 if you have a mapping of BAM -> ID and you want to use it
RENAME=1

#===========================
# Only edit the following variables if you know what you're doing!

if test -z "$1"; then
	echo "Error: you must provide a site as the first argument to this script"
	exit 1
fi

SITE="$1"

# Detect "combination" job based on the first character of the site.
# individual sites are ALWAYS <date(YYYYMMDD)>_<site>
# combined releases are ALWAYS <combined|merged>_<date(YYYYMMDD)>
COMBINED=$(echo $SITE | grep '^[^0-9]' | wc -l)

# Reference that everything was aligned to
# We should ALWAYS use GRCh37!
#REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/ucsc.hg19.fasta"
REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/human_g1k_v37_decoy.fasta"

# number of threads
#N_THREAD=1

# The base directory for the analysis.
# These can be modified with the OUTPUT_DIR and PREFIX variables
BASE_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/input/$SITE"
#BASE_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/96_control/$SITE"

# This is a file of BAM files, one per line.
BAM_LIST="$BASE_DIR/bamlist3"

OUTPUT_DIR="$BASE_DIR/GVCF"

#PBS_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/pbs_output/call_variants/$(echo $BASE_DIR | sed -e 's|/*$||' -e 's|.*/||')"
PBS_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/pbs_output/call_variants/$SITE"

###==========================================================================================================
# DO NOT EDIT ANYTHING BELOW THIS LINE (unless you're brave)!!

if test ! -d "$PBS_DIR"; then
	mkdir -p $PBS_DIR
fi

N_BAM=$(wc -l $BAM_LIST| cut -d' ' -f1)

GVCFLIST="$BASE_DIR/gvcflist"

	
# Assume 90 mins per sample to generate a GVCF
N_min=$(( 90 ))

# get a list of unique IDs for the bamlist
IDLIST=($(sed -e 's|^.*/||' -e 's|\.bam$||' "$BAM_LIST"))
GREP_SUFF=".bam"
if [ -n "$RENAME" ] && [ "$RENAME" -ne 0 ] ; then
	IDLIST=($(sort -k2,2 -u "$BAM_LIST" | cut -f2))
	GREP_SUFF="\$"
	OUTPUT_DIR="${OUTPUT_DIR}/renamed"
	GVCFLIST="$GVCFLIST.renamed"
else
	OUTPUT_DIR="${OUTPUT_DIR}/raw"
fi

if test -f "$GVCFLIST"; then
	echo "Warning: gvcflist exists!  Moving to gvcflist.old"
	mv "$GVCFLIST" "$GVCFLIST.old"
fi


if test ! -d "$OUTPUT_DIR"; then
	mkdir -p "$OUTPUT_DIR"
fi

for id in "${IDLIST[@]}" ; do
	
	BAM_FNS=$(grep "${id}${GREP_SUFF}" "$BAM_LIST")
	N_SAMPLES=1
	N_SAMPLES=$(printf '%s\n' "$BAM_FNS" | wc -l)
	N_min_id=$(( N_min * N_SAMPLES ))
	
	TIME_STR=$(printf "%02d:%02d:00" $((N_min_id / 60)) $((N_min_id % 60)))
	
	echo -e "${OUTPUT_DIR}/${id}.vcf.gz" >> $GVCFLIST

	CALLED=0

	# Add a check for a combination release - we don't need to re-call single files, just symlink them!!
	if [ -n "$COMBINED" ] && [ "$COMBINED" -ne 0 ] && [ "$N_SAMPLES" -eq 1 ]; then
		# Find the gvzf of interest
		GVCF_FNS=$(find /gpfs/group1/m/mdr23/projects/eMERGE-PGX/input/ -type f -name "${id}.vcf.gz")
		N_GVCF=$(printf '%s\n' "$GVCF_FNS" | wc -l)
		
		if [ "$N_GVCF" -eq 1 ]; then
			CALLED=1;
			ln -sf "${GVCF_FNS}"* "${OUTPUT_DIR}"
		fi
	fi
	
	if [ "$CALLED" -eq 0 ] ; then
		#echo "pass"
		qsub -N call_variants3 -l walltime=${TIME_STR} -v PREFIX="${OUTPUT_DIR}/${id}",BAM_LIST="$BAM_FNS",REFERENCE=$REFERENCE,RENAME=$RENAME -w $PBS_DIR /gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/callVars.pbs
	fi
	
done

