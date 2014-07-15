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

# Reference that everything was aligned to
# We should ALWAYS use GRCh37!
#REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/ucsc.hg19.fasta"
REFERENCE="/gpfs/group1/m/mdr23/datasets/GATK/2.5/human_g1k_v37_decoy.fasta"

# number of threads
N_THREAD=1

# The base directory for the analysis.
# These can be modified with the OUTPUT_DIR and PREFIX variables
#BASE_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/input/$SITE"
BASE_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/96_control/$SITE"

# This is a file of BAM files, one per line.
BAM_LIST="$BASE_DIR/bamlist3"

#OUTPUT_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/output/VCF"
OUTPUT_DIR="$BASE_DIR/VCF"

GVCFLIST="$BASE_DIR/gvcflist"

if [ -n "$RENAME" ] && [ "$RENAME" -ne 0 ] ; then
	GVCFLIST="$GVCFLIST.renamed"
	OUTPUT_DIR="${OUTPUT_DIR}/merged"
else
	OUTPUT_DIR="${OUTPUT_DIR}/raw"
fi

PBS_DIR="/gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/pbs_output/genotypeVars/$SITE"

#PREFIX="$(echo $BASE_DIR | sed -e 's|/*$||' -e 's|.*/||')"
PREFIX="$SITE.v3"

# Use this to set specific PBS ARRAYID limits
#ARRAY_ARGS="291"

###==========================================================================================================
# DO NOT EDIT ANYTHING BELOW THIS LINE (unless you're brave)!!

if test ! -d "$BASE_DIR"; then
	echo "ERROR: Could not find base directory.  Do you have the correct site?"
	exit 1
fi

if test ! -d "$PBS_DIR"; then
	mkdir -p $PBS_DIR
fi

N_SAMPLES=$(wc -l $GVCFLIST| cut -d' ' -f1)

# Assume 1 min per sample to genotype, plus a 1/2 hour to annotate/filter
# (NOTE: pretty darn conservative estimates here; I don't want this failing!)
N_min=$(( 1 * N_SAMPLES + 30 ))

# get a list of unique IDs for the bamlist

if test ! -d "$OUTPUT_DIR"; then
	mkdir -p "$OUTPUT_DIR"
fi
	
TIME_STR=$(printf "%02d:%02d:00" $((N_min / 60)) $((N_min % 60)))
	
qsub -N geno_variants -l walltime=${TIME_STR} -l nodes=1:ppn=$N_THREAD -v PREFIX="${OUTPUT_DIR}/${PREFIX}.ALL",GVCF_LIST="$GVCFLIST",BAM_LIST="$BAM_LIST",RENAME="$RENAME",REFERENCE="$REFERENCE",N_THREAD=$N_THREAD, -w $PBS_DIR /gpfs/group1/m/mdr23/projects/eMERGE-PGX/scripts/genotypeVars.pbs

