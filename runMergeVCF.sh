#!/bin/bash


IN_DIR="/gpfs/group1/m/mdr23/datasets/eMERGE-PGX/20140114_nw/e_roden_pgxnorthwestern_seqcustom_102212_833/rawdataset_to_PI_CC/VCF"

SITE="20140114_nw"

REFERENCE="/gpfs/group/mdr23/datasets/GATK/2.5/human_g1k_v37_decoy.fasta"

#============= Variables based on inputs, be careful!! ===================

OUT_DIR="/gpfs/group/mdr23/projects/eMERGE-PGX/input/$SITE/VCF/"

PREFIX="$OUT_DIR/$SITE"

#============= Actual script.  Abandon all hope ye who enter here ========

PBS_OUT="/gpfs/group/mdr23/projects/eMERGE-PGX/scripts/pbs_output/merge_vcf/$SITE"

if test ! -d "$PBS_OUT"; then
	mkdir -p $PBS_OUT
fi

if test ! -d "$OUT_DIR"; then
	mkdir -p $OUT_DIR
fi

qsub -w "$PBS_OUT" -v IN_DIR=$IN_DIR,PREFIX=$PREFIX,REFERENCE=$REFERENCE /gpfs/group/mdr23/projects/eMERGE-PGX/scripts/mergeVCF.pbs
