#!/bin/bash


IN_DIR=""

SITE=""

REFERENCE="/gpfs/group/mdr23/datasets/GATK/2.5/human_g1k_v37_decoy.fasta"

#============= Variables based on inputs, be careful!! ===================

PREFIX="/gpfs/group/mdr23/projects/eMERGE-PGX/input/$SITE/VCF/$SITE"

#============= Actual script.  Abandon all hope ye who enter here ========

PBS_OUT="/gpfs/group/mdr23/projects/eMERGE-PGX/scripts/pbs_output/merge_vcf/$SITE"

if test ! -d "$PBS_OUT"; then
	mkdir -p $PBS_OUT
fi

qsub -w "$PBS_OUT" -v IN_DIR=$IN_DIR,PREFIX=$PREFIX,REFERENCE=$REFERENCE /gpfs/group/mdr23/projects/eMERGE-PGX/scripts/mergeVCF.pbs
