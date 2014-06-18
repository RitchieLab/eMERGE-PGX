#!/bin/bash

if [ -z "$1" ]; then
	SITE="$(date)"
else
	SITE="$1"
fi

# This script takes a VCF file on standard input and outputs an annotated VCF on stdout
PROJECT_DIR=/gpfs/group/mdr23/projects/eMERGE-PGX

snpEff -onlyTr $PROJECT_DIR/files/transcripts.txt -s $PROJECT_DIR/output/annotations/${SITE}_summary.html -canon -c $PROJECT_DIR/files/snpEff.config GRCh37.71 <&0 | cut -f 1-9
