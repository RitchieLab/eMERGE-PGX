#!/usr/bin/env python

import os
import sys

if __name__ == "__main__":
	
	f = file(sys.argv[1],'r')
	
	vcf_ids = {}
	
	for l in f:
		if l.startswith("#CHROM"):
			for (i, v) in enumerate(l.strip().split('\t')):
				if i > 8:
					vcf_ids.setdefault(v,[]).append(i)	
			break
	
	vcf_dup = {k: v for k, v in vcf_ids.iteritems() if len(v) > 1}
	
	bad_pos = []

	geno_list= {}
	
	for l in f:
		bad_ids = []

		curr_line = l.strip().split('\t')
		
		curr_pos = (curr_line[0], curr_line[1])
		fmt = {v : i for (i, v) in enumerate(curr_line[8].split(':'))}
		
		gt_pos = fmt['GT']
				
		for k, v in vcf_dup.iteritems():
			geno_list[k] = [max(-1,sum((int(k) for k in curr_line[g].split(':')[gt_pos].replace('.','-1').split('/')))) for g in v]
			if (any((g==0 for g in geno_list[k])) + any((g==1 for g in geno_list[k])) + any((g==2 for g in geno_list[k])) > 1):
				bad_ids.append(k)

		if len(bad_ids):
			bad_pos.append((curr_pos, tuple(bad_ids)))
	
	for p in bad_pos:
		print "%s\t%s" % ('\t'.join(p[0]), ','.join(p[1]))
