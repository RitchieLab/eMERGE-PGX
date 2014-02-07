#!/usr/bin/env python

import os
import sys

if __name__ == "__main__":
	
	# ALL numbers are a pair, one for raw, and one for filtered
	n_pos = [0, 0]
	n_var = [0, 0]
	n_snp = [0, 0]
	n_new = [0, 0]
	n_samp_inc = [0, 0]
	n_ctrl_inc = [0,0]
	n_single = [0,0]
	n_double = [0,0]
	
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

	geno_dict= {}
	
	for l in f:
	
		curr_line = l.strip().split('\t')
	
		curr_pos = (curr_line[0], curr_line[1])
		fmt = {v : i for (i, v) in enumerate(curr_line[8].split(':'))}
	
		gt_pos = fmt['GT']
	
		bad_ids = []
		
		
		passed = (curr_line[6] == "PASS")
		
		n_pos[0] += 1
		n_pos[1] += passed

		if curr_line[4] != ".":
		
			geno_list = [max(-1,sum((int(k) for k in curr_line[g[0]].split(':')[gt_pos].replace('.','-1').split('/')))) for g in vcf_ids.itervalues()]
			n_var[0] += 1
			n_var[1] += passed
			
			if len(curr_line[4].split(",")) == 1:
				n_snp[0] += 1
				n_snp[1] += passed

		

			count = sum((max(0,g) for g in geno_list))
			
			if curr_line[2] == ".":
				n_new[0] += 1
				n_new[1] += passed

			
			
			if sum((g>0 for g in geno_list)) == 1:
				n_single[0] += 1
				n_single[1] += passed
			if sum((g>0 for g in geno_list)) == 2:
				n_double[0] += 1
				n_double[1] += passed
		
		

			for k, v in vcf_dup.iteritems():
				geno_dict[k] = [max(-1,sum((int(k) for k in curr_line[g].split(':')[gt_pos].replace('.','-1').split('/')))) for g in v]
				if (any((g==0 for g in geno_dict[k])) + any((g==1 for g in geno_dict[k])) + any((g==2 for g in geno_dict[k])) > 1):
					bad_ids.append(k)

			if len(bad_ids):
				bad_pos.append((curr_pos, tuple(bad_ids)))
			
				# eMERGE IDs are always all numbers
				if any((x.isdigit() for x in bad_ids)):
					n_samp_inc[0] += 1
					n_samp_inc[1] += passed
			
				if any((not x.isdigit() for x in bad_ids)):
					n_ctrl_inc[0] += 1
					n_ctrl_inc[1] += passed
			
	print "Base pairs: ", n_pos
	print "Variants: ", n_var
	print "SNPs: ", n_snp
	print "Novel variants: ", n_new
	print "Singletons: ", n_single
	print "Doubletons: ", n_double
	print "Sample Inconsistencies: ", n_samp_inc	
	print "Control Inconsistencies: ", n_ctrl_inc
	
