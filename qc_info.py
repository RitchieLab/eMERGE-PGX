#!/usr/bin/env python

import os
import sys

if __name__ == "__main__":
	
	# ALL numbers are a pair, one for raw, and one for filtered
	n_pos = [0, 0]
	n_var = [0, 0]
	n_snp = [0, 0]
	n_new = [0, 0]
	n_inc = [0, 0]
	n_single = [0,0]
	n_double = [0,0]
	n_obs = [0,0]
	n_new_obs = [0,0]
	
	f = file(sys.argv[1],'r')
	
	vcf_ids = {}
	
	for l in f:
		if l.startswith("#CHROM"):
			for (i, v) in enumerate(l.strip().split('\t')):
				if i > 8:
					vcf_ids.setdefault(v,[]).append(i)	
			break
	
	# calculate depth of coverage
	# we have a list of:
	# avg depth, # depth < 30x, # depth < 10x, # of missing
	id_cov = {k: [0, 0, 0, 0] for k in vcf_ids.iterkeys()}
	
	vcf_dup = {k: v for k, v in vcf_ids.iteritems() if len(v) > 1}
	
	bad_pos = []

	geno_list= {}
	
	n_lines = 0;
	for l in f:
	
		curr_line = l.strip().split('\t')
	
		curr_pos = (curr_line[0], curr_line[1])
		fmt = {v : i for (i, v) in enumerate(curr_line[8].split(':'))}
	
		gt_pos = fmt['GT']
		dp_pos = fmt['DP']
	
		bad_ids = []
		
		
		passed = (curr_line[6] == "PASS")
		
		n_pos[0] += 1
		n_pos[1] += passed

		geno_list = [max(-1,sum((int(k) for k in curr_line[g[0]].split(':')[gt_pos].replace('.','-1').split('/')))) for g in vcf_ids.itervalues()]
	#	print geno_list
	#	print len(geno_list)
		
		if curr_line[4] != ".":
		
			print curr_line[0:8]
			
			n_var[0] += 1
			n_var[1] += passed
			
			if len(curr_line[4].split(",")) == 1:
				n_snp[0] += 1
				n_snp[1] += passed

		

			count = sum((max(0,g) for g in geno_list))
			
			n_obs[0] += count
			n_obs[1] += count*passed
			
			if curr_line[2] == ".":
				n_new[0] += 1
				n_new[1] += passed
				
				n_new_obs[0] += count
				n_new_obs[1] += count*passed
				

			
			
			if sum((g>0 for g in geno_list)) == 1:
				n_single[0] += 1
				n_single[1] += passed
			if sum((g>0 for g in geno_list)) == 2:
				n_double[0] += 1
				n_double[1] += passed
		
		

		i = 0
		for (k, v) in vcf_ids.iteritems():
			dp = int(curr_line[v[0]].split(':')[dp_pos].replace('.','0'))
			id_cov[k][0] += dp
			id_cov[k][1] += (dp < 30)
			id_cov[k][2] += (dp < 10)
			id_cov[k][3] += (geno_list[i] == -1)
			i += 1
		
					
		for k, v in vcf_dup.iteritems():
			geno_list[k] = [max(-1,sum((int(k) for k in curr_line[g].split(':')[gt_pos].replace('.','-1').split('/')))) for g in v]
			if (any((g==0 for g in geno_list[k])) + any((g==1 for g in geno_list[k])) + any((g==2 for g in geno_list[k])) > 1):
				bad_ids.append(k)

		if len(bad_ids):
			bad_pos.append((curr_pos, tuple(bad_ids)))
			n_inc[0] += 1
			n_inc[1] += passed
			
		n_lines +=1


	print "Base pairs: ", n_pos
	print "Variants: ", n_var
	print "SNPs: ", n_snp
	print "Novel variants: ", n_new
	print "Singletons: ", n_single
	print "Doubletons: ", n_double
	print "Inconsistencies: ", n_inc	
	print "Number of Novel Observations", n_new_obs
	print "Number of Observations", n_obs
	print "Percentage of Novel Observations", [n_new_obs[0] / float(n_obs[0]), n_new_obs[1] / float(n_obs[0])]
	
	
	for (k,v) in id_cov.iteritems():
		print k, ":", v[0] / float(n_lines), v[1], v[2], v[3] / float(n_lines)
	#for p in bad_pos:
	#	print "%s\t%s" % ('\t'.join(p[0]), ','.join(p[1]))
