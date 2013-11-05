#! /usr/bin/env python

import sys
import pprint
import operator


def print_dict(d):
	keys = [(-len(k), k) for k in d.keys()]
	keys.sort()
	for k in keys:
		if k[1]:
			print k[1], ":", d[k[1]]
			
def print_concord(c):
	r = [100*v / float(sum(c)) for v in c]
	print "Concordance:", "%.3f" % r[0]
	print "Heterozygous Discord:", "%.3f" % r[1]
	print "Homozygous Discord:", "%.3f" % r[2]
	print "Missing Discord:", "%.3f" % r[-1]

if __name__ == "__main__":


	vcf_list = [file(f,'r') for f in sys.argv[1:]]
	
	vcf_ids = []
	
	for f in vcf_list:
		for l in f:
			if l.startswith("#CHROM"):
				vcf_ids.append({v: i for (i, v) in enumerate(l.strip().split('\t')) if i > 8})
				break
	
	# OK, we've now got all the headers, let's double check
	
	common_ids = set.intersection(*(set(d.keys()) for d in vcf_ids))
	all_ids = set.union(*(set(d.keys()) for d in vcf_ids))
		
	id_count = {}
	for i in all_ids:
		s = sum((i in d for d in vcf_ids))
		id_count[s] = id_count.get(s,0) + 1
	
	print "Number of IDs in number of files:"
	
	keys = id_count.keys()
	keys.sort()
	for k in keys:
		print k, "files: ", id_count[k], "IDs"

	# a true/false of whether or not to advance the VCF position on the iteration of the loop
	to_advance = [1 for f in vcf_list]
	curr_pos = [('0',0) for f in vcf_list]
	curr_line = ['' for f in vcf_list]
	at_eof = [0 for f in vcf_list]
	
	prev_pos = ('0',0)
	seen_chrom = set()
	
	r_site_count = {}
	f_site_count = {}
	
	r_referent_count = {}
	f_referent_count = {}
	
	r_polyallelic_count = {}
	f_polyallelic_count = {}
	
	r_alt_discord = 0
	f_alt_discord = 0
	
	r_concord = [0, 0, 0, 0]
	f_concord = [0, 0, 0, 0]
	
	r_single_concord = [0, 0, 0, 0]
	f_single_concord = [0, 0, 0, 0]
	
	r_all_concord = [0, 0, 0, 0]
	f_all_concord = [0, 0, 0, 0]
	
	ref = ['' for f in vcf_list]
	alt = ['' for f in vcf_list]
	
	r_comp_sites = 0
	f_comp_sites = 0
	r_comp_extra = 0
	f_comp_extra = 0
	r_single = 0
	f_single = 0
	
	nonref = set((frozenset([i]) for i in ['.','']))
	
	while not all(at_eof):
		for i,f in enumerate(vcf_list):
			if to_advance[i]:
				try:
					curr_line[i] = f.next().strip().split('\t')
					curr_pos[i] = (curr_line[i][0], int(curr_line[i][1]))
					ref[i] = curr_line[i][3]
					alt[i] = curr_line[i][4].split(',')
				except StopIteration:
					at_eof[i] = 1
					to_advance[i] = 0
					curr_pos[i] = ('',-1)
					ref[i] = ''
					alt[i] = ['']
		
		cand_pos = set(curr_pos)
		
		if len(cand_pos) == 1:
			working_pos = cand_pos.pop()
		else:
			# First find all current chrom positions
			chrom_pos = [(c, p) for c,p in cand_pos if c == prev_pos[0]]
			chrom_pos.sort()
			if len(chrom_pos) >= 1:
				working_pos = chrom_pos[0]
			else:
				all_pos = [(c,p) for c,p in cand_pos]
				all_pos.sort()
				working_pos = all_pos[0]
				if working_pos[0] in seen_chrom:
					print >> sys.stderr, "WARNING: Chromosomes in VCF files may not be in identical order: please sort VCF files by chromosome and then position"
				
		seen_chrom.add(working_pos[0])
		prev_pos = working_pos
			
		if working_pos == ('', -1):
			continue
			
		#print working_pos
				
		to_advance = [c == working_pos for c in curr_pos]
		filter_status = [(v and curr_line[i][6] == "PASS") for i, v in enumerate(to_advance)]
			
		in_raw = tuple((i for i, v in enumerate(to_advance) if v))
		in_filtered = tuple((i for i, v in enumerate(filter_status) if v))
		
		r_site_count[in_raw] = r_site_count.get(in_raw,0)+1
		f_site_count[in_filtered] = f_site_count.get(in_filtered,0)+1
		
		raw_referent = tuple((j for j,x in enumerate([(v and alt[i][0] == '.') for i,v in enumerate(to_advance)]) if x))
		filt_referent = tuple((j for j,x in enumerate([(v and alt[i][0] == '.') for i,v in enumerate(filter_status)]) if x))
		
		r_referent_count[raw_referent] = r_referent_count.get(raw_referent,0)+1
		f_referent_count[filt_referent] = f_referent_count.get(filt_referent,0)+1
		
		raw_poly = tuple((j for j,x in enumerate([(v and len(alt[i]) > 1) for i,v in enumerate(to_advance)]) if x))
		filt_poly = tuple((j for j,x in enumerate([(v and len(alt[i]) > 1) for i,v in enumerate(filter_status)]) if x))
		
		r_polyallelic_count[raw_poly] = r_polyallelic_count.get(raw_poly,0)+1
		f_polyallelic_count[filt_poly] = f_polyallelic_count.get(filt_poly,0)+1
		
		r_alts = set((frozenset(alt[i]) for i,v in enumerate(to_advance) if v)) - nonref
		f_alts = set((frozenset(alt[i]) for i,v in enumerate(filter_status) if v)) - nonref
		
		r_alt_discord += (len(r_alts) > 1)
		f_alt_discord += (len(f_alts) > 1)
			
		if len(r_alts) == 1 and sum(len(a) for a in r_alts) == 1:
			r_comp_extra += (len(set((frozenset(alt[i]) for i,v in enumerate(to_advance) if v))) > 1)
			f_comp_extra += (sum(filter_status) > 1 and len(set((frozenset(alt[i]) for i,v in enumerate(filter_status) if v))) > 1) 
			r_comp_sites += 1
			f_comp_sites += (sum(filter_status) > 1)
			# list is [# homo. ref., # hetero., # homo alt, # missing]
			r_geno_count = {p : [0,0,0,0] for p in common_ids}
			f_geno_count = {p : [0,0,0,0] for p in common_ids}
			r_called = [0 for i in curr_line]
			f_called = [0 for i in curr_line]
			for i, g in enumerate(curr_line):
				if to_advance[i]:
					for p in common_ids:
						# geno will be either 0,1,2 (normal encoding), or -1 (missing)
						geno = max(-1,sum((int(k) for k in g[vcf_ids[i][p]].split(':')[0].replace('.','-1').split('/'))))
						r_called[i] += (geno > 0)
						r_geno_count[p][geno] += 1
							
						if f_alts and filter_status[i]:
							f_called[i] += (geno > 0)
							f_geno_count[p][geno] += 1
						# NOTE: do we want to treat filtered as "missing" or "not present"?
						#else:
						#	f_geno_count[p][-1] += 1
				
			r_curr_concord = [sum((sum((n*(n-1)/2 for n in v)) for v in r_geno_count.itervalues())),
							  sum((v[0]*v[1] + v[1]*v[2] for v in r_geno_count.itervalues())),
							  sum((v[0]*v[2] for v in r_geno_count.itervalues())),
							  sum((v[-1]*(v[0] + v[1] + v[2]) for v in r_geno_count.itervalues()))]
			
			prev_all = sum(f_all_concord)
			
			r_all_concord[2] += sum((v[0]>0 and v[2]>0) for v in r_geno_count.itervalues())
			r_all_concord[1] += sum(((not (v[0]>0 and v[2]>0)) and ((v[0]>0 and v[1]>0) or (v[1]>0 and v[2]>0))) for v in r_geno_count.itervalues())
			r_all_concord[-1] += sum((sum(n>0 for n in v[0:3]) == 1 and (not all(to_advance) or v[3]>0)) for v in r_geno_count.itervalues())
			r_all_concord[0] += sum(((all(to_advance) and sum(n>0 for n in v) == 1) or sum(n>0 for n in v[0:3]) == 0) for v in r_geno_count.itervalues())
			
			f_all_concord[2] += sum((v[0]>0 and v[2]>0) for v in f_geno_count.itervalues())
			f_all_concord[1] += sum(((not (v[0]>0 and v[2]>0)) and ((v[0]>0 and v[1]>0) or (v[1]>0 and v[2]>0))) for v in f_geno_count.itervalues())
			f_all_concord[-1] += sum((sum(n>0 for n in v[0:3]) == 1 and (not all(filter_status) or v[3]>0)) for v in f_geno_count.itervalues())
			f_all_concord[0] += sum(((all(filter_status) and sum(n>0 for n in v) == 1) or sum(n>0 for n in v[0:3]) == 0) for v in f_geno_count.itervalues())	
							  
			r_concord = map(operator.add, r_concord, r_curr_concord)
			
			f_curr_concord = [sum((sum((n*(n-1)/2 for n in v)) for v in f_geno_count.itervalues())),
							  sum((v[0]*v[1] + v[1]*v[2] for v in f_geno_count.itervalues())),
							  sum((v[0]*v[2] for v in f_geno_count.itervalues())),
							  sum((v[-1]*(v[0] + v[1] + v[2]) for v in f_geno_count.itervalues()))]
			
			f_concord = map(operator.add, f_concord, f_curr_concord)
			
			if any((v == 1 for v in r_called)):
				r_single += 1
				r_single_concord = map(operator.add, r_single_concord, r_curr_concord)
			
			if any((v == 1 for v in f_called)):
				f_single += 1
				f_single_concord = map(operator.add, f_single_concord, f_curr_concord)
		
	
	print "Raw Results"
	print "==========="
	print "Sites by file:"
	
	print_dict(r_site_count)
	
	print "Referent sites by file:"
	
	print_dict(r_referent_count)

	print "Polyallelic sites by file:"
	
	print_dict(r_polyallelic_count)
	
	print "Number of discordant ALT:", r_alt_discord
	
	print "\nPairwise Concordance Results"
	print "------------------"
	print "Total sites:", r_comp_sites
	print "Total matching sites", r_comp_sites - r_comp_extra
	
	print_concord(r_concord)
	
	print "\nSingleton Concordance"
	print "---------------------"
	print "Total Singletons:", r_single
		
	print_concord(r_single_concord)
	
	print "\nOverall Concordance"
	print "-------------------"
	
	print_concord(r_all_concord)
	
	print "\n\n"
	
	print "Filtered Results"
	print "================"
	print "Sites by file:"
	
	print_dict(f_site_count)
	
	print "Referent sites by file:"
	
	print_dict(f_referent_count)

	print "Polyallelic sites by file:"
	
	print_dict(f_polyallelic_count)
	
	print "Number of discordant ALT:", f_alt_discord
	
	print "\nPairwise Concordance Results"
	print "------------------"
	print "Total sites:", f_comp_sites
	print "Total matching sites", f_comp_sites - f_comp_extra
	
	print_concord(f_concord)
	
	print "\nSingleton Concordance"
	print "---------------------"
	print "Total Singletons:", f_single
	
	print_concord(f_single_concord)
	
	print "\nOverall Concordance"
	print "-------------------"
	
	print_concord(f_all_concord)
	
