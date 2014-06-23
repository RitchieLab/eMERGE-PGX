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

	NoneSet = set([None])

	vcf_list = [file(f,'r') for f in sys.argv[1:3]]
	vcf_ids = []
	
	for f in vcf_list:
		for l in f:
			if l.startswith("#CHROM"):
				vcf_ids.append({v: i for (i, v) in enumerate(l.strip().split('\t')) if i > 8})
#				print l.strip().split('\t')
#				print len(l.strip().split('\t'))
#				print [(i, v) for (i, v) in enumerate(l.strip().split('\t')) if i > 8]
				break				

	
	# OK, we've now got all the headers, let's double check
	
	common_ids = set.intersection(*(set(d.keys()) for d in vcf_ids))
	
#	print [i for i in common_ids]
#	print [(i, [(x, vcf_ids[i][x]) for x in common_ids]) for i in xrange(len(vcf_ids))]
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
	
	ref = [None for f in vcf_list]
	alt = [None for f in vcf_list]
	allele = [None for f in vcf_list]
	
	gt_idx = [None for f in vcf_list]
	ft_idx = [None for f in vcf_list]

	r_concord = [0, 0, 0, 0]
	f_concord = [0, 0, 0, 0]
	
	r_snp_concord = [0, 0, 0, 0]
	f_snp_concord = [0, 0, 0, 0]
	

	
	while not all(at_eof):
		for i,f in enumerate(vcf_list):
			if to_advance[i]:
				try:
					curr_line[i] = f.next().strip().split('\t')
					curr_pos[i] = (curr_line[i][0], int(curr_line[i][1]))
					ref[i] = curr_line[i][3]
					alt[i] = curr_line[i][4].replace('.','').split(',')
					if len(alt[i]) == 1 and len(alt[i][0]) == 0:
						alt[i] = []
					allele[i] = [ref[i]] + alt[i] + [None]
					try:
						ft_idx[i] = curr_line[i][8].split(':').index('FT')
					except ValueError:
						ft_idx[i] = None
					
					gt_idx[i] = curr_line[i][8].split(':').index('GT')
				except StopIteration:
					at_eof[i] = 1
					to_advance[i] = 0
					curr_pos[i] = ('',-1)
					ref[i] = None
					alt[i] = None
					gt_idx[i] = None
					ft_idx[i] = None
							
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
		
		# only perform concordance check if called in both sites!
		if all(to_advance):
			is_snp = all( (all( (len(s)==1 for s in allele_list[:-1]) ) for allele_list in allele ) )
#			print is_snp, allele
	
			filter_status = [(v and curr_line[i][6] == "PASS") for i, v in enumerate(to_advance)]
			# get a mapping of id : set of alleles
			# We'll assume everyone here is diploid!
			
#			print len(curr_line)
			
			r_allele_set = {p : [NoneSet for v in vcf_list] for p in common_ids}
			f_allele_set = {p : [NoneSet for v in vcf_list] for p in common_ids}
			
			for i, g in enumerate(curr_line):
#				print len(g), g

				for p in common_ids:
					# geno will be either 0,1,2 (normal encoding), or -1 (missing)
#					if p == '38913308':
#						print i, vcf_ids[i][p], g[vcf_ids[i][p]], g[vcf_ids[i][p]].split(':')[gt_idx[i]].replace('.','-1').split('/')
					
					geno_set = set(g[vcf_ids[i][p]].split(':')[gt_idx[i]].replace('.','-1').split('/'))
					passed = filter_status[i] and (ft_idx[i] is None or g[vcf_ids[i][p]].split(':')[ft_idx[i]] == "PASS")
					
					r_allele_set[p][i] = set(allele[i][int(g_idx)] for g_idx in geno_set)
					f_allele_set[p][i] = (NoneSet, r_allele_set[p][i])[passed]
				
			
			# Now check the allele_set
			for p in common_ids:
#				print curr_pos, p, r_allele_set[p]
			
				base_set = set.intersection(*r_allele_set[p])
				all_set = set.union(*r_allele_set[p])
				
#				print base_set, all_set
					
				# Uh-oh! there's a discordance!	
				if base_set != all_set:
					no_miss_set = all_set - NoneSet
					
					# in this case, we have missing discord
					if len(no_miss_set) == 1:
						r_concord[-1] += 1
						r_snp_concord[-1] += is_snp
						
						# otherwise, let's check for homozygous discord
					elif sum(len(s-NoneSet)==1 for s in r_allele_set[p]) > 1:
						r_concord[2] += 1
						r_snp_concord[2] += is_snp
					else:
						r_concord[1] += 1
						r_snp_concord[1] += is_snp
				else:
					r_concord[0] += 1
					r_snp_concord[0] += is_snp
					
					
				base_set = set.intersection(*f_allele_set[p])
				all_set = set.union(*f_allele_set[p])
					
				# Uh-oh! there's a discordance!	
				if base_set != all_set:
					no_miss_set = all_set - NoneSet
					
#					print working_pos, is_snp, p, base_set, all_set, f_allele_set[p]
					# in this case, we have missing discord
					if sum(len(s-NoneSet)==1 for s in f_allele_set[p]) > 1:
						f_concord[2] += 1
						f_snp_concord[2] += is_snp
					elif sum(len(s-NoneSet)>1 for s in f_allele_set[p]) > 1:
						print working_pos, is_snp, p, base_set, all_set, f_allele_set[p]
						f_concord[1] += 1
						f_snp_concord[1] += is_snp
					
						# otherwise, let's check for homozygous discord
					else:
						
						f_concord[-1] += 1
						f_snp_concord[-1] += is_snp
				else:
					f_concord[0] += 1
					f_snp_concord[0] += is_snp
						
	print "Raw Results"
	print "==========="
	print "Total:"
	print_concord(r_concord)
	print "\nSNPs Only:"
	print_concord(r_snp_concord)

	print "\n\nFiltered Results"
	print "==========="
	print "Total:"
	print_concord(f_concord)
	print "\nSNPs Only:"
	print_concord(f_snp_concord)
	
					
	exit(0)
	
	r_site_count = {}
	f_site_count = {}
	
	r_referent_count = {}
	f_referent_count = {}
	
	r_polyallelic_count = {}
	f_polyallelic_count = {}
	
	r_alt_disagree = 0
	f_alt_disagree = 0
	r_alt_discord = 0
	f_alt_discord = 0
	
	r_single_concord = [0, 0, 0, 0]
	f_single_concord = [0, 0, 0, 0]
	
	r_all_concord = [0, 0, 0, 0]
	f_all_concord = [0, 0, 0, 0]

	
	r_comp_sites = 0
	f_comp_sites = 0
	r_comp_extra = 0
	f_comp_extra = 0
	r_single = 0
	f_single = 0
	
	nonref = set((frozenset([i]) for i in ['.','']))
	
	while True:
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
		
		r_alt_list = [(len(s), s) for s in r_alts]
		r_alt_list.sort(reverse=True)
		
		f_alt_list = [(len(s), s) for s in f_alts]
		f_alt_list.sort(reverse=True)

		
		# check for a superset of alleles
		r_alt_discord += (sum((len(v[1] - r_alt_list[0][1]) for v in r_alt_list[1:])) > 0)
		f_alt_discord += (sum((len(v[1] - f_alt_list[0][1]) for v in f_alt_list[1:])) > 0)
		
		r_alt_disagree += (len(r_alts) > 1)
		f_alt_disagree += (len(f_alts) > 1)
		
		

			
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

						# NOTE: do we want to treat filtered as "missing" or "not present"?
						#else:
						#	f_geno_count[p][-1] += 1
				
			r_curr_concord = [sum((sum((n*(n-1)/2 for n in v)) for v in r_geno_count.itervalues())),
							  sum((v[0]*v[1] + v[1]*v[2] for v in r_geno_count.itervalues())),
							  sum((v[0]*v[2] for v in r_geno_count.itervalues())),
							  sum((v[-1]*(v[0] + v[1] + v[2]) for v in r_geno_count.itervalues()))]
#			print working_pos  
#			print r_curr_concord
			
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
		
	
	
	print "Sites by file:"
	
	print_dict(r_site_count)
	
	print "Referent sites by file:"
	
	print_dict(r_referent_count)

	print "Polyallelic sites by file:"
	
	print_dict(r_polyallelic_count)
	
	print "Number of disagreeing ALT:", r_alt_disagree
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
	
	print "Number of disagreeing ALT:", f_alt_disagree
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
	
