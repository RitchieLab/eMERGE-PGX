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
	print "Total Concordance:", "%.3f" % (r[0] + r[3] + r[4] + r[5])
	print "  Ref/Ref Concordance:", "%.3f" % r[3]
	print "  NonRef/NonRef Concordance:", "%.3f" % (r[0] + r[4] + r[5])
	print "    Het/Het Concordance:", "%.3f" % r[4]
	print "    Missing Concordance:", "%.3f" % r[5]

	print "Total Discord:", "%.3f" % (r[1] + r[2])
	print "  Heterozygous Discord:", "%.3f" % r[1]
	print "  Homozygous Discord:", "%.3f" % r[2]
	print "Missing Discord:", "%.3f" % r[-1]

def print_samp_concord(samp_id, c):
	r = [100*v / float(sum(c)) for v in c]
	print samp_id, ":", "%.3f / %.3f" % (r[1], r[-1])

def get_MAF(x):
	"""
	Gets the MAF given the info string
	"""
	for a in x.split(';'):
		if a.startswith("MLEAF"):
			key, val = a.split("=")
			afs = [float(v) for v in val.split(',')]
			afs.append(1-sum(afs))
			afs.sort(reverse=True)
			return afs[1]
	
	return 0
	
def to_allele_str(x):
	if len(x) == 1:
		if None in x:
			return './.'
		else:
			s = x.pop()
			x.add(s)
			return s + '/' + s
	else:
		s1 = x.pop()
		s2 = x.pop()
		x.add(s1)
		x.add(s2)
		return s1 + '/' + s2
		

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
	is_variant = [0 for f in vcf_list]
	curr_pos = [('0',0) for f in vcf_list]
	curr_line = ['' for f in vcf_list]
	at_eof = [0 for f in vcf_list]
	oparen = ['(' for f in vcf_list]
	cloparen = [')' for f in vcf_list]
	
	prev_pos = ('0',0)
	seen_chrom = set()
	
	ref = [None for f in vcf_list]
	alt = [None for f in vcf_list]
	allele = [None for f in vcf_list]
	maf = [0 for f in vcf_list]
	
	gt_idx = [None for f in vcf_list]
	ft_idx = [None for f in vcf_list]

	# concordance is:
	# Hom. Alt concordance
	# Het. Discord
	# Hom. Discord
	# Ref/Ref Concord
	# Het. Concord
	# Missing Concord
	# missing discord
	r_concord = [0, 0, 0, 0, 0, 0, 0]
	f_concord = [0, 0, 0, 0, 0, 0, 0]
	
	r_snp_concord = [0, 0, 0, 0, 0, 0, 0]
	f_snp_concord = [0, 0, 0, 0, 0, 0, 0]
	
	r_novel_concord = [0, 0, 0, 0, 0, 0, 0]
	f_novel_concord = [0, 0, 0, 0, 0, 0, 0]
	
	r_snp_novel_concord = [0, 0, 0, 0, 0, 0, 0]
	f_snp_novel_concord = [0, 0, 0, 0, 0, 0, 0]
	
	id_concord = {p : [[0,0,0],[0,0,0]] for p in common_ids}
	id_snp_concord = {p : [[0,0,0],[0,0,0]] for p in common_ids}
	
	discord_details = []
	
	while not all(at_eof):
		for i,f in enumerate(vcf_list):
			if to_advance[i]:
				try:
					curr_line[i] = f.next().strip().split('\t')
					curr_pos[i] = (curr_line[i][0], int(curr_line[i][1]))
					ref[i] = curr_line[i][3]
					alt[i] = curr_line[i][4].replace('.','').split(',')
					maf[i] = get_MAF(curr_line[i][7])
					if len(alt[i]) == 1 and len(alt[i][0]) == 0:
						alt[i] = []
					allele[i] = [ref[i]] + alt[i] + [None]
					is_variant[i] = (len(alt[i]) > 0)
					try:
						ft_idx[i] = curr_line[i][8].split(':').index('FT')
					except ValueError:
						ft_idx[i] = None
					
					gt_idx[i] = curr_line[i][8].split(':').index('GT')
				except StopIteration:
					at_eof[i] = 1
					to_advance[i] = 0
					is_variant[i] = 0
					maf[i] = 0
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
		
		# only perform concordance check if called in both sites and at least one site has a variant!
		if all(to_advance) and any(is_variant):
		
			is_novel = all( (curr_line[i][2] == "." for i,v in enumerate(to_advance)) )
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
					
					geno_set = set(g[vcf_ids[i][p]].split(':')[gt_idx[i]].replace('.','-1').replace('|', '/').split('/'))
					passed = filter_status[i] and (ft_idx[i] is None or g[vcf_ids[i][p]].split(':')[ft_idx[i]] == "PASS")
					
					r_allele_set[p][i] = set(allele[i][int(g_idx)] for g_idx in geno_set)
					f_allele_set[p][i] = (NoneSet, r_allele_set[p][i])[passed]
				
			
			# Now check the allele_set
						
			for p in common_ids:
#				print curr_pos, p, r_allele_set[p]
			
				base_set = set.intersection(*r_allele_set[p])
				all_set = set.union(*r_allele_set[p])
				
#				print base_set, all_set

				r_discord = True
				f_discord = True
					
				# Uh-oh! there's a discordance!	
				if base_set != all_set:
					no_miss_set = all_set - NoneSet
					
					
					# we know that it's heterozygous if everybody agrees on one allele!
					if len(base_set) == 1:
						r_concord[1] += 1
						r_snp_concord[1] += is_snp
						r_novel_concord[1] += is_novel
						r_snp_novel_concord[1] += is_snp and is_novel
						id_concord[p][0][1] += 1
						id_snp_concord[p][0][1] += is_snp						
					
						# otherwise, let's check for homozygous discord
					elif sum((len(s-NoneSet) - len(s-base_set)) == 0 for s in r_allele_set[p]) > 1:
						#if(is_snp):
						#	print "homo. raw concordance"
						#	print base_set
						#	print r_allele_set[p]
						#	print [len(s-NoneSet) - len(s-base_set) for s in r_allele_set[p]]
						#	print all_set
						#	print no_miss_set
						#	print base_set
						#	sys.exit(1)
						r_concord[2] += 1
						r_snp_concord[2] += is_snp
						r_novel_concord[2] += is_novel
						r_snp_novel_concord[2] += is_snp and is_novel
						id_concord[p][0][1] += 1
						id_snp_concord[p][0][1] += is_snp						
					else:
						#if(is_snp):
						#	print "missing raw concordance"
						#	print base_set
						#	print r_allele_set[p]
						#	print [len(s) - len(s-base_set-NoneSet) for s in r_allele_set[p]]
						#	print all_set
						#	print no_miss_set
						#	print base_set
						#	sys.exit(1)
						r_concord[-1] += 1
						r_snp_concord[-1] += is_snp
						r_novel_concord[-1] += is_novel
						r_snp_novel_concord[-1] += is_snp and is_novel
						id_concord[p][0][-1] += 1
						id_snp_concord[p][0][-1] += is_snp						
				else:
					# OK, so we know it's concordant - let's check if it's ref/ref
					r_discord = False
					isref = (len(base_set) == 1 and len(base_set - set(ref)) == 0)
					ishet = (len(base_set) > 1)
					ismiss = (len(base_set) == 1 and len(base_set - NoneSet) == 0)
					r_concord[0 + 3*isref + 4*ishet + 5*ismiss] += 1
					r_snp_concord[0 + 3*isref + 4*ishet + 5*ismiss] += is_snp
					r_novel_concord[0 + 3*isref + 4*ishet + 5*ismiss] += is_novel
					r_snp_novel_concord[0 + 3*isref + 4*ishet + 5*ismiss] += is_snp and is_novel
					id_concord[p][0][0] += 1
					id_snp_concord[p][0][0] += is_snp						
					
					
				base_set = set.intersection(*f_allele_set[p])
				all_set = set.union(*f_allele_set[p])
					
				# Uh-oh! there's a discordance!	
				if base_set != all_set:
					no_miss_set = all_set - NoneSet
					
#					print working_pos, is_snp, p, base_set, all_set, f_allele_set[p]
					if len(base_set) == 1:
#						if(is_snp):
#							print working_pos, is_snp, p, base_set, all_set, f_allele_set[p]
						f_concord[1] += 1
						f_snp_concord[1] += is_snp
						f_novel_concord[1] += is_novel
						f_snp_novel_concord[1] += is_snp and is_novel
						id_concord[p][1][1] += 1
						id_snp_concord[p][1][1] += is_snp						
					elif sum((len(s-NoneSet) - len(s-base_set)) == 0 for s in f_allele_set[p]) > 1:
											#if(is_snp):
						#	print "homo. raw concordance"
						#	print base_set
						#	print r_allele_set[p]
						#	print [len(s-NoneSet) - len(s-base_set) for s in r_allele_set[p]]
						#	print all_set
						#	print no_miss_set
						#	print base_set
						#	sys.exit(1)
						#if(is_snp):
						#	print working_pos, is_snp, p, base_set, all_set, f_allele_set[p]
						f_concord[2] += 1
						f_snp_concord[2] += is_snp
						f_novel_concord[2] += is_novel
						f_snp_novel_concord[2] += is_snp and is_novel
						id_concord[p][1][1] += 1
						id_snp_concord[p][1][1] += is_snp						
					else:
						f_concord[-1] += 1
						f_snp_concord[-1] += is_snp
						f_novel_concord[-1] += is_novel
						f_snp_novel_concord[-1] += is_snp and is_novel
						id_concord[p][1][-1] += 1
						id_snp_concord[p][1][-1] += is_snp						
				else:
					# Again, check for ref/ref concordance
					f_discord = False
					isref = (len(base_set) == 1 and len(base_set - set(ref)) == 0)
					ishet = (len(base_set) > 1)
					ismiss = (len(base_set) == 1 and len(base_set - NoneSet) == 0)					
					f_concord[0 + isref * 3+ 4*ishet + 5*ismiss] += 1
					f_snp_concord[0 + isref * 3+ 4*ishet + 5*ismiss] += is_snp
					f_novel_concord[0 + isref * 3+ 4*ishet + 5*ismiss] += is_novel
					f_snp_novel_concord[0 + isref * 3+ 4*ishet + 5*ismiss] += is_snp and is_novel
					id_concord[p][1][0] += 1
					id_snp_concord[p][1][0] += is_snp
					
				
				# OK, if there was a discordance, print everything we know!
				if r_discord or f_discord:
					raw_alleles = [to_allele_str(al_s) for al_s in r_allele_set[p]]
					filt_allele = [to_allele_str(al_s) for al_s in f_allele_set[p]]
					
					allele_str = ' '.join(':'.join(s) for s in zip(raw_alleles, filt_allele)) 
					pos_str = working_pos[0] + ":" + str(working_pos[1])
					allele_type = ("INDEL","SNP")[is_snp]
					novel_type = ("dbSNP", "NOVEL")[is_novel]
					
					to_add = " ".join([pos_str, p, allele_type, novel_type, allele_str, " ".join((str(m) for m in maf))])
					discord_details.append(to_add)
		
					
	print "Raw Results"
	print "==========="
	print "Total:"
	print_concord(r_concord)
	print "\nSNPs Only:"
	print_concord(r_snp_concord)
	print "\nNovel:"
	print_concord(r_novel_concord)
	print "\nNovel SNPs Only:"
	print_concord(r_snp_novel_concord)

	print "\n\nFiltered Results"
	print "==========="
	print "Total:"
	print_concord(f_concord)
	print "\nSNPs Only:"
	print_concord(f_snp_concord)
	print "\nNovel:"
	print_concord(f_novel_concord)
	print "\nNovel SNPs Only:"
	print_concord(f_snp_novel_concord)
	
	print "\n\nBy Sample(Raw, All):"
	print "============="
	for p in common_ids:
		print_samp_concord(p, id_concord[p][0]);

	print "\n\nBy Sample(Raw, SNP):"
	print "============="
	for p in common_ids:
		print_samp_concord(p, id_snp_concord[p][0]);
		
	print "\n\nBy Sample(Filtered, All):"
	print "============="
	for p in common_ids:
		print_samp_concord(p, id_concord[p][1]);
		
	print "\n\nBy Sample(Filtered, SNP):"
	print "============="
	for p in common_ids:
		print_samp_concord(p, id_snp_concord[p][1]);				
	
	
	print "\n\nRaw Concordance details"
	print "--------------------"
	print "chr:pos ID Type Novelty Allele1(raw:filt) Allele2(raw:filt) MAF1 MAF2"
	
	print '\n'.join(discord_details)
