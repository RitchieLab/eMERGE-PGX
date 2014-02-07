#!/usr/bin/env python

import csv
import sys

if __name__ == "__main__":
	raw_f = file(sys.argv[1],'rb')
	csv_f = csv.DictReader(raw_f)
	
	# Read the summary missing rate
	miss_f = file(sys.argv[2], 'rb')
	miss_csv = csv.DictReader(miss_f)
	
	miss_dict = {l["SMtag"] : l["MissRate"] for l in miss_csv}
	miss_f.close()
	
	# Get the summary data we need - here we're just storing sums
	# - Proportion of target @ 0 / min. coverage
	# - # SNVs on target
	# - Total Reads
	# - Mapped Reads
	# - Total Read Bases
	# - Usable Read Bases
	# - On-Target Read Bases
	# - On-Target Usable Bases
	# - Proportion On-Target
	# - SNVs - % SNVs in dbSNP37
	# - # Var. per unit on-target bases
	# - Novel SNVs
	sdata = [0] * 12
	ontarget_base = []
	
	print '\t'.join(["Local_ID","Cov>10x", "Cov>30x", "Mean_Cov", "Mean_Missing", "Num_SNV_On_Target"])
	
	n = 0
	for f in csv_f:
		print '\t'.join([f["Local_Id"],f["PCT_TARGET_BASES_10X"],f["PCT_TARGET_BASES_30X"],f["MEAN_TARGET_COVERAGE"],miss_dict[f["SM_tag"]],f["COUNT_SNV_ON_TARGET"]])
		n+=1
		
		sdata[0] += float(f["ZERO_CVG_TARGETS_PCT"])
		sdata[1] += int(f["COUNT_SNV_ON_TARGET"])
		sdata[2] += int(f["TOTAL_READS"])
		sdata[3] += int(f["TOTAL_READS"]) - int(f["UNMAPPED_READS"])
		sdata[4] += int(f["TOTAL_READS"]) * float(f["MEAN_READ_LENGTH_PAIR"])
		sdata[5] += int(f["PF_UQ_READS_ALIGNED"]) * float(f["MEAN_READ_LENGTH_PAIR"])
		sdata[6] += int(f["ON_TARGET_BASES"])
		sdata[7] += int(f["ON_TARGET_BASES"]) * float(f["PCT_USABLE_BASES_ON_TARGET"])
		sdata[8] += float(f["PCT_SELECTED_BASES"])
		sdata[9] += float(f["PERCENT_SNV_ON_TARGET_SNP137"])
		sdata[10] += int(f["COUNT_SNV_ON_TARGET"]) / float(f["TARGET_TERRITORY"])
		sdata[11] += int(f["COUNT_SNV_ON_TARGET"]) * (1 - float(f["PERCENT_SNV_ON_TARGET_SNP137"])/100)
		
		ontarget_base.append(int(f["ON_TARGET_BASES"]) * float(f["PCT_USABLE_BASES_ON_TARGET"]))
		
	raw_f.close()
	
	def median(lst):
		return sum(sorted(lst)[(len(lst)-1)/2:(len(lst)-1)/2+(len(lst)-1)%2+1])/float((len(lst)-1)%2 + 1)
	
	print >> sys.stderr, "Summary level statistics:"
	print >> sys.stderr, "Proportion of target @ 0 / min. coverage:", sdata[0] / n * 100
	print >> sys.stderr, "Mean Missing Call Rate:", sum((float(v) for v in miss_dict.itervalues())) / len(miss_dict)
	print >> sys.stderr, "# SNVs on target:", sdata[1] / float(n)
	print >> sys.stderr, "Total Reads:", sdata[2] / float(n)
	print >> sys.stderr, "Mapped Reads:", sdata[3] / float(n)
	print >> sys.stderr, "Total Read Bases:", sdata[4] / float(n)
	print >> sys.stderr, "Usable Read Bases:", sdata[5] / float(n)
	print >> sys.stderr, "On-Target Read Bases:", sdata[6] / float(n)
#	print >> sys.stderr, "On-Target Usable Bases (Total):", sum(ontarget_base)
	print >> sys.stderr, "On-Target Usable Bases (Mean):", sum(ontarget_base) / float(len(ontarget_base))
#	print >> sys.stderr, "On-Target Usable Bases (Median):", median(ontarget_base)
	print >> sys.stderr, "Proportion On-Target:", sdata[8] / n
	print >> sys.stderr, "SNVs - % SNVs in dbSNP37:", sdata[9] / n
	print >> sys.stderr, "Novel SNVs:", int(round(sdata[11],0)) # sdata[1] * (1 - sdata[9] / (100 * n))
	print >> sys.stderr, "Var. per unit on-target bases:", sdata[10] / n


	
