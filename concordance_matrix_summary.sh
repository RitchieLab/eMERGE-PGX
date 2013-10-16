#!/bin/bash

echo -e "Concord\tHetero.\tHomo.\tMissing"
cat $1 | sed '1d' | cut -f2- | sed 's/\t/\n/g' | awk '{if (NR %5 == 1) good += $1; else if (NR % 4 == 0 || NR > 12) ugly += $1; else if (NR % 3 == 0) vbad += $1; else bad += $1  } END {sum=good+vbad+bad+ugly; printf "%.3f\t%.3f\t%.3f\t%.3f\n", good/sum*100, bad/sum*100, vbad/sum*100, ugly/sum*100}'
