#! /usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = 'pangshuai' 

import sys
import os

finput = open(sys.argv[1],"r")
fout = open(sys.argv[2],"w")

line_num = 0
for line in finput:
	if line.startswith("#"):
		continue
	if line.startswith("Total_pairs_processed"):
		Total_pairs_processed = line.strip("_")
	if line.startswith("Unmapped_pairs"):
		Unmapped_pairs = line.strip("_")
	if line.startswith("Low_qual_pairs"):
		Low_qual_pairs = line.strip("_")
	if line.startswith("Unique_paired_alignments"):
		Unique_paired_alignments = line.strip("_")
	if line.startswith("Multiple_pairs_alignments"):
		Multiple_pairs_alignments = line.strip("_")
	if line.startswith("Pairs_with_singleton"):
		Pairs_with_singleton = line.strip("_")
	
fout.write("Type"+"\t"+"Pairs Count"+"\t"+"Percent(%)"+"\n"+Total_pairs_processed+Unmapped_pairs+Low_qual_pairs+Unique_paired_alignments+Multiple_pairs_alignments+Pairs_with_singleton)
	
finput.close()
fout.close()
