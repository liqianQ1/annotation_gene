#! /usr/bin/env python
# -*- coding:UTF-8 -*-

__author__ = 'pangshuai' 

import sys
import os

finput1 = open(sys.argv[1],"r")
finput2 = open(sys.argv[2],"r")
finput3 = open(sys.argv[3],"r")
foutput = open(sys.argv[4],"w")

for line in finput1:
	if line.startswith("Total_pairs_processed"):
		Total_pairs_processed = float(line.strip().split("\t")[1])
	if line.startswith("Unique_paired_alignments"):
		Unique_paired_alignments = float(line.strip().split("\t")[1])

for line in finput2:
	if line.strip().split("\t")[0]=="valid_interaction":
		valid_interaction = int(line.strip().split("\t")[1])
		valid_interaction_total_percent = round((valid_interaction*100/Total_pairs_processed),2)
		valid_interaction_uniq_percent = round((valid_interaction*100/Unique_paired_alignments),2)
	if line.startswith("valid_interaction_rmdup"):
		valid_interaction_rmdup = int(line.strip().split("\t")[1])
		valid_interaction_rmdup_total_percent = round((valid_interaction_rmdup*100/Total_pairs_processed),2)
		valid_interaction_rmdup_uniq_percent = round((valid_interaction_rmdup*100/Unique_paired_alignments),2)
		valid_interaction_dup = valid_interaction-valid_interaction_rmdup
		valid_interaction_dup_total_percent = round((valid_interaction_dup*100/Total_pairs_processed),2)
		valid_interaction_dup_uniq_percent = round((valid_interaction_dup*100/Unique_paired_alignments),2)

for line in finput3:
	if line.startswith("Dangling_end_pairs"):
		Dangling_end_pairs = int(line.strip().split("\t")[1])
		Dangling_end_pairs_total_percent = round((Dangling_end_pairs*100/Total_pairs_processed),2)
		Dangling_end_pairs_uniq_percent = round((Dangling_end_pairs*100/Unique_paired_alignments),2)
	if line.startswith("Religation_pairs"):
		Religation_pairs = int(line.strip().split("\t")[1])
		Religation_pairs_total_percent = round((Religation_pairs*100/Total_pairs_processed),2)
		Religation_pairs_uniq_percent = round((Religation_pairs*100/Unique_paired_alignments),2)
	if line.startswith("Self_Cycle_pairs"):
		Self_Cycle_pairs = int(line.strip().split("\t")[1])
		Self_Cycle_pairs_total_percent = round((Self_Cycle_pairs*100/Total_pairs_processed),2)
		Self_Cycle_pairs_uniq_percent = round((Self_Cycle_pairs*100/Unique_paired_alignments),2)
	if line.startswith("Single-end_pairs"):
		Single_end_pairs = int(line.strip().split("\t")[1])
		Single_end_pairs_total_percent = round((Single_end_pairs*100/Total_pairs_processed),2)
		Single_end_pairs_uniq_percent = round((Single_end_pairs*100/Unique_paired_alignments),2)
	if line.startswith("Dumped_pairs"):
		Dumped_pairs = int(line.strip().split("\t")[1])
		Dumped_pairs_total_percent = round((Dumped_pairs*100/Total_pairs_processed),2)
		Dumped_pairs_uniq_percent = round((Single_end_pairs*100/Unique_paired_alignments),2)

foutput.write("Type"+"\t"+"Pairs Count"+"\t"+"Percent of total(%)"+"\t"+"Percent of unique(%)"+"\n")
foutput.write("Valid interaction pairs"+"\t"+str(valid_interaction)+"\t"+str(valid_interaction_total_percent)+"\t"+str(valid_interaction_uniq_percent)+"\n")
foutput.write("Valid interaction rmdup"+"\t"+ str(valid_interaction_rmdup)+"\t"+str(valid_interaction_rmdup_total_percent)+"\t"+str(valid_interaction_rmdup_uniq_percent)+"\n")
foutput.write("Valid interaction duplicate"+"\t"+ str(valid_interaction_dup)+"\t"+str(valid_interaction_dup_total_percent)+"\t"+str(valid_interaction_dup_uniq_percent)+"\n")
foutput.write("Dangling end pairs"+"\t"+str(Dangling_end_pairs)+"\t"+str(Dangling_end_pairs_total_percent)+"\t"+str(Dangling_end_pairs_uniq_percent)+"\n")
foutput.write("Religation pairs"+"\t"+str(Religation_pairs)+"\t"+str(Religation_pairs_total_percent)+"\t"+str(Religation_pairs_uniq_percent)+"\n")
foutput.write("Self Cycle pairs"+"\t"+str(Self_Cycle_pairs)+"\t"+str(Self_Cycle_pairs_total_percent)+"\t"+str(Self_Cycle_pairs_uniq_percent)+"\n")
foutput.write("Single-end pairs"+"\t"+str(Single_end_pairs)+"\t"+str(Single_end_pairs_total_percent)+"\t"+str(Single_end_pairs_uniq_percent)+"\n")
foutput.write("Dumped pairs"+"\t"+str(Dumped_pairs)+"\t"+str(Dumped_pairs_total_percent)+"\t"+str(Dumped_pairs_uniq_percent)+"\n")


finput1.close()
finput2.close()
finput3.close()
foutput.close()
