#!/usr/bin/env python
#coding=utf-8

import os
import argparse
import ConfigParser
parser=argparse.ArgumentParser(prog='Fastqc',description='fastqc')
parser.add_argument('-RD','--rootdir',help='the directory of the project',metavar='')
parser.add_argument('-PP','--pipline',help='the directory of pipline',metavar='')
argv=vars(parser.parse_args())

if argv['pipline'] == None:
	pipline='/share/public/pipeline/SS/denovo/Hic_Pro/1.2.1'
if argv['pipline'] != None:
	pipline=argv['pipline'].strip()
if argv['rootdir'] == None:
	raise Exception('You should provide the directory of the project!')
if argv['rootdir'] != None:
	rootdir=argv['rootdir'].strip()
argv=vars(parser.parse_args())

samplelist=[]
dict1={}
dict2={}

#read raw_fastq
raw=open(rootdir+"/config/sample_info",'r')
for eachline in raw:
	eachline=eachline.strip().split('\t')	
	samplelist.append(eachline[0])
	dict1[eachline[0]]=eachline[1]
	dict2[eachline[0]]=eachline[2]

samplelists=list(set(samplelist))
#print samplelists
samplelists.sort(key=samplelist.index)
raw.close()

for eachS in samplelists:
	fastqc=open(rootdir+'/00_qc/Fastqc_'+eachS+'.sh','w')
#read1	
	fastqc.write('%s/bin/BerryTools FqCheck -f %s -o %s/00_qc/%s_1.check\n\n' %(pipline,dict1[eachS],rootdir,eachS)) 
#read2
	fastqc.write('%s/bin/BerryTools FqCheck -f %s -o %s/00_qc/%s_2.check\n\n' %(pipline,dict2[eachS],rootdir,eachS))	
	fastqc.write('%s/bin/plot_base_qual.pl -1 %s/00_qc/%s_1.check -2 %s/00_qc/%s_2.check -o %s/00_qc/%s\n\n' %(pipline,rootdir,eachS,rootdir,eachS,rootdir,eachS))
	fastqc.close()	

fastqc_merge=open(rootdir+'/00_qc/Merge_fastqc.sh','w')
fastqc_merge.write('%s/bin/ReadsIn.pl -i %s/00_qc' %(pipline,rootdir))
fastqc_merge.close()
