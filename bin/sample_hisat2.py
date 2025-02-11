#!/usr/bin/env python
#coding=utf-8

import os
import argparse
import ConfigParser
parser=argparse.ArgumentParser(prog='Fastqc',description='fastqc')
parser.add_argument('-S','--sampleinfo',help='the samplelist of the project',metavar='')
parser.add_argument('-RD','--rootdir',help='the directory of the project',metavar='')
parser.add_argument('-PP','--pipeline',help='the directory of pipeline',metavar='')
argv=vars(parser.parse_args())
if	argv['sampleinfo'] == None:
	raise Exception('You should provide the samplelist of the project!')
if  argv['sampleinfo'] != None:
	sampleinfo=argv['sampleinfo'].strip()
if argv['pipeline'] == None:
	raise Exception('You should provide the directory of pipeline!')
if argv['pipeline'] != None:
	pipeline=argv['pipeline'].strip()
if argv['rootdir'] == None:
	raise Exception('You should provide the directory of the project!')
if argv['rootdir'] != None:
	rootdir=argv['rootdir'].strip()
argv=vars(parser.parse_args())

samplelist=[]
dict1={}
dict2={}

#read raw_fastq
#raw=open(rootdir+'/'+'sample_info','r')
raw=open(sampleinfo,'r')
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
	fastqc=open(rootdir+'/02_hisat2/hisat2_'+eachS+'.sh','w')
#hisat2 #!/bin/bash	
#	fastqc.write('#!/bin/bash\n') 
	fastqc.write('/home/biohuaxing/Personal/liqian/Anaconda3/envs/ann/bin/hisat2 -x %s/01_index/ref -p 4 -1 %s -2 %s --pen-noncansplice 1000000 |/home/biohuaxing/Software/samtools-1.19.2/bin/samtools view -@ 4 -bS > %s/02_hisat2/%s.bam\n' %(rootdir,dict1[eachS],dict2[eachS],rootdir,eachS))
#	fastqc.write('%s/software/samtools sort -@ 8 %s/02_hisat2/%s.bam -o %s/02_hisat2/%s.bam.sort.bam \n' %(pipeline,rootdir,eachS,rootdir,eachS)) 
	fastqc.write('touch %s/02_hisat2/hisat2_%s\n' %(rootdir,eachS)) 
	fastqc.close()	
