#!/usr/bin/env python
#coding=utf-8

import os
import argparse
import ConfigParser
parser=argparse.ArgumentParser(prog='Fastqc',description='fastqc')
parser.add_argument('-RD','--rootdir',help='the directory of the project',metavar='')
parser.add_argument('-PP','--pipeline',help='the directory of pipeline',metavar='')
argv=vars(parser.parse_args())

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
####raw=open(rootdir+'/'+'sample_info','r')
####for eachline in raw:
####	eachline=eachline.strip().split('\t')	
####	samplelist.append(eachline[0])
####	dict1[eachline[0]]=eachline[1]
####	dict2[eachline[0]]=eachline[2]

#### samplelists=list(set(samplelist))
#####samplelists.sort(key=samplelist.index)
#####raw.close()

#for eachS in samplelists:
#	fastqc=open(rootdir+'/02_hisat2/hisat2_'+eachS+'.sh','w')
##hisat2	
#	fastqc.write('%s/software/hisat2 -x ../01_index/ref -p 60 -1 %s -2 %s -S %s.sam\n\n' %(pipeline,dict1[eachS],dict2[eachS],eachS)) 
##samtobam
#	fastqc.write('%s/software/samtools view -@ 40 -bS %s.sam > %s.bam \n\n' %(pipeline,eachS,eachS))	
#	fastqc.close()	

fastqc_merge=open(rootdir+'/03_stringtie/Merge_hisat2.sh','w')
fastqc_merge.write('#!/bin/bash\n')
fastqc_merge.write('/home/biohuaxing/Software/samtools-1.19.2/bin/samtools merge %s/03_stringtie/merge.bam %s/02_hisat2/*.bam \n\n' %(rootdir,rootdir))
fastqc_merge.write('/home/biohuaxing/Software/samtools-1.19.2/bin/samtools sort -@ 8  %s/03_stringtie/merge.bam -o %s/03_stringtie/merge.sort.bam \n\n' %(rootdir,rootdir))
fastqc_merge.write('/home/biohuaxing/Software/samtools-1.19.2/bin/samtools index %s/03_stringtie/merge.sort.bam \n\n' %(rootdir))
fastqc_merge.write('/home/biohuaxing/Software/samtools-1.19.2/bin/samtools flagstat %s/03_stringtie/merge.sort.bam > %s/03_stringtie/merge.flagstat.txt \n\n' %(rootdir,rootdir))
fastqc_merge.close()
