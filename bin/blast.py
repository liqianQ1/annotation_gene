#!/usr/bin/env python
#coding=utf-8

import os
import argparse
import ConfigParser
parser=argparse.ArgumentParser(prog='Blast',description='blast')
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
	blast=open(rootdir+'/03_blast/Blast_'+eachS+'.sh','w')
#read1	
	blast.write('/share/public/software/seqkit/seqkit sample -n 10000 -s 11 %s -o %s/03_blast/%s.blast.fastq\n\n' %(dict1[eachS],rootdir,eachS))
	blast.write('perl %s/bin/fq2fa.pl %s.blast.fastq %s.blast.fasta\n\n'%(pipline,eachS,eachS))
	blast.write('/share/public/software/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastn -db /share/service/denovo/test/NT/20170822/nt  -outfmt  11  -num_threads  32  -evalue  1e-5  -num_alignments  10  -max_hsps  1  -query %s.blast.fasta  -out  %s.blast.outfmt11\n\n'%(eachS,eachS))
	blast.write('/share/public/software/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blast_formatter -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname scomname stitle sblastname sskingdom" -archive %s.blast.outfmt11  -num_alignments 1 > %s.blast.outfmt11to6\n\n'%(eachS,eachS))
	blast.write('perl %s/bin/query_tax.pl  %s.blast.outfmt11to6 %s\n\n'%(pipline,eachS,eachS))
	blast.close()	
