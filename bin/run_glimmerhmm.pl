#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

die "perl $0 <glimmer> <fasta> <db> <outgff> <outdir>\n" unless(@ARGV == 5);

my ($glimmer, $fasta, $db, $gff, $odir) = @ARGV;
my $inseq  = Bio::SeqIO->new( -file   => "<$fasta", -format => 'fasta');
system("echo -n '' >$gff");
while(my $seq_obj = $inseq->next_seq){
	open OUT, ">$odir/temp.fa" or die "$!\n";
	print OUT ">".$seq_obj->id."\n";
	print OUT $seq_obj->seq."\n";
	close OUT;
	!system("$glimmer $odir/temp.fa -d $db -g  -f  -o temp.gff && cat temp.gff >>$gff && rm temp.gff && rm $odir/temp.fa") or die "Run glimmer error!\n";
}

