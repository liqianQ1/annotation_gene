#!/home/biohuaxing/Personal/liqian/Anaconda3/envs/ann/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../lib/";
use Fasta;
use ANNOTATION;
use Data::Dumper;

my $usage=<<USAGE;
	Usage:	$0 <input.gff>
USAGE

#print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

die $usage unless @ARGV;
my (@iput) = @ARGV;
# 处理每个输入文件

my %gff;
foreach (@iput) {
	ANNOTATION::Read_GFF ($_, \%gff);
}

my %stats;
foreach my $scaffold_id (keys %gff) {
	my $gff = $gff{$scaffold_id};
	ANNOTATION::stats_gene ($gff, \%stats);
}

$stats{'intron_total_number'} = $stats{'exon_total_number'} - $stats{'mRNA_total_number'};
$stats{'intron_total_length'} = $stats{'mRNA_total_length'} - $stats{'exon_total_length'};


printf "%s\t%d\n",   'Number of genes',         $stats{'gene_total_number'};
printf "%s\t%d\n",   'Total genic length',      $stats{'gene_total_length'};
printf "%s\t%d\n",   'Mean gene length',        $stats{'gene_total_length'} / $stats{'gene_total_number'};

printf "%s\t%d\n",   'Number of transcripts',   $stats{'mRNA_total_number'};
printf "%s\t%.1f\n", 'Transcripts per gene',    $stats{'mRNA_total_number'} / $stats{'gene_total_number'};
printf "%s\t%d\n",   'Total transcript length', $stats{'exon_total_length'};
printf "%s\t%d\n",   'Mean transcript length',  $stats{'exon_total_length'} / $stats{'mRNA_total_number'};

printf "%s\t%d\n",   'Number of exons',         $stats{'exon_total_number'};
printf "%s\t%.1f\n", 'Exons per transcript',    $stats{'exon_total_number'} / $stats{'mRNA_total_number'};
printf "%s\t%d\n",   'Mean exon length',        $stats{'exon_total_length'} / $stats{'exon_total_number'};
printf "%s\t%d\n",   'Number of coding exons',  $stats{'CDS_total_number'};

printf "%s\t%d\n",   'Number of introns',       $stats{'intron_total_number'};
printf "%s\t%d\n",   'Mean intron length',      $stats{'intron_total_length'} / $stats{'intron_total_number'};

printf "%s\t%d\n",   'Total cds length',        $stats{'CDS_total_length'};
printf "%s\t%d\n",   'Mean CDS length',         $stats{'CDS_total_length'}  / $stats{'mRNA_total_number'};


__END__
