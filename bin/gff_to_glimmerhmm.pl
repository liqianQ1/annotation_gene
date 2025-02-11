#!/home/biohuaxing/Personal/liqian/Anaconda3/envs/ann/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;

#use FindBin qw($Bin);
#use lib "$Bin/../lib/";
#use FASTA;
#use ANNOTATION;

my $usage=<<USAGE;
	Usage:	$0 <scaffold.fa>  <pasa.gff> [-o <output_prefix>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput);
GetOptions ("o:s"=>\$oput);
die $usage unless (@ARGV>=2);
my ($iput, $pasa) = @ARGV;
$oput ||= 'cds';


my $oput_cds = "$oput.fa";
my $oput_lst = "$oput.lst";


my %gff;
Read_GFF ($pasa, \%gff);


my $i=0;
open my $OCDS, ">$oput_cds" or die;
open my $OLST, ">$oput_lst" or die;
open my $IPUT, "<$iput" or die;
while (!eof $IPUT) {
	my $scaffold = read_fasta_one_seq ($IPUT);
	my $gff = $gff{$$scaffold{id}};
	$gff or next;
	extract_cds ($gff);
	foreach (@$gff) {
		my $start = $$_[0]{start}-1;
		my $end   = $$_[0]{end};
		my $len   = $end - $start;
		my %cds;
		$cds{id} = 'seq'.++$i;
		$cds{seq} = substr ($$scaffold{seq}, $start, $len);
		$cds{len} = $len;
		next if($cds{seq} =~ /NNNN/i); #fix by zhanghk
		write_fasta_one_seq ($OCDS, \%cds);
		for (my $j=1;$j<=$#$_;$j++) {
			if ($$_[0]{strand} eq '+') {
				print $OLST join(" ", $cds{id},$$_[$j]{start}-$start,$$_[$j]{end}-$start), "\n";
			}
			else {
				my $k = $#$_ - $j + 1;
				print $OLST join(" ", $cds{id},$$_[$k]{end}-$start,$$_[$k]{start}-$start), "\n";
			}
		}
		print $OLST "\n";
	}
}
close $IPUT;
close $OCDS;
close $OLST;

sub Read_GFF {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^\w/ or next;
		chomp;
		my @c = split /\t/;
		my $scaffold_id = $c[0];
		my %h;
		$h{program} = $c[1];
		$h{type}    = $c[2];
		$h{start}   = $c[3];
		$h{end}     = $c[4];
		$h{score}   = $c[5];
		$h{strand}  = $c[6];
		$h{phase}   = $c[7];
		$h{ann}     = $c[8]||'';
		($h{start}, $h{end}) = ($h{end}, $h{start}) if ($h{start}>$h{end});
		if ($h{ann} =~ /^ID=([^;]+);?(.*)$/) {
			$h{id}  = $1;
			$h{ann} = $2||'';
		}
		if ($h{ann} =~ /^Parent=([^;]+);?(.*)$/) {
			$h{parent} = $1;
			$h{ann}    = $2||'';
		}

		if (!$h{parent}) {
			my $i = $#{$$gff{$scaffold_id}} + 1;
			$$gff{$scaffold_id}[$i][0] = \%h;
		}
		elsif ($h{parent} eq $$gff{$scaffold_id}[-1][0]{id} && $#{$$gff{$scaffold_id}[-1]}>1 && $h{type} eq 'mRNA') {
			my $i = $#{$$gff{$scaffold_id}} + 1;
			%{$$gff{$scaffold_id}[$i][0]} = %{$$gff{$scaffold_id}[$i-1][0]};
			$$gff{$scaffold_id}[$i][0]{type} = '';
			$$gff{$scaffold_id}[$i][1] = \%h;
		}
		else {
			push @{$$gff{$scaffold_id}[-1]}, \%h;
		}
	}
	close IPUT;
}

sub extract_cds {
	my $gff = shift;
	foreach my $g (@$gff) {
		my @g;
		$g[0] = $$g[0];
		foreach ( sort {$$a{start}<=>$$b{start} || $$a{end}<=>$$b{end}} @$g) {
			$$_{type} eq 'CDS' and push @g, $_;
		}
		@$g = @g;
	}
}

sub read_fasta_one_seq {
	my $IPUT = shift;
	my $fasta;
	$/="\n";
	<$IPUT> =~ /\>?(.*?)\s+(.*)?$/ or die;
	$$fasta{id}  = $1;
	$$fasta{ann} = $2||'';
	$/='>';
	$$fasta{seq} = <$IPUT>;
	chomp $$fasta{seq};
	$$fasta{seq} =~ s/\n//g;
	$$fasta{seq} =~ s/\*$//;
	$$fasta{len} = length $$fasta{seq};
	$/="\n";
	return $fasta;
}

sub write_fasta_one_seq {
	my ($OPUT, $fasta) = @_;
	my $cn = 50;
	if ($$fasta{ann}) {
		print $OPUT ">$$fasta{id} $$fasta{ann}\n";
	}
	else {
		print $OPUT ">$$fasta{id}\n";
	}
	$$fasta{len} ||= length $$fasta{seq};
	for (my $i=0;$i<$$fasta{len};$i+=$cn) {
		print $OPUT (substr($$fasta{seq},$i,$cn)), "\n";
	}
}

__END__
