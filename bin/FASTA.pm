package FASTA;
use strict;
use warnings;
use POSIX;

use FindBin qw($Bin);
use lib "$Bin/";
use TOOLS;

############################################################

sub read_fasta_all_seq {
	my $iput = shift;
	my $fasta;
	open my $IPUT, "<$iput" or die;
	while (!eof $IPUT) {
		push @$fasta, read_fasta_one_seq ($IPUT);
	}
	close $IPUT;
	return $fasta;
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

############################################################

sub write_fasta_all_seq {
	my ($oput, $fasta) = @_;
	open my $OPUT, ">$oput" or die;
	if (ref $fasta eq 'ARRAY') {
		foreach (@$fasta) {
			write_fasta_one_seq ($OPUT, $_);
		}
	}
	else {
			write_fasta_one_seq ($OPUT, $fasta);
	}
	close $OPUT;
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

############################################################

sub contig_of_scaffold {
	my ($scaffold, $n) = @_;
	my $contig = ();
	my $gap    = ();
	$n ||= 1;
	$$scaffold{gap_num}=0;
	$$scaffold{gap_len}=0;
	my $i=0;
	my $j=0;
	my $start=0;
	my $end=0;
	while ($$scaffold{seq} =~ /([Nn\.]{$n,})?(((?![Nn\.]{$n,}).)*)?/g) {
		if ($1) {
			my $length = length $1;
			$$scaffold{gap_num}++;
			$$scaffold{gap_len}+= $length;
			$start = $end + 1;
			$end   = $start + $length - 1;
			$$gap[$j]{start} = $start;
			$$gap[$j]{end}   = $end;
			$$gap[$j]{len}   = $length;
			$j++;
		}
		if ($2) {
			my $length = length $2;
			$$contig[$i]{seq} = $2;
			$start = $end + 1;
			$end   = $start + $length - 1;
			$$contig[$i]{start} = $start;
			$$contig[$i]{end}   = $end;
			$$contig[$i]{len}   = $length;
			$i++;
		}
	}
	return ($contig, $gap);
}

############################################################

sub sort_by_strand {
	my ($fasta, $strand_tab, $sort) = @_;
	my ($key, $tab) = TOOLS::read_table($strand_tab);
	my %id;
	for (my $i=0;$i<@$fasta;$i++) {
		exists $id{$$fasta[$i]{id}} and die;
		$id{$$fasta[$i]{id}} = $i;
	}
	foreach (@$tab) {
		exists $id{$$_{id}} or die;
		if ($$_{strand} eq '-') {
			$$fasta[$id{$$_{id}}]{seq} = reverse $$fasta[$id{$$_{id}}]{seq};
			$$fasta[$id{$$_{id}}]{seq} =~ tr/ACGTacgt/TGCAtgca/;
		}
		push @$sort, $$fasta[$id{$$_{id}}];
		delete $id{$$_{id}};
	}
	foreach (sort {$$b{len}<=>$$a{len}} @$fasta) {
		exists $id{$$_{id}} or next;
		push @$sort, $_;
	}

}

############################################################

sub length_stats {
	my $lens = shift;
	my %stats;
	@{$lens} = sort {$b<=>$a} @{$lens};
	$stats{Num} = @{$lens};
	$stats{Len} += $_ foreach (@{$lens});
	$stats{MAX} = $$lens[0];
	$stats{MIN} = $$lens[-1];

	my %length;
	for (my $i=10; $i<100; $i+=10) {
		$length{"N$i"} = $stats{Len} * $i /100;
	}

	my $length = 0;
	my $i = 10;
	foreach (@{$lens}) {
		$i > 90 and last;
		$length += $_;
		if ($length >= $length{"N$i"}) {
			$stats{"N$i"} = $_;
			$i+=10;
		}
	}

	return \%stats;
}


sub scaffold_stats {
	my ($stats, $oput) = @_;
	open OPUT, ">$oput" or die;
	my @key = qw (Contig_Num  Contig_Len  GC  Scaffold_Num  Scaffold_Len  Gap_Num  Gap_Len);
	$$stats{Contig_Num}   = $$stats{contig}{Num};
	$$stats{Contig_Len}   = $$stats{contig}{Len};
	$$stats{GC}           = sprintf "%.2f", $$stats{GC}/$$stats{contig}{Len}*100;
	$$stats{Scaffold_Num} = $$stats{scaffold}{Num};
	$$stats{Scaffold_Len} = $$stats{scaffold}{Len};
	print OPUT (TOOLS::print_title \@key);
	print OPUT (TOOLS::print_value \@key, $stats), "\n";

	@key = qw(- MAX N10 N20 N30 N40 N50 N60 N70 N80 N90 MIN);
	$$stats{contig}{'-'}   = 'Contig';
	$$stats{scaffold}{'-'} = 'Scaffold';
	print OPUT (TOOLS::print_title \@key);
	print OPUT (TOOLS::print_value \@key, $$stats{contig});
	print OPUT (TOOLS::print_value \@key, $$stats{scaffold});
	close OPUT;
}

############################################################

sub join_fasta_by_nnn {
	my $fasta1 = shift;
	my $fasta2;
	$$fasta2{id}  = 'joined';
	$$fasta2{seq} = '';
	foreach (@$fasta1) {
		$$fasta2{seq} .= $$_{seq}.'N'x2000;
	}
	$$fasta2{len} = length $$fasta2{seq};
	return $fasta2;
}

############################################################

sub filter_fasta_by_len {
	my ($fasta1, $len) = @_;
	my $fasta2;
	foreach (@$fasta1) {
		push @$fasta2, $_ if ($$_{len} >= $len);
	}
	return $fasta2;
}

############################################################

sub read_fasta_id {
	my $iput = shift;
	my $id;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/\>(\S+)/ or next;
		$$id{$1} and die;
		$$id{$1} = 1;
	}
	close IPUT;
	return $id;
}

############################################################

sub split_fasta_by_number {
	my ($iput, $S, $oput) = @_;
	my @oput;
	my @ls = `ls $oput.*`;
	if (@ls) {
		foreach (@ls) {
			chomp;
			push @oput, $_;
		}
		return @oput;
	}
	$S ||= 1;
	if ($S==1) {
		$oput[0] = "$oput.00";
		`ln -s $iput $oput[0]`;
	}
	else {
		my $fasta1 = read_fasta_all_seq ($iput);
		my $T = @$fasta1;
		my $N = POSIX::ceil $T/$S;
		my $i = 0;
		my $j = 0;
		my $k = 0;
		my @fasta2;
		foreach (@$fasta1) {
			$i++;
			$j++;
			push @fasta2, $_;
			if ($j>=$N or $i==@$fasta1) {
				$oput[$k] = sprintf "$oput.%02d", $k;
				write_fasta_all_seq ($oput[$k], \@fasta2);
				@fasta2 = ();
				$k++;
				$j = 0;
			}
		}
	}
	return @oput;
}


sub split_fasta_by_length {
	my ($iput, $S, $oput) = @_;
	my @oput;
	my @ls = `ls $oput.*`;
	if (@ls) {
		foreach (@ls) {
			chomp;
			push @oput, $_;
		}
		return @oput;
	}
	$S ||= 1;
	if ($S==1) {
		$oput[0] = "$oput.00";
		`ln -s $iput $oput[0]`;
	}
	else {
		my $fasta1 = read_fasta_all_seq ($iput);
		my $T = 0;
		$T += $$_{len} foreach (@$fasta1);
		my $N = POSIX::ceil $T/$S;
		my $i = 0;
		my $j = 0;
		my $k = 0;
		my @fasta2;
		foreach (@$fasta1) {
			$i++;
			$j+= $$_{len};
			push @fasta2, $_;
			if ($j>=$N or $i==@$fasta1) {
				$oput[$k] = sprintf "$oput.%02d", $k;
				write_fasta_all_seq ($oput[$k], \@fasta2);
				@fasta2 = ();
				$k++;
				$j = 0;
			}
		}
	}
	return @oput;
}

############################################################

1;

__END__
