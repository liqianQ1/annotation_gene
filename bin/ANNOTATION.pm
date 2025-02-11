package ANNOTATION;
use strict;
use warnings;
use 5.010;
use Data::Dumper;

use FindBin qw($Bin);
use lib "$Bin/";
use FASTA;
use TOOLS;

############################################################

sub sort_scaffold_id {
	my $scaffold = shift;
	my %scaffold_id;
	my $flag = 100000;
	foreach (keys %$scaffold) {
		#/(\d+)/ or die;
		if(/(\d+)/){
			$scaffold_id{$_} = $1;
		}else{
			$scaffold_id{$_} = $flag;
			$flag++;
		}
	}
	my @scaffold_id = sort {$scaffold_id{$a}<=>$scaffold_id{$b}} keys %scaffold_id;
	return @scaffold_id;
}


sub sort_gff {
	my $gff = shift;
	@$gff = sort {$$a[0]{start}<=>$$b[0]{start} || $$b[0]{end}<=>$$a[0]{end} || $$a[0]{score}=~/\d/ && $$b[0]{score}=~/\d/ && $$b[0]{score}<=>$$a[0]{score}} @$gff;
}


sub sort_gene {
	my $gff = shift;
	foreach my $g (@$gff) {
		@$g = sort {$$a{start}<=>$$b{start} || $$b{end}<=>$$a{end}} @$g;
	}
}

############################################################

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


sub Print_GFF {
	my ($OPUT, $gff, $scaffold_id, $scaffold_len) = @_;
	print $OPUT "##sequence-region $scaffold_id 1 $scaffold_len\n" if ($scaffold_len);
	foreach my $g (@$gff) {
		$$g[0]{filter} and next;
		foreach (@$g) {
			$$_{type} or next;
			my $txt = '';
			$txt .= "$scaffold_id\t$$_{program}\t$$_{type}\t$$_{start}\t$$_{end}\t$$_{score}\t$$_{strand}\t$$_{phase}\t";
			$txt .= "ID=$$_{id};" if (defined $$_{id});
			$txt .= "Parent=$$_{parent};" if (defined $$_{parent});
			$txt .= "$$_{ann}" if ($$_{ann});
			print $OPUT "$txt\n";
		}
	}
}


sub Print_FFN {
	my ($OPUT, $gff) = @_;
	foreach my $g (@$gff) {
		$$g[0]{filter} and next;
		$$g[0]{fasta}{id} ||= $$g[0]{id};
		FASTA::write_fasta_one_seq ($OPUT, $$g[0]{fasta});
	}
}

############################################################

sub Read_RepeatMasker {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^(\s+)?\d/ or next;
		s/^\s+//;
		chomp;
		my @c = split /\s+/;
		my $scaffold_id = $c[4];
		my %h;
		$h{program} = 'RepeatMasker';
		$h{type}    = 'similarity';
		$h{start}   = $c[5];
		$h{end}     = $c[6];
		$h{score}   = $c[0];
		$h{strand}  = ($c[8] eq '+') ? '+' : '-';
		$h{phase}   = '.';
		my $target  = $c[9];
		my $class   = $c[10];
		my ($target_start, $target_end) = ($c[11]=~/[\(\)]/)? ($c[13],$c[12]) : ($c[11],$c[12]);
		$h{ann}     = "Target=$target $target_start $target_end;Class=$class;PercDiv=$c[1];PercDel=$c[2];PercIns=$c[3];";
		($h{start}, $h{end}) = ($h{end}, $h{start}) if ($h{start}>$h{end});
		my $i = $#{$$gff{$scaffold_id}} + 1;
		$$gff{$scaffold_id}[$i][0] = \%h;
	}
	close IPUT;
}


sub Read_ProteinMask {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^\d/ or next;
		chomp;
		my @c = split /\s+/;
		my $scaffold_id = $c[3];
		my %h;
		$h{program} = 'RepeatProteinMask';
		$h{type}    = 'TEprotein';
		$h{start}   = $c[4];
		$h{end}     = $c[5];
		$h{score}   = $c[1];
		$h{strand}  = $c[6];
		$h{phase}   = '.';
		my $target  = $c[7];
		my $class   = $c[8];
		my ($target_start, $target_end) = ($c[9]<$c[10])? ($c[9],$c[10]):($c[10],$c[9]);
		$h{ann}     = "Target=$target $target_start $target_end;Class=$class;";
		($h{start}, $h{end}) = ($h{end}, $h{start}) if ($h{start}>$h{end});
		my $i = $#{$$gff{$scaffold_id}} + 1;
		$$gff{$scaffold_id}[$i][0] = \%h;
	}
	close IPUT;
}


sub Read_TRF {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	my $scaffold_id;
	while (<IPUT>) {
		$scaffold_id = $1 if(/^Sequence:\s+(\S+)/);
		/^\d+/ or next;
		chomp;
		my @c = split /\s+/;
		my %h;
		$h{program} = 'TRF';
		$h{type}    = 'TandemRepeat';
		$h{start}   = $c[0];
		$h{end}     = $c[1];
		$h{score}   = $c[7];
		$h{strand}  = '+';
		$h{phase}   = '.';
		my $class   = 'TandemRepeat';
		$h{ann}     = "Class=$class;PeriodSize=$c[2];CopyNumber=$c[3];PercentMatches=$c[5];PercentIndels=$c[6];Consensus=$c[13];";
		($h{start}, $h{end}) = ($h{end}, $h{start}) if ($h{start}>$h{end});
		my $i = $#{$$gff{$scaffold_id}} + 1;
		$$gff{$scaffold_id}[$i][0] = \%h;
	}
	close IPUT;
}

############################################################


sub Read_Exonerate {
	my ($iput, $gff, $program) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^\S+\texonerate/ or next;
		chomp;
		my @c = split /\t/;
		my $scaffold_id = $c[0];
		my %h;
		$h{program} = $program;
		$h{type}    = $c[2];
		$h{start}   = $c[3];
		$h{end}     = $c[4];
		$h{score}   = $c[5];
		$h{strand}  = $c[6];
		$h{phase}   = $c[7];
		$h{ann}     = '';
		($h{start}, $h{end}) = ($h{end}, $h{start}) if ($h{start}>$h{end});

		if ($h{type} eq 'gene') {
			my $i = $#{$$gff{$scaffold_id}} + 1;
			$$gff{$scaffold_id}[$i][0] = \%h;
			%{$$gff{$scaffold_id}[$i][1]} = %h;
			$$gff{$scaffold_id}[$i][1]{type} = 'mRNA';
		}
		elsif ($h{type} =~ /cds/) {
			$h{type} = 'CDS';
			my $i = $#{$$gff{$scaffold_id}[-1]} + 2;
			$$gff{$scaffold_id}[-1][$i] = \%h;
		}
		elsif ($h{type} =~ /exon/) {
			my $i = $#{$$gff{$scaffold_id}[-1]} - 1;
			$$gff{$scaffold_id}[-1][$i] = \%h;
		}
	}
	close IPUT;
}


sub Read_GTF {
	my ($iput, $gff) = @_;
	my $transcript_id = '';
	my $gene_id = '';
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
		$h{score}   = $c[5]||'.';
		$h{strand}  = $c[6];
		$h{phase}   = $c[7];
		($h{start}, $h{end}) = ($h{end}, $h{start}) if ($h{start}>$h{end});

		$c[8] =~ /transcript_id\s"(\S+)"/ or die;
		if ($1 ne $transcript_id) {
			$transcript_id = $1;
			my $i = $#{$$gff{$scaffold_id}} + 1;
			my %h0 = %h;
			my %h1 = %h;
			$c[8] =~ /gene_id\s"(\S+)"/ or die;
			$h0{type} = ($gene_id eq $1)? '':'gene';
			$h0{id}   = $1;
			$gene_id = $1;
			$h1{type} = 'mRNA';
			$h1{id}   = $transcript_id;
			$h1{parent} = $gene_id;
			if ($h{type} eq 'transcript') {
				$$gff{$scaffold_id}[$i][0] = \%h0;
				$$gff{$scaffold_id}[$i][1] = \%h1;
			}
			else {
				$h0{end}  = 0;
				$h1{end}  = 0;
				$h{parent} = $transcript_id;
				$$gff{$scaffold_id}[$i][0] = \%h0;
				$$gff{$scaffold_id}[$i][1] = \%h1;
				$$gff{$scaffold_id}[$i][2] = \%h;
			}
		}
		else {
			$h{parent} = $transcript_id;
			push @{$$gff{$scaffold_id}[-1]}, \%h;
		}
	}
	foreach my $scaffold_id (keys %$gff) {
		foreach (@{$$gff{$scaffold_id}}) {
			$$_[0]{end} ||= $$_[-1]{end};
			$$_[1]{end} ||= $$_[-1]{end};
		}
	}
	close IPUT;
}


sub Read_Glimmer {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	my $scaffold_id;
	while (<IPUT>) {
		/^>(\S+)/ and $scaffold_id = $1 and next;
		my ($id, $start, $end, $strand, $score) = (split /\s+/, $_)[0,1,2,3,4];
		$strand =~ s/\d//;
		if ($strand eq '-') {
			($start, $end) = ($end, $start);
		}
		$end > $start or next;
		my @g = ();
		$g[0]{program} = 'GLIMMER';
		$g[0]{type}    = 'gene';
		$g[0]{start}   = $start;
		$g[0]{end}     = $end;
		$g[0]{score}   = '.';
		$g[0]{strand}  = $strand;
		$g[0]{phase}   = '.';
		$g[0]{id}      = $id;
		$g[1]{program} = 'GLIMMER';
		$g[1]{type}    = 'CDS';
		$g[1]{start}   = $start;
		$g[1]{end}     = $end;
		$g[1]{score}   = $score;
		$g[1]{strand}  = $strand;
		$g[1]{phase}   = 0;
		$g[1]{parent}  = $id;
		push @{$$gff{$scaffold_id}}, \@g;
	}
	close IPUT;
}


sub Read_tRNA {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (! eof IPUT) {
		<IPUT> =~ /^(.*).trna(\d+).*?(\d+)-(\d+).*?Length: (\d+)/ or next;
		my $scaffold_id = $1;
		my %h;
		($h{id}, $h{start}, $h{end}) = ("tRNA$2", $3, $4);
		<IPUT> =~ /^Type: (\S+)\s+Anticodon: (\S\S\S) .*?Score: (\S+)/ or die;
		($h{tRNA_tp}, $h{tRNA_ac}, $h{score}) = ($1, $2, $3);
		<IPUT> =~ /^Possible/ and next;
		<IPUT> =~ /Seq: (\S+)$/ or die;
		my $seq = $1;
		<IPUT> =~ /Str: (\S+)$/ or die;
		$h{tRNA_str} = $1;

		if ($h{end} > $h{start}) {
			$h{strand} = '+';
		}
		else {
			$h{strand} = '-';
			($h{start}, $h{end}) = ($h{end}, $h{start});
		}
		$h{program} = 'tRNAscan-SE';
		$h{type}    = 'tRNA';
		$h{phase}   = '.';
		$h{ann} = "product=tRNA-$h{tRNA_tp};anticodon=$h{tRNA_ac};";
		$h{fasta}{seq} = uc $seq;
		$h{fasta}{ann} = "locus=$scaffold_id:$h{start}:$h{end}:$h{strand}; $h{ann}";
		my $i = $#{$$gff{$scaffold_id}} + 1;
		$$gff{$scaffold_id}[$i][0] = \%h;
	}
	close IPUT;
}


sub Read_Rfam {
	my ($iput, $gff) = @_;
	state $rfam = read_rfam_tab ($iput);
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^#/ and next;
		chomp;
		my $scaffold_id;
		my %h;
		my $ann;
		($scaffold_id, $h{start}, $h{end}, $h{score}, $h{strand}, $ann) = (split /\t/, $_)[0,3,4,5,6,8];

		$ann =~ /;id=([^;]+);.*?rfam-acc=([^;]+);rfam-id=([^;]+)/ or die;
		($h{id}, $h{rfam_acc}, $h{rfam_id}) = ($1, $2, $3);
		$h{rfam_id} =~ /$$rfam{$h{rfam_acc}}{ID}(\.\d+)?/ or die;

		$h{program} = 'Rfam';
		$h{type}    = 'ncRNA';
		$h{phase}   = '.';
		$h{ann}     = "ncrna_class=$$rfam{$h{rfam_acc}}{Class};rfam_id=$h{rfam_id};rfam_acc=$h{rfam_acc};rfam_tp=$$rfam{$h{rfam_acc}}{TP};rfam_de=\"$$rfam{$h{rfam_acc}}{DE}\";";
		my $i = $#{$$gff{$scaffold_id}} + 1;
		$$gff{$scaffold_id}[$i][0] = \%h;
	}
	close IPUT;
}

#### 
#AUTHOR = Wu Di 
####



sub Read_rRNA2{
	my ($iput, $gff) = @_;
	open FH,"cat $iput | tr -s \" \" | sed  \"s\/ \/\t\/g\" | " 
	|| die "Can't open file";
	my $cmfetch = "/public/software/Infernal/bin/cmfetch";
	my $Rfam = "/public/database/Rfam/CMs/Rfam.cm";
	while (<FH>) {
		next unless !/^#/;
		chomp;
		my %h;
		my $coverage;
		my ($target_name,$Taccession,$query_name,
		$Qaccession,$mdl,$mdl_from,$mdl_to,
		$seq_from,$seq_to,$strand,$trunc,
		$pass,$gc,$bias,$score,$Evalue,$inc,$description_of_target)  = split "\t";
		my $target_len = `$cmfetch $Rfam  $Taccession | grep CLEN | tr -s ' ' | cut -d ' ' -f2`;
		($h{strand}, $h{start}, $h{end}, $h{score}) = ($strand , $seq_from , $seq_to ,$score);

		$coverage = sprintf "%.2f", 100*($mdl_to-$mdl_from+1)/$target_len; ###
		$h{program} = 'INFERNAL';
		$h{type}    = 'rRNA';
		$h{phase}   = '.';
		$h{ann}     = "Target=$target_name $mdl_from $mdl_to;coverage=$coverage;";
		my $i = $#{$$gff{$query_name}} + 1;
		$$gff{$query_name}[$i][0] = \%h;
	}
	close FH;
}

sub Read_rRNA {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^#/ and next;
		chomp;
		my $scaffold_id;
		my %h;
		my ($target, $target_len, $target_start, $target_end, $coverage);
		($target, $target_len, $target_start, $target_end, $h{strand}, $scaffold_id, $h{start}, $h{end}, $h{score}) = (split /\t/, $_)[0,1,2,3,4,5,7,8,10];
		$target =~ /(rRNA_[\d\.]+S)/ or die;
		$target = $1;
		$coverage = sprintf "%.2f", 100*($target_end-$target_start+1)/$target_len;
		$h{program} = 'BLASTN';
		$h{type}    = 'rRNA';
		$h{phase}   = '.';
		$h{ann}     = "Target=$target $target_start $target_end;coverage=$coverage;";
		my $i = $#{$$gff{$scaffold_id}} + 1;
		$$gff{$scaffold_id}[$i][0] = \%h;
	}
	close IPUT;
}


sub Read_Solar {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^#/ and next;
		chomp;
		my ($target, $target_len, $target_start, $target_end, $strand, $scaffold_id, $start, $end, $score) = (split /\t/, $_)[0,1,2,3,4,5,7,8,10];
		my $coverage = sprintf "%.2f", 100*($target_end-$target_start+1)/$target_len;
		$target_end-$target_start+1 == $target_len or next;
		my @g = ();
		$g[0]{program} = 'BLASTN';
		$g[0]{type}    = 'gene';
		$g[0]{start}   = $start;
		$g[0]{end}     = $end;
		$g[0]{score}   = '.';
		$g[0]{strand}  = $strand;
		$g[0]{phase}   = '.';
		$g[0]{ann}     = "Target=$target $target_start $target_end;coverage=$coverage;";
		$g[1]{program} = 'BLASTN';
		$g[1]{type}    = 'CDS';
		$g[1]{start}   = $start;
		$g[1]{end}     = $end;
		$g[1]{score}   = $score;
		$g[1]{strand}  = $strand;
		$g[1]{phase}   = 0;
		$g[1]{ann}     = '';
		push @{$$gff{$scaffold_id}}, \@g;
	}
	close IPUT;
}


sub Read_Blast {
	my ($iput, $blast) = @_;
	open my $IPUT, "<$iput" or die;
	my @key;
	while (<$IPUT>) {
		/^# Fields: (.*)$/ or next;
		@key = split ', ', $1;
		last;
	}
	@key or die;
	my %title = ('query id'          => 'Query_ID',
	             'subject id'        => 'Subject_ID',
	             'query length'      => 'Query_len',
	             'subject length'    => 'Subject_len',
	             'evalue'            => 'Evalue',
	             '% identity'        => '%Ident',
	             '% positives'       => '%Pos',
	             'alignment length'  => 'Align_len',
	             'mismatches'        => 'Mismatch',
	             'gap opens'         => 'Gap_opening',
	             'q. start'          => 'Query_start',
	             'q. end'            => 'Query_end',
	             's. start'          => 'Subject_start',
	             's. end'            => 'Subject_end',
	             'bit score'         => 'Score',
	             'subject tax ids'   => 'Tax_ID',
	             'subject sci names' => 'Sci_Name',
	             'subject title'     => 'Subject_Title');
	for (my $i=0;$i<@key;$i++) {
		$key[$i] = $title{$key[$i]} || $key[$i];
	}
	while (<$IPUT>) {
		/^#/ and next;
		my $h = TOOLS::read_value(\@key, $_);
		push @$blast, $h;
	}
	close $IPUT;
}

############################################################

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


sub extract_fasta {
	my ($scaffold, $gff) = @_;
	foreach my $g (@$gff) {
		my %s;
		my ($len, $seq);
		if (@$g > 1) {
			for (my $i=1;$i<@$g;$i++) {
				$$g[$i]{type} eq 'CDS' or next;
				$len  = $$g[$i]{end} - $$g[$i]{start} + 1;
				$seq .= substr $$scaffold{seq}, $$g[$i]{start}-1, $len;
			}
		}
		elsif (@$g == 1) {
			$len = $$g[0]{end} - $$g[0]{start} + 1;
			$seq = substr $$scaffold{seq}, $$g[0]{start}-1, $len;
		}
		if ($$g[0]{strand} eq '-') {
			$seq = reverse $seq;
			$seq =~ tr/ACGTRYMK/TGCAYRKM/;
			$seq =~ tr/acgtrymk/tgcayrkm/;
		}
		$s{seq} = uc $seq;
		my $N = $seq =~ tr/Nn//;
		$s{gap} = $N;
		$s{ann} = "locus=$$scaffold{id}:$$g[0]{start}:$$g[0]{end}:$$g[0]{strand};";
		$s{ann}.= " $$g[0]{ann}" if ($$g[0]{ann});
		$$g[0]{fasta} = \%s;
	}
}


sub extract_fa {
	my ($gff, $scaffold) = @_;
	foreach my $g (@$gff) {
		my %s;
		my ($len, $seq);
		$len = $$g[0]{end} - $$g[0]{start} + 1;
		$seq = substr $$scaffold{seq}, $$g[0]{start}-1, $len;
		if ($$g[0]{strand} eq '-') {
			$seq = reverse $seq;
			$seq =~ tr/ACGTRYMK/TGCAYRKM/;
			$seq =~ tr/acgtrymk/tgcayrkm/;
		}
		$s{seq} = uc $seq;
		my $N = $seq =~ tr/Nn//;
		$s{gap} = $N;
		$s{ann} = "locus=$$scaffold{id}:$$g[0]{start}:$$g[0]{end}:$$g[0]{strand};";
		$s{ann}.= " $$g[0]{ann}" if ($$g[0]{ann});
		$$g[0]{fasta} = \%s;
	}
}


sub unique_gene_id {
	my ($gff, $count, $scaffold_id, $program) = @_;
	$program ||= $$gff[0][0]{program};
	my $number;
	my $parent;
	my $m = 1;
	foreach my $g (@$gff) {
		$$g[0]{filter} and next;
		if ($$g[0]{type} eq 'gene') {
			$number = sprintf "%05d", ++$$count;
			$$g[0]{program} = $program;
			$$g[0]{id} = "$scaffold_id.g$number";
			$parent = $$g[0]{id};
			$m = 1;
		}
		elsif (!$$g[0]{type}) {
			$$g[0]{program} = $program;
			$$g[0]{id} = "$scaffold_id.g$number";
			$parent = $$g[0]{id};
		}
		else {
			die;
		}
		my %n;
		$n{exon} = 0;
		$n{utr5p} = 0;
		$n{utr3p} = 0;
		for (my $i=1;$i<@$g;$i++) {
			my $type = $$g[$i]{type};
			if ($type eq 'mRNA') {
				$$g[$i]{program} = $program;
				$$g[$i]{id} = "$parent.m$m";
				$$g[$i]{parent}  = $$g[0]{id};
				$parent = $$g[$i]{id};
				$m++;
				$n{$_} = 0 foreach (keys %n);
			}
			else {
				if ($type eq 'CDS') {
					$type = 'cds';
				}
				elsif ($type eq 'five_prime_UTR') {
					$type = 'utr5p';
				}
				elsif ($type eq 'three_prime_UTR') {
					$type = 'utr3p';
				}
				$$g[$i]{program} = $program;
				$$g[$i]{parent} = $parent;
				$$g[$i]{id} = "$parent.$type";
				if ($type =~ /^exon$|^utr5p$|^utr3p$/) {
					$n{$type}++;
					$$g[$i]{id} .= $n{$type};
				}
			}
		}
	}
}


sub unique_id {
	my ($gff, $count) = @_;
	foreach my $g (@$gff) {
		$$g[0]{filter} and next;
		my $type = $$g[0]{type};
		my $number = ++$$count{$type};
		if ($type eq 'mRNA') {
			$number  = sprintf "%05d", $number;
			$$g[0]{id} = "$type$number";
		}
		elsif ($type =~ /^gene|^rRNA/) {
			$number  = sprintf "%04d", $number;
			$$g[0]{id} = "$type$number";
		}
		elsif ($$g[0]{program} =~ /^tRNA|^Rfam/) {
			$number  = sprintf "%04d", $number;
			$$g[0]{id} = "$type$number";
		}
		if (@$g > 1) {
			for (my $i=1;$i<@$g;$i++) {
				$type = $$g[$i]{type};
				++$$count{$type};
				$$g[$i]{parent} = $$g[0]{id};
			}
		}
		$$g[0]{fasta}{id} = $$g[0]{id} if ($$g[0]{fasta});
	}
}


sub stats_gene {
	my ($gff, $stats) = @_;
	foreach my $g (@$gff) {
		foreach (@$g) {
			my $type = $$_{type};
			$type or next;
			$$stats{"${type}_total_number"}++;
			$$stats{"${type}_total_length"}+= $$_{end} - $$_{start} + 1;
		}
	}
}


############################################################


sub filter_gap {
	my ($gff, $N) = @_;
	$N ||= 0;
	foreach my $g (@$gff) {
		$$g[0]{filter} = 1 if ($$g[0]{fasta}{gap} > $N);
	}
}


sub filter_rfam {
	my ($OPUT, $gff, $scaffold_id, $filter_rfam_id) = @_;
	foreach (@$gff) {
		if ($$_[0]{program} eq 'Rfam') {
			$$_[0]{ann} =~ /rfam_id=(.*?);/ or die;
			my $rfam_id = $1;
			if ($$filter_rfam_id{$rfam_id}) {
				my @gff = ($_);
				$$_[0]{ann} = '#'.$$_[0]{ann};
				Print_GFF ($OPUT, \@gff, $scaffold_id);
				$$_[0]{filter} = 1;
			}
		}
	}
}


sub filter_overlap {
	my ($OPUT, $gff, $scaffold_id, $program_level, $per) = @_;
	$per ||= 0.8;
	my ($old, $old_len);
	foreach (@$gff) {
		$$_[0]{filter} and next;
		if (!$old) {
			$old = $_;
			$old_len  = $$_[0]{end} - $$_[0]{start} + 1;
			next;
		}
		my $new = $_;
		my $new_len  = $$_[0]{end} - $$_[0]{start} + 1;

		my ($new_filter, $old_filter) = (0,0);
		my $over_lap = $$old[0]{end} - $$new[0]{start} + 1;
		if ($over_lap>0 && ($over_lap/$old_len>=$per || $over_lap/$new_len>=$per)) {
			my ($old_prg, $new_prg) = ($$program_level{$$old[0]{program}}||99, $$program_level{$$new[0]{program}}||99);
			if ($old_prg < $new_prg) {
				$old_filter = 1;
			}
			elsif ($old_prg > $new_prg) {
				$new_filter = 1;
			}
			elsif ($$old[0]{score} <= $$new[0]{score}) {
				$old_filter = 1;
			}
			else {
				$new_filter = 1;
			}
		}

		if ($old_filter) {
			$$old[0]{ann} = '*'.$$old[0]{ann};
			my @gff = ($old, $new);
			Print_GFF ($OPUT, \@gff, $scaffold_id);
			$$old[0]{filter} = 1;
			$old = $new;
			$old_len = $$old[0]{end} - $$old[0]{start} + 1;
		}
		elsif ($new_filter) {
			$$new[0]{ann} = '*'.$$new[0]{ann};
			my @gff = ($old, $new);
			Print_GFF ($OPUT, \@gff, $scaffold_id);
			$$new[0]{filter} = 1;
		}
		else {
			$old = $new;
			$old_len = $new_len;
		}
	}
}


sub filter_miscodon {
	my ($gff) = @_;
	my (%stop_codon, %start_codon);
	$stop_codon{TAG} = 1;
	$stop_codon{TAA} = 1;
	$stop_codon{TGA} = 1;
	$start_codon{ATG} = 1;
	$start_codon{CTG} = 1;
	$start_codon{GTG} = 1;
	$start_codon{TTG} = 1;
	foreach my $g (@$gff) {
		my $stop  = substr ($$g[0]{fasta}{seq}, $$g[0]{fasta}{len}-3, 3);
		my $start = substr ($$g[0]{fasta}{seq}, 0, 3);
		$$g[0]{filter} = 1 if (!$stop_codon{$stop});
		$$g[0]{filter} = 1 if (!$start_codon{$start});
	}
}


############################################################

sub read_rfam_tab {
	my $iput_gff = shift;
	my $rfam_tab;
	open IPUT, "<$iput_gff" or die;
	while (<IPUT>) {
		/^# CM file:\s+(\S+)$/ or next;
		$rfam_tab = $1;
		$rfam_tab =~ s/(\w+)$/tab/;
		last;
	}
	close IPUT;
	my $rfam;
	my $tab = TOOLS::read_table ($rfam_tab);
	foreach (@{$tab}) {
		$$rfam{$$_{AC}} = $_;
	}
	return $rfam;
}

############################################################

sub Read_SIM4 {
	my ($iput, $gff) = @_;
	state %target;
	open IPUT, "<$iput" or die;
	my ($scaffold_id, $scaffold_len, $target, $target_len, $identity, $coverage, @exon);
	my $strand = '+';
	<IPUT>;
	while (<IPUT>) {
		if (/seq1 = (\S+), (\d+) bp/) {
			($target, $target_len) = ($1, $2);
		}
		elsif (/seq2 = (\S+), (\d+) bp/) {
			($scaffold_id, $scaffold_len) = ($1, $2);
			<IPUT>;
		}
		elsif (/\(complement\)/) {
			$strand = '-';
		}
		elsif (/(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d+)\%/) {
			if (@exon && $3-$exon[-1][3] <= 10) {
				($exon[-1][1], $exon[-1][3]) = ($2, $4);
				$exon[-1][4] = ($exon[-1][4]+$5)/2;
			}
			else {
				push @exon, [$1,$2,$3,$4,$5];
			}
			$identity += ($2 - $1 + 1) * $5 / 100;
			$coverage += ($2 - $1 + 1);
		}
		if (/^\n/ || eof IPUT) {
			$identity = $identity / $coverage;
			$coverage = $coverage / $target_len;
			if ($identity>=0.95 and $coverage>=0.95) {
				my @g = ();
				$g[0]{program} = 'blat';
				$g[0]{type}    = 'mRNA';
				$g[0]{start}   = $exon[0][2];
				$g[0]{end}     = $exon[-1][3];
				$g[0]{score}   = sprintf "%.4f", $identity * $coverage;
				$g[0]{strand}  = $strand;
				$g[0]{phase}   = '.';
				$g[0]{ann}     = "Target=$target $exon[0][0] $exon[-1][1];";
				$g[0]{id}      = $target;
				$target{$target}++;
				if ($target{$target}>1) {
					$g[0]{id} .= "_D$target{$target}";
				}
				for (my $i=0;$i<@exon;$i++) {
					my $j = $i+1;
					$g[$j]{program} = 'blat';
					$g[$j]{type}    = 'exon';
					$g[$j]{start}   = $exon[$i][2];
					$g[$j]{end}     = $exon[$i][3];
					$g[$j]{score}   = $exon[$i][4];
					$g[$j]{strand}  = $strand;
					$g[$j]{phase}   = '.';
					$g[$j]{parent}  = $g[0]{id};
					$g[$j]{ann}     = "Target=$target $exon[$i][0] $exon[$i][1];";
				}
				push @{$$gff{$scaffold_id}}, \@g;
			}
			($scaffold_id, $scaffold_len, $target, $target_len, $identity, $coverage, @exon) = ();
			$strand = '+';
		}
	}
	close IPUT;
}


sub cluster {
	my $gff = shift;
	my $gff_cluster = ();
	my ($old, $old_end);
	foreach (@$gff) {
		if (!$old) {
			$old = $_;
			$old_end = $$old[0]{end};
			$$gff_cluster[0][0] = $old;
			next;
		}
		my $new = $_;
		my $over_lap = $old_end - $$new[0]{start} + 1;
		if ($over_lap > 0) {
			push @{$$gff_cluster[-1]}, $new;
		}
		else {
			my $i = $#{$gff_cluster} + 1;
			$$gff_cluster[$i][0] = $new;
		}
		$old = $new;
		$old_end = $$new[0]{end} if ($old_end < $$new[0]{end});
	}
	return $gff_cluster;
}

############################################################

1;

__END__

# http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo?revision=1.263
# http://www.sequenceontology.org/resources/gff3.html
