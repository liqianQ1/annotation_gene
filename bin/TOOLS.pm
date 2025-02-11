package TOOLS;
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;

############################################################

sub absolute_path {
	my ($path, $pwd) = @_;
	my ($abs_path, $name, $parent_path);

	$pwd ||= `pwd`;
	chomp ($path, $pwd);
	$pwd  =~ /^\// or die "Error: $pwd is not an absolute path";
	$path =~ /^\// or $path = "$pwd/$path";

	my @old_path = split /\/+/, $path;
	my @new_path;
	foreach (@old_path) {
		if    ($_ eq '' || $_ eq '.') {next;}
		elsif ($_ eq '..')            {pop  @new_path;}
		else                          {push @new_path, "$_";}
	}

	$abs_path = "/".(join '/', @new_path);
	$name = pop @new_path or die "Error: $path";
	$parent_path = "/".(join '/', @new_path);

	return $abs_path, $name, $parent_path;
}

############################################################

sub read_table {
	my $iput = shift;
	my ($key, $tab);
	open IPUT, "<$iput" or die;
	my $txt = <IPUT>;
	$key = read_title ($txt);
	while (<IPUT>) {
		/^#/ and next;
		/\w/ or next;
		push @$tab, read_value ($key, $_);
	}
	close IPUT;
	return $key, $tab;
}


sub read_title {
	my $txt = shift;
	my $key;
	chomp $txt;
	$txt =~ s/ +//g;
	$txt =~ s/^\#//;
	@$key = split "\t", $txt;
	return $key;
}


sub read_value {
	my ($key, $txt) = @_;
	my $hash;
	chomp $txt;
#	$txt =~ s/ //g;
	my @val = split "\t", $txt;
	for (my $i=0; $i<@{$key}; $i++) {
		$$hash{$$key[$i]} = (defined $val[$i])? $val[$i]:'';
	}
	return $hash;
}


sub read_hash {
	my ($iput, $hash) = @_;
	open IPUT, "<$iput" or die;
	my $txt;
	while ($txt=<IPUT>) {
		$txt =~ /^#/ or next;
		my $key = read_title ($txt);
		$txt = <IPUT>;
		my $tab = read_value ($key, $txt);
		$$hash{$_} = $$tab{$_} foreach (keys %$tab);
	}
	close IPUT;
}

############################################################

sub print_table {
	my ($oput, $key, $tab) = @_;
	open OPUT, ">$oput" or die;
	print OPUT (print_title ($key));
	if (ref $tab eq 'ARRAY') {
		foreach (@{$tab}) {
			print OPUT (print_value ($key, $_));
		}
	}
	else {
			print OPUT (print_value ($key, $tab));
	}
	close OPUT;
}


sub print_title {
	my $key = shift;
	my $txt = '#'.(join "\t", @$key)."\n";
	return $txt;
}


sub print_value {
	my ($key, $hash) = @_;
	my @val = ();
	foreach (@{$key}) {
		if (defined $$hash{$_}) {
			push @val, $$hash{$_};
		}
		else {
			push @val, '';
		}
	}
	my $txt = (join "\t", @val)."\n";
#	$txt =~ s/ //g;
	return $txt;
}

############################################################

sub read_ini {
	my $iput = shift;
	my $config;
	open IPUT, "<$iput" or die;
	my $key;
	while (<IPUT>) {
		/^#/ and next;
		/\w/ or next;
		if (/^(\S+)\s*=\s*(.*?)\s+(\#.*)?$/) {
			$$config{$1} = $2;
		}
		elsif (/^(\S+):/) {
			$key = $1;
		}
		elsif (/^\t(.*?)\s+(\#.*)?$/){
			push @{$$config{$key}}, $1;
		}
		else {
			die;
		}
	}
	close IPUT;
	return $config;
}


sub cmd {
	my @c = @_;
	my $config = shift @c;
	if (ref ($$config{$c[0]}) eq 'ARRAY') {
		my @cmd;
		foreach (@{$$config{$c[0]}}) {
			while ($_ =~ s/\$([A-Z][A-Z_\d]+)/$$config{$1}/g){}
			my $cmd = $_;
			$cmd =~ s/\$(\d+)/$c[$1]/g;
			push @cmd, $cmd;
		}
		return \@cmd;
	}
	else {
		while ($$config{$c[0]} =~ s/\$([A-Z][A-Z_\d]+)/$$config{$1}/g){}
		my $cmd = $$config{$c[0]};
		$cmd =~ s/\$(\d+)/$c[$1]/g;
		return $cmd;
	}
}

############################################################

sub fq2fq {
	my ($iput_fq, $oput_fq, $optr, $opts, $optl) = @_;
	$optr ||= 0;
	$opts ||= 0;
	$optl ||= 100 if ($opts);
	if ($iput_fq=~/gz$/) {open IPUT, "gzip -cd $iput_fq |" or die;}
	else                 {open IPUT, "<$iput_fq" or die;}
	open OPUT, ">$oput_fq" or die;
	my $i=0;
	while (my $id=<IPUT>) {
		my $seq = <IPUT>;
		<IPUT>;
		my $qvl = <IPUT>;
		chomp $seq;
		chomp $qvl;
		$i++;
		if ($opts or $optl) {
			$seq = substr $seq, $opts, $optl;
			$qvl = substr $qvl, $opts, $optl;
		}
		print OPUT "$id$seq\n+\n$qvl\n";
		last if ($i==$optr);
	}
	close IPUT;
	close OPUT;
}


sub fq2fa {
	my ($iput_fq, $oput_fa, $optr, $opts, $optl) = @_;
	$optr ||= 0;
	$opts ||= 0;
	$optl ||= 100 if ($opts);
	if ($iput_fq=~/gz$/) {open IPUT, "gzip -cd $iput_fq |" or die;}
	else                 {open IPUT, "<$iput_fq" or die;}
	open OPUT, ">$oput_fa" or die;
	my $i=0;
	while (my $id=<IPUT>) {
		$id =~ s/^\@/\>/;
		my $seq = <IPUT>;
		<IPUT>;
		my $qvl = <IPUT>;
		chomp $seq;
		$i++;
		if ($opts or $optl) {
			$seq = substr $seq, $opts, $optl;
		}
		print OPUT "$id$seq\n";
		last if ($i==$optr);
	}
	close IPUT;
	close OPUT;
}

############################################################

1;
__END__
