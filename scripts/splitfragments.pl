#!/usr/bin/env perl

use strict;
use warnings;

my %atoms = ();
my %fragments = ();

while (<>) {
	if ( $_ =~ m/^MODEL\s+(\d+)/x ) {
		%atoms = ();
	} elsif ( $_ =~ m/^HETATM\s*(\d+)/x ) {
		$atoms{$1} = $_;
	} elsif ( $_ =~ m/^REMARK\s+8\sRIGID\s(\d+)/x) {

		next if ( defined ($fragments{$1}) || $1 == -1);
		
		$fragments{$1} = [];

		chomp;

		my @frag_atoms = split /\s[c|j]/, $_;
		my $junk = shift @frag_atoms;
		foreach my $i (@frag_atoms) {
			push @{$fragments{$1}}, $atoms{$i};
		}
	}
}

print "REMARK   4 NRPDB\n";

foreach my $frag_num ( sort{ $a <=> $b } (keys %fragments) ) {
	printf "MODEL %8d\n", $frag_num;
	
	my @atoms_in_frag = @{$fragments{$frag_num}};

	foreach my $atom_num (0 .. $#atoms_in_frag) {
		printf "HETATM%5d", $atom_num + 1;
		print  substr ($atoms_in_frag[$atom_num], 11, 69), "\n";
	}

	print "ENDMDL\n";
}

print "REMARK   4 END\n";

