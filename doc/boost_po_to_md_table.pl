#!/usr/bin/env perl

use warnings;
use strict;

my %variable_group = ();
my $current_variable;
my $current_option;

while (<>) {

	chomp;

	if ( $_ =~ m/^$/ ) {
		next;
	}

	if ( $_ =~ m/^\S/ ) {
		$current_variable = $_;
		$variable_group{ $_ } = {};
		next;
	}

        if ( $current_variable eq "Generic options:" ) {
		next;
	}

        if ( $_ =~ m/^\s+\-\-(\S+)\s\[=arg\(=(\S+)\)\](?:\s+(\S.*))?/x ) {
                $current_option = $1;
                $variable_group{ $current_variable }->{$current_option} = {def_value=>$2.' (Implicit)', desc=>$3};
                next;
        }

        if ( $_ =~ m/^\s+\-\-(\S+)\sarg\s\(=(\S+)\)(?:\s+(\S.*))?/x ) {
                $current_option = $1;
                $variable_group{ $current_variable }->{$current_option} = {def_value=>$2, desc=>$3};
		next;
        }

	if ( $_ =~ m/^\s+\-\-(\S+)\sarg\s+(\S.*)/x ) {
		$current_option = $1;
		$variable_group{ $current_variable }->{$current_option} = {def_value=>'~None~', desc=>$2};
		next;
	}

	if ( $_ =~ m/^\s+(\S.*)/) {
		$variable_group{ $current_variable }->{$current_option}->{desc}.=$1;
	}
}

while ( (my $key, my $value) = each %variable_group ) {

	next if $key eq "Generic options:";
	
	print $key, "\n\n";

	print "|Option Name | Default | Description |\n";
	print "|------------|--------| ----------------------------- |\n";
	
	while ( (my $key2, my $value2) = each %{$value} ) {
		print "| $key2 | $value2->{def_value} | $value2->{desc} |\n"
	}

	print "\n";
}

