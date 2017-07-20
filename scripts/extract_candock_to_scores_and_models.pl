#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $number=10;
my $no_make_files = 0;
my $print_title = 0;
my %variables = ();

GetOptions("n|number=i" => \$number,
           "f|nofiles"  => \$no_make_files,
           "t|title" => \$print_title,
           "v|variables=s" => \%variables);

my $currkey;
my $currlig;
my %poses = ();

if ($print_title) {

    print 'number, molecule, score';

    foreach my $variable (keys %variables) {
        print ", $variable";
    }

    print "\n";

}

while (my $line = <>) {
    if ( $line =~ /REMARK  20 non-binder (\S+)/) {
        $currlig = $1;
        $currkey = "non";

        unless ( defined $poses{$currlig} ) {
            $poses{"non"} = {};
        }

        $poses{$currlig}->{"non"} = {PROTEIN=>[], LIGAND=>[], SCORE=>10000};
    }

    if ( $line =~ /^REMARK.*COMPLEX OF (\S*).* SCORE OF (\S*)/) {
        $currkey = $.;
        $currlig = $1;

        unless ( defined $poses{$currlig} ) {
            $poses{$currlig} = {};
        }
        
        $poses{$currlig}->{$currkey} = {PROTEIN=>[], LIGAND=>[], SCORE=>$2};
        <>;
    }

    elsif ( ! $no_make_files && $line =~ /^ATOM  / ) {
        push @{$poses{$currlig}->{$currkey}->{PROTEIN}}, substr( $line, 0, 80 );
    }

    elsif ( ! $no_make_files && $line =~ /^HETATM/ ) {
        push @{$poses{$currlig}->{$currkey}->{LIGAND}}, substr( $line, 0, 80 );
    }
}

my $CREATE_PSE = undef;

unless ( $no_make_files ) {
    open $CREATE_PSE, '>', 'create_pse.py' or die $!;
}

my @ordered_ligands = ();

while ( my ($ligandkey, $ligand_hash) = each %poses ) {

    my $count = 1;

    foreach my $key ( sort  {$ligand_hash->{$a}->{SCORE} <=> $ligand_hash->{$b}->{SCORE}}
                    ( keys %{$ligand_hash} ) ) {

        unless ( $no_make_files ) {
            open my $COMBINED, '>', "$ligandkey"."_combo_$count.pdb" || die $!;

            print $COMBINED ($_, "\n") foreach ( @{$ligand_hash->{$key}->{PROTEIN}} );
            print $COMBINED ($_, "\n") foreach ( @{$ligand_hash->{$key}->{LIGAND}} );

            close $COMBINED;

            my $pymol_key = $ligandkey;
            $pymol_key =~ s/[\(|\)]/_/g;

            print $CREATE_PSE "cmd.load( \"" .$ligandkey. "_combo_$count.pdb\")\n";
            print $CREATE_PSE "cmd.create( \"$pymol_key\", \"" .$pymol_key."_combo_$count\", 1,";
            $count == 1 ? print $CREATE_PSE " 0)\n" : print $CREATE_PSE " -1)\n";
            print $CREATE_PSE "cmd.delete( \"".$pymol_key."_combo_$count\" )\n";

            push @ordered_ligands, { ligand_name => $pymol_key, score => $ligand_hash->{$key}->{SCORE} };
        }

        print "$count, $ligandkey, $ligand_hash->{$key}->{SCORE}";
      
        foreach my $variable (values %variables ) {
            print ", $variable"
        }

        print "\n";

        $count++;

        last if $count > $number;
    }

}

unless ( $no_make_files ) {

    foreach my $sorted_ligands (  sort {$b->{score} <=> $a->{score}  } 
                              ( @ordered_ligands)  ) {
        print $CREATE_PSE "cmd.order( \"$sorted_ligands->{ligand_name}\", location=\"top\" )\n";

    }

    print $CREATE_PSE "cmd.hide(\"everything\")\n";
    print $CREATE_PSE "cmd.show(\"cartoon\")\n";
    print $CREATE_PSE "cmd.show(\"sticks\",\"organic\")\n";

    print $CREATE_PSE "cmd.save( \"all.pse\")\n";

    close $CREATE_PSE;
}
