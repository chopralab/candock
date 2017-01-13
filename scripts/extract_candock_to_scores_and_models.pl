#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $number=10;
my $no_make_files = 0;
my @ligands = ();

GetOptions("n|number=i" => \$number,
           "l|ligand=s" => \@ligands,
           "f|nofiles"  => \$no_make_files);

my $currkey;
my $currlig;
my %poses = ();

print "number, molecule, score\n";

while (my $line = <>) {
    if ( $line =~ /^REMARK.*COMPLEX OF (\S*).* SCORE OF (\S*)/) {
        $currkey = $.;
        $currlig = $1;

        if ( $#ligands != -1 and not ($currlig ~~ @ligands) ) {
                next;
        }

        unless ( defined $poses{$currlig} ) {
            $poses{$currlig} = {};
        }
        
        $poses{$currlig}->{$currkey} = {PROTEIN=>[], LIGAND=>[], SCORE=>$2};
        <>;
    }

    if ( $#ligands != -1 and not ($currlig ~~ @ligands) ) {
            next;
    }
    
    elsif ( $line =~ /^ATOM  / ) {
        push @{$poses{$currlig}->{$currkey}->{PROTEIN}}, substr( $line, 0, 80 );
    }

    elsif ( $line =~ /^HETATM/ ) {
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

      print "$count, $ligandkey, $ligand_hash->{$key}->{SCORE}\n";

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
