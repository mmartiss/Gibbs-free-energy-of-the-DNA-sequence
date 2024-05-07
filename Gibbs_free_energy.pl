#!/usr/bin/perl
use strict;
use warnings;

=head
M. A.
This program reads DNA sequences in FASTA format.
Calculates the free energy of duplex formation for the 
given sequence and its complementary sequence.
It uses the SantaLuca algorithm.
=cut

my $parameters_ref;

my %parameters;

$parameters_ref = \%parameters;

sub pairEnergy($$){
    #The subroutine first tries to find the energy of the pair 
    #in the 'pairwise' part of the $parameters hash. 
    #If it can't find the energy, it calculates the complement of the pair.
    #IT tries to find the energy of the complement pair.

    my ($arg1, $parameters_ref) = @_;
    my $energy = $parameters_ref->{'pairwise'}->{$arg1};

    if (!defined $energy) {
        my %complement = ('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C');

        my $complement_pair = $complement{substr($arg1, 1, 1)} . 
            $complement{substr($arg1, 0, 1)};

        $energy = $parameters_ref->{'pairwise'}->{$complement_pair};
    }

    return $energy;
}

sub calculateSantaLuciaScore($$) {
    # The subroutine iterates over the DNA sequence, extracting each pair 
    #of nucleotides and calculating their energy using the pairEnergy subroutine.
    # It adds the energy of each pair to a total energy counter.
    my ($arg1, $parameters_ref) = @_;
    my $total_free_energy = 0;
    for (my $i = 0; $i < length($arg1) - 1; $i++) {
        my $pair = substr($arg1, $i, 2);
        $total_free_energy += pairEnergy($pair, $parameters_ref);
    }
    $total_free_energy += $parameters_ref->{'terminal'}->{substr($arg1, 0, 1)};
    $total_free_energy += $parameters_ref->{'terminal'}->{substr($arg1, -1 , 1)};
    
    #If you want to use getTerminalNucleotideScore, uncomment down below
    #my $first_nucleotide_energy = getTerminalNucleotideScore(substr($arg1, 0 , 1), $parameters);
    #printf $first_nucleotide_energy . "\n";
    #my $last_nucleotide_energy = getTerminalNucleotideScore(substr($arg1, -1 , 1), $parameters);
    #printf $last_nucleotide_energy . "\n";

    my $symmetry_correction = isSelfComplementarySequence($arg1) ? 0.43 : 0;
    $total_free_energy += $symmetry_correction;
    return $total_free_energy;
}

sub getTerminalNucleotideScore($;$) {
    # The subroutine returns given nucleotides energy
    my ($nucleotide, $parameters_ref) = @_;
    my $nucleotide_free_energy = $parameters_ref->{'terminal'}->{$nucleotide};
    return $nucleotide_free_energy;
}

sub printResults($$) {
    # The subroutine prints out sequence total free energy with given format
    my ($arg1, $arg2) = @_;
    printf ">%s |dG:%.2f|\n", $arg1, sprintf("%.2f", $arg2);
}

sub isSelfComplementarySequence($) {
    # The subroutine checks if the sequence is the same as its reverse complement.
    my ($arg1) = @_;
    my $reverse_complement = reverse($arg1);
    $reverse_complement =~ tr/ACGT/TGCA/;
    return $arg1 eq $reverse_complement ? 1 : 0;
}

%parameters = (
    'terminal' => { 
        'A' => 1.03, 
        'T' => 1.03, 
        'G' => 0.98, 
        'C' => 0.98 
    },
    'pairwise' => {
        'AA' => -1.00, #'TT' => -1.00,
        'AT' => -0.88, #'TA' => -0.88,
        'TA' => -0.58, #'AT' => -0.58,
        'CA' => -1.45, #'GT' => -1.45,
        'TG' => -1.45, #'AC' => -1.45,
        'GT' => -1.44, #'CA' => -1.44,
        'CT' => -1.28, #'GA' => -1.28,
        'GA' => -1.30,
        'CG' => -2.17, #'GC' => -2.17,
        'GC' => -2.24,
        'GG' => -1.84 #'CC' => -1.84,
    }
);

foreach my $input (@ARGV) {
    if (-e $input) { #Checks if file exists
    }
    else {
        die "Input '$input' is not a valid file\n";
    }
}


# Main
my $line_number = 0;
my $header = '';
my $sequence = '';
my $previous_line_was_header = 0;
while (my $line = <>) {
    $line_number++;
    chomp $line;
    if ($line =~ /^>(.*)/) {
        if ($previous_line_was_header) {
            $line_number--;
            warn "$0: WARNING, empty sequence found, line $line_number -- " .
                "Skipping sequence.\n";
            $header = $1;
            $line_number++;
        } 
        else {
            if ($sequence ne '') {
                my $free_energy = calculateSantaLuciaScore(
                    $sequence, $parameters_ref
                );
                printResults($header, $free_energy);
            }
            $header = $1;
            $sequence = '';
        }
        $previous_line_was_header = 1;
    } 
    else {
        $sequence .= uc($line);
        $previous_line_was_header = 0;
    }
    if (eof) {
        $line_number = 0;
    }
}
if ($sequence ne '') {
    my $free_energy = calculateSantaLuciaScore(
        $sequence, $parameters_ref
    );
    printResults($header, $free_energy);
}

foreach my $input (@ARGV){
    close $input or die "Failed to close file: $!"; #$! - Contains error message
}
