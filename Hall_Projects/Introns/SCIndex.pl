#!/usr/bin/perl

use strict; 
use warnings;
use autodie;

die "Usage: $0 Filename\n" if @ARGV != 1;

my $file = shift;

open my $infh, '<', $file;
open my $outfh, '>', 'AtPTCindex.txt';

while (my $line = <$infh>) {
    chomp($line);

    my $result = '';

    for (my $i = 0; $i < (length($line) - 2); $i += 3) {
        my $codon = substr($line, $i, 3);
        if ($codon =~ m/TAG|TGA|TAA/) {
            # m added before /TAG... above
            my $stopseq = substr($line, $i, 4);
               $result .= "$i,$stopseq";
            last;
        }
    }

    print "$result\n";
#   print $outfh "$result\n";
#   print $outfh "$result $.\n";
}

close $infh;
close $outfh;