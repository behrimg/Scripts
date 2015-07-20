#!/usr/bin/perl

#use strict; 
use warnings;

# A program to find the first inframe stop codon of non-spliced intron containing genes
print "ENTER THE FILENAME FOR DNA SEQUENCES:= ";
# Asks for Sequence file and if file does not exist prints error message
my $filename = <STDIN>;
#my $sequence;
my @sequence;
chomp $filename;
unless ( open(DNAFILE, $filename) ) {
    print "Cannot open file \"$filename\"\n\n";
}
@sequence = <DNAFILE>;
close DNAFILE;
open (FILE, ">AltCeInt5Pri.txt");
my $j;
#Change $j<(375) to n=number of sequences
for($j=0;$j<(23289);$j+=1)
{
my $string = $sequence[$j];
my $char = D;
my $result = index($string,$char);
print FILE "$result\n" 
}
#}
close FILE;
exit;

#Giving last inframe stop codon 

