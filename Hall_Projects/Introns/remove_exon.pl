#!/bin/perl
#use strict;
#use warnings;

open(DNAFILE, "<", "Phase0Seq.txt");
my @sequence=<DNAFILE>;
close DNAFILE;
open(INDEXFILE, "<", "Phase0Splice5pri.txt");
my @intronStart=<INDEXFILE>;
close INDEXFILE;
open (FILE, ">Frame0NoExon.txt");

my $j;
my $i;
my $observed;

#Adjust j for number of sequences
for($j=0;$j<27000;$j+=1)
{
for($i=0)
{
$observed=substr($sequence[$j],$i+$intronStart[$j],-1);
    #Unhash below to shift for phase and frame of interest for SCindex.pl
    
    #$observed=substr($sequence[$j],$i+$intronStart[$j]+1,-1);
    #$observed=substr($sequence[$j],$i+$intronStart[$j]+2,-1);
}
print FILE "$observed \n";
}
close FILE;
