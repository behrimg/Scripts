#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;

open my $PromoterFile, '<', 'PromoterNew.uniq.txt'     or die $!;
open my $SNPSFile,     '<', 'SNP_Position_File.txt'        or die $!;
open my $output,       ">", "PromoterMatch.txt" or die $!;


my @data;

while (<$PromoterFile>) {
    chomp;
    my @CordFile    = split;
    my $Lend        = $CordFile[0];
    my $Rend        = $CordFile[1];
    my $Gene		= $CordFile[2];
    my $Sigma24 	= $CordFile[6];
    my $Sigma28		= $CordFile[7];
    my $Sigma32		= $CordFile[8];
    my $Sigma38		= $CordFile[9];
    my $Sigma70		= $CordFile[10];

    push(
        @data,
        {   lend        => $CordFile[0],
            rend        => $CordFile[1],
            gene		=> $CordFile[2],
            sigma24  	=> $CordFile[6],
            sigma28  	=> $CordFile[7],
            sigma32  	=> $CordFile[8],
            sigma38  	=> $CordFile[9],
            sigma70  	=> $CordFile[10]
        }
    );
}

print Dumper \@data;

foreach my $value (<$SNPSFile>) {
    chomp $value;
    my $found = 0;
    foreach my $element (@data) {
        if (    $value >= $element->{lend}
            and $value <= $element->{rend} )
        {
            #print "Found $value\n";
            print {$output} join( "\t",
                $value, $element->{gene}, $element->{sigma24}, $element->{sigma28}, $element->{sigma32}, $element->{sigma38}, $element->{sigma70} ),
                "\n";
            $found++;
            last;
        }

    }
    if ( not $found ) {
        print {$output} $value,"\n";
    }
}

close $output;
close $PromoterFile;
close $SNPSFile;
