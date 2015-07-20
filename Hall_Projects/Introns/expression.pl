use warnings;
use strict;

## Check arguments.
die qq[Usage: perl $0 <pattern-file> <input-file>\n] unless @ARGV == 2;

## Open input files.
open my $pattern_fh, qq[<], shift @ARGV or die qq[Cannot open pattern file\n];
open my $input_fh, qq[<], shift @ARGV or die qq[Cannot open input file\n];

## Hash to save patterns.
my (%pattern, %input);

## Read each pattern and save how many times appear in the file.
while ( <$pattern_fh> ) { 
    chomp;
    if ( exists $pattern{ $_ } ) { 
        $pattern{ $_ }->[1]++;
    }   
    else {
        $pattern{ $_ } = [ $., 1 ];
    }   
}

## Read file with data and save them in another hash.
while ( <$input_fh> ) { 
    chomp;
    my @f = split;
    $input{ $f[1] } = $_; 
}

## For each pattern, search it in the data file. If it appears, print line those
## many times saved previously, otherwise print a blank line.
for my $p ( sort { $pattern{ $a }->[0] <=> $pattern{ $b }->[0] } keys %pattern ) { 
    if ( $input{ $p } ) { 
        printf qq[%s\n], $input{ $p } for ( 1 .. $pattern{ $p }->[1] );
    }   
    else {
         # Old behaviour.
         # printf qq[\n];

         # New requirement.
         printf qq[\n] for ( 1 .. $pattern{ $p }->[1] );
    }   
}
