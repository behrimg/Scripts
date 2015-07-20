Description of scripts used for Intron project

--------------
expression.pl

perl script to match affymetrix expression data to genes of interest. 
pattern file = genes (or string) of interest
input file = array (or data) files of interest

pattern file is matched to the 2nd column of tab delimited data 
(the gene id column of array data)

if pattern is matched entire line of array data is printed
if pattern is not matched a blank line is printed so data files 
can easily be pasted together

------------------------
intron_index.pl

perl script that determines the position of the first 5' splice site encountered using 
data output from feature_extract. http://www.cbs.dtu.dk/services/FeatureExtract/

script must be altered to update number of sequences to analyze

------------
remove_exon.pl

perl script that uses the output from intron_index.pl to remove the first exon from 
sequences in preparation for SCIndex.pl

Script is currently written to look at phase 0 frame 0 stop codons, but can be adjusted 
for different phase/frame combinations by changing the hashed $observed = line. 

script must be altered to update number of sequences to analyze

--------
SCIndex.pl

perl script to determine the position of the first inframe stopcodon encountered 
then print that position followed by the stop codon sequence +1

frame is determined by input file from remove_exon.pl

script must be altered to update number of sequences to analyze
__________