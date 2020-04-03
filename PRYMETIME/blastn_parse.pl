#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO;

#print usage statement if blast output filename is not provided on teh command line
my $usage = "\nUSAGE: $0 blast_file\n\n";
my $file = shift or die ($usage);

#Import Blast output file as a BioPerl object
my $SearchIO_obj = new Bio::SearchIO(-format => 'blast',-file => $file);

#Print header line of the output
print "Query\tHit\tQuery_Start\tQuery_End\tHit_Start\tHit_End\tStrand\tLength\tIdentity\tEvalue\n";

#loop through the blast file, going through each query ("Result"), each hit within each query, and each high-scoring pair (HSP) within each hit.
#Extraction and print key information including hit location, length, identity, and e-value.
while( my $result_obj = $SearchIO_obj->next_result ) {
        my $query_name = $result_obj->query_name;
        my $query_desc = $result_obj->query_description;
        while( my $hit_obj = $result_obj->next_hit ) {
                my $hit_name = $hit_obj->name;
                my $hit_desc = $hit_obj->description;
                while (my $hsp_obj = $hit_obj->next_hsp){
                        my $evalue = $hsp_obj->evalue;
                        my $id = $hsp_obj->percent_identity;
                        my $length = $hsp_obj->length('total');
                        my $query_start = $hsp_obj->start('query');
                        my $query_end = $hsp_obj->end('query');
                        my $hit_start = $hsp_obj->start('hit');
                        my $hit_end = $hsp_obj->end('hit');
                        my $strand = $hsp_obj->strand('hit');
                        print "$query_name $query_desc\t$hit_name $hit_desc\t$query_start\t$query_end\t$hit_start\t$hit_end\t$strand\t$length\t$id\t$evalue\n"
                }
        }
}


exit;
