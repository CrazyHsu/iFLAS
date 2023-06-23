#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::gpeParser;
use pm::bedParser;

my ($bed, $geneCol) = (12);
sub usage{
    my $scriptName = basename $0;
    print <<HELP;
Usage: perl $scriptName INPUT.bed >OUTPUT.gpe
    If INPUT.bed not specified, input from STDIN
Options:
    -b --bed   INT  The bed format([12], 6, 8)
    -g --gene  INT  Set the 'gene name' in OUTPUT.gpe with the INT column in INPUT.bed.
                    'Last' can be used specified to the last column in INPUT.bed.
    -h --help       Print this help information
HELP
    exit(-1);
}

GetOptions(
            'b|bed=i'   => \$bed,
            'g|gene=s'  => \$geneCol,
            'h|help'    => sub{&usage()}
          )||usage();

$ARGV[0]='-' unless defined $ARGV[0];
open BED,"$ARGV[0]" or die "Can't open file ($ARGV[0]): $!";

while(<BED>){
    chomp;
    my @fields = split "\t";
    my ($chr, $start, $end, $name, $strand, $cdsStart, $cdsEnd, $sizes, $relStarts) = @fields[0..3, 5..7, 10, 11];
    my ($blockStarts, $blockEnds, $exonFrames);
    if($bed == 6){
        $blockStarts = [$start];
        $blockEnds = [$end];
        $exonFrames = 'unk';
    }elsif($bed == 8){
        $blockStarts = [$start];
        $blockEnds = [$end];
        $exonFrames = join ',', @{&gpeParser::getExonFrames($cdsStart, $cdsEnd, $strand, (join ',', @$blockStarts), (join ',', @$blockEnds))};
    }elsif($bed == 12){
        my @sizes = split ',', $sizes;
        my @relStarts = split ',', $relStarts;
        ($blockStarts, $blockEnds) = bedParser::getAbsLoc($start, \@sizes, \@relStarts);
        $exonFrames = join ',', @{&gpeParser::getExonFrames($cdsStart, $cdsEnd, $strand, (join ',', @$blockStarts), (join ',', @$blockEnds))};
    }
    
    my $geneName;
    if(defined $geneCol){
        $geneName = $geneCol eq "Last" ? $fields[$#fields] : $fields[$geneCol-1];
    }else{
        $geneName = 'gene name';
    }
    say join "\t",($name,
                   $chr,
                   $strand,
                   $blockStarts->[0],
                   $blockEnds->[-1],
                   $cdsStart,
                   $cdsEnd,
                   scalar (@$blockStarts),
                   (join ',', @$blockStarts),
                   (join ',', @$blockEnds),
                   0,
                   $geneName,
                   "unk",
                   "unk",
                   $exonFrames);
}
