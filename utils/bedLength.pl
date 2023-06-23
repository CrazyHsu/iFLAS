#!/usr/bin/env perl

use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::bedParser;
use pm::gpeParser;

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Options:
    -h --help       Print this help information
HELP
    exit(-1);
}

GetOptions(
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
while(<IN>){
    chomp;
    my ($chr, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, undef, $blockCount, $sizes, $relStarts) = split "\t";
    my $length;
    if(defined $relStarts){
        my @sizes = split ',', $sizes;
        my @relStarts = split ',', $relStarts;
        my ($starts, $ends) = &bedParser::getAbsLoc($start, \@sizes, \@relStarts);
        $length = &gpeParser::getExonsLength($starts, $ends);
    }else{
        $length = $end - $start;
    }
    say join "\t", ($_, $length);
}

