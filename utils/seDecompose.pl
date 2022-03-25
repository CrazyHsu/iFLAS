#!/bin/env perl

use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Options:
    -h  --help       Print this help information
HELP
    exit(-1);
}

GetOptions(
            'h|help' => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
while(<IN>){
    chomp;
    my ($chr, $gene, $locus, $strand) = split ":";
    my ($leftSite, $SEs, $rightSite) = split '@', $locus;
    my @SEs = split ';', $SEs;
    my ($firstSeStart, $firstSeEnd) = split '-', $SEs[0];
    my ($lastSeStart, $lastSeEnd) = split '-', $SEs[-1];
    say join "\t", ($chr, $leftSite, $firstSeStart, $strand);
    for(my $i = 1; $i <= $#SEs; $i++){
        my $preSeEnd = (split '-', $SEs[$i-1])[1];
        my $seStart = (split '-', $SEs[$i])[0];
        say join "\t", ($chr, $preSeEnd, $seStart, $strand);
    }
    say join "\t", ($chr, $lastSeEnd, $rightSite, $strand);
    say STDERR join "\t", ($chr, $leftSite, $rightSite, $strand); # exclusion
}

