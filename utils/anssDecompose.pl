#!/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

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
    -n  --anss  STR 5 for A5SS, 3 for A3SS
    -h  --help      Print this help information
HELP
    exit(-1);
}

my ($anss);
GetOptions(
            'n|anss=s'  => \$anss,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
while(<IN>){
    chomp;
    my ($chr, $gene, $locus, $strand) = split ":";
    my ($start, $end) = split '-', $locus;
    my ($excStart, $excEnd, $incStart, $incEnd);
    if($anss == 5 && $strand eq '+' || $anss == 3 && $strand eq '-'){
        say join "\t", ($chr, $end, "dummy", $strand);
        say STDERR join "\t", ($chr, $start, "dummy", $strand); # exclusion
    }else{
        say join "\t", ($chr, "dummy", $start, $strand);
        say STDERR join "\t", ($chr, "dummy", $end, $strand); # exclusion
    }
}

