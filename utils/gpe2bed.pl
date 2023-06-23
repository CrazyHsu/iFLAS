#!/usr/bin/env perl
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::gpeParser;

my ($bin, $gene, $plus);
my ($bedType, $itemRgb) = (12, "0,0,0");

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.gpe >OUTPUT.bed
    if INPUT.gpe isn't specified, input from STDIN
Options:
    -b --bin            With bin column
    -t --bedType    INT Bed type. It can be 3, 6, 9 or 12[$bedType]
    -i --itemRgb    STR RGB color[$itemRgb]
    -g --gene           Output 'gene name' in INPUT.gpe as bed plus column
    -p --plus           Output bed plus when there are additional columns in gpe
    -h --help           Print this help information
HELP
    exit(-1);
}

GetOptions(
            'b|bin'         => \$bin,
            't|bedType=i'   => \$bedType,
            'i|itemRgb=s'   => \$itemRgb,
            'g|gene'        => \$gene,
            'p|plus'        => \$plus,
            'h|help'        => sub{usage()}
        )||usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN,"$ARGV[0]" or die "Can't open file ($ARGV[0]): $!";

while(<IN>){
    chomp;
    my @fields = split "\t";
    shift @fields if defined $bin;
    if($bedType == 3){
        print join "\t", @fields[1, 3, 4];
    }elsif($bedType == 6){
        print join "\t", @fields[1, 3, 4, 0, 10, 2];
    }elsif($bedType == 9){
        print join "\t",(@fields[1, 3, 4, 0, 10, 2, 5, 6], $itemRgb);
    }elsif($bedType == 12){
        my @blockStarts = split ",", $fields[8];
        my @blockEnds = split ",", $fields[9];
        my @blockSizes = gpeParser::getSizes(\@blockStarts, \@blockEnds);
        my @relativeStarts = map{$blockStarts[$_] - $blockStarts[0]}0..$#blockStarts;
        print join "\t", (@fields[1, 3, 4, 0, 10, 2, 5, 6],$ itemRgb, $fields[7], (join ",", @blockSizes) , (join ",", @relativeStarts));
    }
    if(defined $gene){
        print "\t$fields[11]";
    }
    if(defined $plus && @fields > 15){
        print "\t";
        print join "\t", @fields[15..$#fields];
    }
    print "\n";
}
