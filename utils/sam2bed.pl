#!/bin/env perl

use 5.012;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::samParser;
use pm::gpeParser;
use pm::bedParser;

my ($tags, $getXS);
GetOptions(
            't|tags=s'  => \$tags,
            's|XS'      => \$getXS,
            'h|help'    => sub{usage()}
        ) || usage();

if(defined $ARGV[0]){
    if (-B $ARGV[0]){
        open IN, "samtools view -h $ARGV[0]|" or die "Can't open $ARGV[0]: $!";
    }else{
        open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";
    }
}else{
    open IN, "-";
}

while(<IN>){
    next if /^@/;
    unless(&samParser::isMapped($_)){
        say STDERR $_;
        next;
    }
    chomp;
    my @fields = split "\t";
    my ($name, $flag, $chr, $start, $mapQ, $cigar) = @fields[0..5];
    my $strand = $flag & 0b10000 ? '-' : '+';
    if(defined $getXS){
        my $xs = &samParser::getTagValue($_, "XS");
        if(defined $xs && $xs =~ /[+-]/){
            $strand = $xs;
        }
    }
    $start--;
    my ($blockStarts, $blockEnds) = &samParser::cigar2Blocks($start, $cigar);
    my @blockSizes = &gpeParser::getSizes($blockStarts, $blockEnds);
    my @blockRelStarts = &bedParser::getRelStarts($blockStarts);
    print join "\t", ($chr, $blockStarts->[0], $blockEnds->[-1], $name, $mapQ, $strand, $blockStarts->[0], $blockEnds->[-1], '255,0,0',
                    (scalar @$blockStarts),
                    (join ',', @blockSizes),
                    (join ',', @blockRelStarts)
                    );
    if(defined $tags){
        for my $tag(split ',', $tags){
            my $tagValue = &samParser::getTagValue($_, $tag);
            print "\t", defined $tagValue ? $tagValue : 'NA';
        }
    }
    print "\n";
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.sam/bam >OUTPUT.bed12+ 2>unmapped.sam
    If INPUT isn't specified, input from STDIN
Option:
    -t --tags   STR[]   Comma-separated tag names, eg. MD,XS
                        Values of these tags will be append to additional columns in output in order
    -s --XS             Get the strand from XS tag. Default: the strand inferred from flag
    -h --help           Print this help information
HELP
    exit(-1);
}