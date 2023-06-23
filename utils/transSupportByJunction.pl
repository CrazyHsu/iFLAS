#!/bin/env perl

use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::gpeParser;
use pm::bedParser;

my ($junctionFile, $bin);

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.gpe >OUTPUT.tsv
    If INPUT.gpe isn't specified, input from STDIN
Option:
    -j  --junction  FILE    The junction file in bed
    -b  --bin               With bin in gpe
    -h  --help              Print this help information
Output:
    STDOUT: tsv file with
        Col1: transcript name
        Col2: intron count of the transcript
        Col3: supported intron count of the transcript
    STDERR: tsv file with
        Col1: intron count
        Col2: supported intron count
        Col3: count of transcripts with Col1 introns
        Col4: count of transcripts with Col2 supported introns
HELP
    exit(-1);
}

GetOptions(
            'j|junc=s'  => \$junctionFile,
            'b|bin'     => \$bin,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
open JUNC, "$junctionFile" or die "Can't read file ($junctionFile): $!";

my %juncHash;
while(<JUNC>){
    chomp;
    my ($chr, $start, $strand, $blockSizes, $blockStarts) = (split "\t")[0, 1, 5, 10, 11];
    my @blockSizes = split ',', $blockSizes;
    my @blockStarts = split ',', $blockStarts;
    my ($starts, $ends) = &gpeParser::getIntrons(&bedParser::getAbsLoc($start, \@blockSizes, \@blockStarts));
    $juncHash{"$chr:$strand:$starts->[0]-$ends->[0]"} = '';
}

my %summaryHash;
while(<IN>){
    chomp;
    my @fields = split "\t";
    shift @fields if defined $bin;
    my ($name, $chr, $strand, $exonCount, $starts, $ends) = @fields[0..2, 7..9];
    my @starts = split ',', $starts;
    my @ends = split ',', $ends;
    my ($intronStarts, $intronEnds) = &gpeParser::getIntrons(\@starts, \@ends);
    my $supportedIntronCount = 0;
    for(my $i = 0; $i < @$intronStarts; $i++){
        $supportedIntronCount++ if exists $juncHash{"$chr:$strand:$intronStarts->[$i]-$intronEnds->[$i]"};
    }
    say join "\t", ($name, $exonCount - 1, $supportedIntronCount);
    $summaryHash{$exonCount-1}{transCount}++;
    $summaryHash{$exonCount-1}{supportedTrans}{$supportedIntronCount}++;
}

for my $juncCount(sort{$a<=>$b} keys %summaryHash){
    my $juncCountH = $summaryHash{$juncCount};
    my $supportedTransH = $juncCountH->{supportedTrans};
    for my $supportedIntronCount(sort{$a<=>$b} keys %$supportedTransH){
        say STDERR join "\t", ($juncCount, $supportedIntronCount, $juncCountH->{transCount},
                                $supportedTransH->{$supportedIntronCount});
    }
}



