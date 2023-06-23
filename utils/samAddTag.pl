#!/bin/env perl

use 5.012;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::samParser;

my ($coverage, $match, $misMatch, $nBase, $identity, $indel, $list, $listTag, $unmappedFile, $checkHardClip);

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.sam/bam >OUTPUT.sam 2>inconsistent.bam
    If INPUT isn't specified, input from STDIN
Option:
    -c --coverage           Add coverage as CV:f tag
    -m --match              Add number of matched bases as MC:i tag
    -x --mismatch           Add number of mismatch bases as XM:i tag
    -n --nBase              Add number of ambiguous bases (N in general) as XN:i tag
    -i --identity           Add identity as ID:f tag
       --indel              Add indels length as IN:i tag
       --list       TSV     Add tag from a tsv file (with read name and tag value as columns), matching with read name
       --listTag    STR     The prefix of tag adding by --list
    -u --unmapped   SAM     Output unmapped read to SAM file
       --checkHardClip      Whether to check consist of cigar length (consider hard clip) with sequence length
    -h --help               Print this help information
HELP
    exit(-1);
}

GetOptions(
            'c|coverage'    => \$coverage,
            'm|match'       => \$match,
            'x|mismatch'    => \$misMatch,
            'n|nBase'       => \$nBase,
            'i|identity'    => \$identity,
            'indel'         => \$indel,
            'list=s'        => \$list,
            'listTag=s'     => \$listTag,
            'u|unmapped=s'  => \$unmappedFile,
            'checkHardClip' => \$checkHardClip,
            'h|help'        => sub{usage()}
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

my %readName2TagValue;
if(defined $list){
    open LIST , "$list" or die "Can't open $list: $!";
    my ($listLine);
    while($listLine = <LIST>){
        chomp $listLine;
        my ($readName, $tagValue) = split "\t", $listLine;
        $readName2TagValue{$readName} = $tagValue;
    }
}

if(defined $unmappedFile){
    open UM, ">$unmappedFile" or die "Can't open $unmappedFile: $!";
}

while(<IN>){
    chomp;
    if (/^@/){
        say;
        next;
    }
    my @fields = split "\t";
    my ($readName, $flag, $chr, $cigar, $seq) = @fields[0, 1, 2, 5, 9];
    if(defined $unmappedFile && &samParser::isUnmapped($flag)){
        say UM $_;
        next;
    }
    my @addedTags;
    my $seqLength = length $seq;
    if(defined $checkHardClip && &samParser::getReadHardClippedLength($cigar) != $seqLength){
        say STDERR;
        next;
    }
    my ($alignedLength, $matchLength);
    my $md;
    $md  = &samParser::getTagValue($_, "MD") if (defined $match || defined $misMatch || defined $nBase || defined $identity);
    if(defined $coverage){
        my $readLength = &samParser::getReadLength($cigar);
        $alignedLength = &samParser::getAlignedLength($cigar);
        $coverage = sprintf "%.6f", $alignedLength / $readLength * 100;
        push @addedTags, "CV:f:$coverage";
    }
    if(defined $match){
        $matchLength = &samParser::getMatchLength($md);
        push @addedTags, "MC:i:$matchLength";
    }
    if(defined $misMatch){
        my $misMatchLength = &samParser::getMisMatchLength($md);
        push @addedTags, "XM:i:$misMatchLength";
    }
    if(defined $nBase){
        my $nBaseCount = &samParser::getNBaseCount($md);
        push @addedTags, "XN:i:$nBaseCount";
    }
    if(defined $identity){
        $matchLength = &samParser::getMatchLength($md) unless defined $matchLength;
        $alignedLength = &samParser::getAlignedLength($cigar) unless defined $alignedLength;
        my $identity = sprintf "%.6f", $matchLength / $alignedLength * 100;
        push @addedTags, "ID:f:$identity";
    }
    if(defined $indel){
        my $indelLength= &samParser::getIndelLength($cigar);
        push @addedTags, "IN:i:$indelLength";
    }
    if(defined $list){
        push @addedTags, "$listTag:$readName2TagValue{$readName}";
    }
    say join "\t", ($_, @addedTags);
}
