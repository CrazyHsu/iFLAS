#!/usr/bin/env perl

use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::gpeParser;
use pm::bedParser;

my ($bin, $tss, $tts, $intron, $exon, $cds, $utr, $prime, $cmpl, $upstream, $dnstream, $chrSizeFile, $chain, $addIndex);
sub usage{
    my $scriptName = basename $0;
    print <<HELP;
Usage: perl $scriptName OPTION INPUT.gpe >OUTPUT.bed
    if INPUT.gpe isn't specified, input from STDIN
Example: perl $scriptName -b -g /data/genome/chr.size/hg19.size --upstream 1000 hg19.refGene.gpe >hg19.refGene.bed
Option:
    -b --bin                With bin column
    -i --intron             Fetch introns in each transcript
    -e --exon               Fetch exons in each transcript
    -c --cds                Fetch CDS in each transcript
    -u --utr                Fetch UTRs in each transcript, 5'UTR then 3'UTR (or 3' first)
    -p --prime      INT     5 for 5'UTR, 3 for 3'UTR(force -u)
       --complete           Only fetch UTR for completed transcripts
       --upstream   INT     Fetch upstream INT intergenic regions(force -g)
       --downstream INT     Fetch downstream INT intergenice regions(force -g)
    -g --chrSize    FILE    Tab-separated file with two columns: chr name and its length
    -s --single             Bundle all features into single line for each transcript
       --addIndex           Add exon/intron/CDS/UTR index as suffix of name in the 4th column
    -h --help               Print this help information
HELP
    exit(-1);
}

GetOptions(
            'b|bin'         => \$bin,
            'tss'           => \$tss,
            'tts'           => \$tts,
            'i|intron'      => \$intron,
            'e|exon'        => \$exon,
            'c|cds'         => \$cds,
            'u|utr'         => \$utr,
            'p|prime=i'     => \$prime,
            'complete'      => \$cmpl,
            'upstream=i'    => \$upstream,
            'downstream=i'  => \$dnstream,
            'g|chrSize=s'   => \$chrSizeFile,
            's|single'      => \$chain,
            'addIndex'      => \$addIndex,
            'h|help'        => sub{&usage()}
        )||usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open GPE, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";

my %chrSize;
if(defined $chrSizeFile){
    open chrSize, "$chrSizeFile" or die "Can't open $chrSizeFile: $!";
    while(<chrSize>){
        chomp;
        my ($chr, $length) = split "\t";
        $chrSize{$chr} = $length;
    }
}

while(<GPE>){
    chomp;
    my @fields = split "\t";
    shift @fields if defined $bin;
    my ($name, $chr, $strand, $start, $end, $cdsStart, $cdsEnd, $exonStarts, $exonEnds, $gene, $cdsStartStat, $cdsEndStat) = @fields[0..6, 8, 9, 11..13];
    my @exonStarts = split ',', $exonStarts;
    my @exonEnds = split ',', $exonEnds;
    if(defined $tss){
        my ($tssStart, $tssEnd);
        if($strand eq '+'){
            $tssStart = $start;
            $tssEnd = $tssStart + 1;
        }else{
            $tssEnd = $end;
            $tssStart = $tssEnd - 1;
        }
        say join "\t", ($chr, $tssStart, $tssEnd, $name, 0, $strand);
    }elsif(defined $tts){
        my ($ttsStart, $ttsEnd);
        if($strand eq '-'){
            $ttsStart = $start;
            $ttsEnd = $ttsStart + 1;
        }else{
            $ttsEnd = $end;
            $ttsStart = $ttsEnd - 1;
        }
        say join "\t", ($chr, $ttsStart, $ttsEnd, $name, 0, $strand);
    }elsif(defined $intron){
        my ($intronStarts, $intronEnds) = &gpeParser::getIntrons(\@exonStarts, \@exonEnds);
        next if @$intronStarts == 0;
        if(defined $chain){
            say join "\t", ($chr, $intronStarts->[0],
                            $intronEnds->[-1],
                            $name,
                            0,
                            $strand,
                            $intronStarts->[0],
                            $intronEnds->[-1],
                            '0,0,0',
                            scalar @$intronStarts,
                            (join ',', &gpeParser::getSizes($intronStarts, $intronEnds)),
                            (join ',', &bedParser::getRelStarts($intronStarts))
                            );
        }else{
            for (my $i = 0; $i < @$intronStarts; $i++){
                my ($intronStart, $intronEnd) = ($intronStarts->[$i], $intronEnds->[$i]);
                my $intronID = $name;
                $intronID .= '.' . ($strand eq '+' ? $i+1 : @$intronStarts - $i) if defined $addIndex;
                say join "\t", ($chr,
                                $intronStart,
                                $intronEnd,
                                $intronID,
                                0,
                                $strand,
                                $name,
                                $intronEnd - $intronStart);
            }
        }
    }elsif(defined $exon){
        if(defined $chain){
            say join "\t", ($chr, $exonStarts[0],
                                $exonEnds[-1],
                                $name,
                                0,
                                $strand,
                                $exonStarts[0],
                                $exonEnds[-1],
                                '0,0,0',
                                $#exonStarts + 1,
                                (join ',', &gpeParser::getSizes(\@exonStarts, \@exonEnds)),
                                (join ',', &bedParser::getRelStarts(\@exonStarts))
                                );
        }else{
            for (my $i = 0; $i < @exonStarts; $i++){
                my $exonID = $name;
                $exonID .= '.' . ($strand eq '+' ? $i+1 : @exonStarts - $i) if defined $addIndex;
                say join "\t", ($chr,
                                $exonStarts[$i],
                                $exonEnds[$i],
                                $exonID,
                                0,
                                $strand,
                                $name,
                                $exonEnds[$i] - $exonStarts[$i],
                                $gene
                                );
            }
        }
    }elsif(defined $cds){
        my ($cdsStarts, $cdsEnds) = &gpeParser::getCDSExons($cdsStart, $cdsEnd, \@exonStarts, \@exonEnds);
        next if @$cdsStarts == 0;
        if(defined $chain){
            say join "\t", ($chr, $cdsStarts->[0],
                                $cdsEnds->[-1],
                                $name,
                                0,
                                $strand,
                                $cdsStarts->[0],
                                $cdsEnds->[-1],
                                '0,0,0',
                                scalar @$cdsStarts,
                                (join ',', &gpeParser::getSizes($cdsStarts, $cdsEnds)),
                                (join ',', &bedParser::getRelStarts($cdsStarts)),
                                $gene);
        }else{
            for (my $i = 0; $i < @$cdsStarts; $i++){
                my $cdsID = $name;
                $cdsID .= '.' . ($strand eq '+' ? $i+1 : @$cdsStarts - $i) if defined $addIndex;
                say join "\t", ($chr,
                                $cdsStarts->[$i],
                                $cdsEnds->[$i],
                                $cdsID,
                                0,
                                $strand,
                                $name,
                                $cdsEnds->[$i] - $cdsStarts->[$i]
                                );
            }
        }
    }elsif(defined $utr){
        my ($utrStarts, $utrEnds);
        next if defined $cmpl && ($cdsStartStat !~ /cmpl/ || $cdsEndStat !~ /cmpl/);
        if(defined $prime){
            if($cdsStart != $start && ($strand eq '+' && $prime == 5 || $strand eq '-' && $prime == 3)){
                ($utrStarts, $utrEnds) = &gpeParser::getLeftUTRExons($cdsStart, \@exonStarts, \@exonEnds);
            }elsif($cdsEnd != $end && ($strand eq '+' && $prime == 3 || $strand eq '-' && $prime == 5)){
                ($utrStarts, $utrEnds) = &gpeParser::getRightUTRExons($cdsEnd, \@exonStarts, \@exonEnds);
            }else{
                next;
            }
            my @blockSizes = &gpeParser::getSizes($utrStarts, $utrEnds);
            my $blockRelStarts = join ',', &bedParser::getRelStarts($utrStarts);
            if(defined $chain){
                say join "\t", ($chr, $utrStarts->[0], $utrEnds->[-1], $name, 0, $strand,
                                $utrEnds->[-1], $utrEnds->[-1], '0,0,0', $#blockSizes + 1,
                                (join ',', @blockSizes),
                                $blockRelStarts);
            }else{
                for (my $i = 0; $i < @$utrStarts; $i++){
                    my $utrID = $name;
                    $utrID .= '.' . ($strand eq '+' ? $i+1 : @$utrStarts - $i) if defined $addIndex;
                    say join "\t", ($chr,
                                    $utrStarts->[$i],
                                    $utrEnds->[$i],
                                    $utrID,
                                    0,
                                    $strand,
                                    $name,
                                    $blockSizes[$i]
                                    );
                }
            }
        }else{
            if($cdsStart != $start){
                ($utrStarts, $utrEnds) = &gpeParser::getLeftUTRExons($cdsStart, \@exonStarts, \@exonEnds);
                my @blockSizes = &gpeParser::getSizes($utrStarts, $utrEnds);
                my $blockRelStarts = join ',', &bedParser::getRelStarts($utrStarts);
                if(defined $chain){
                    say join "\t", ($chr, $utrStarts->[0], $utrEnds->[-1], $name, 0, $strand,
                                    $utrEnds->[-1], $utrEnds->[-1], '0,0,0', $#blockSizes + 1,
                                    (join ',', @blockSizes),
                                    $blockRelStarts);
                }else{
                    for (my $i = 0; $i < @$utrStarts; $i++){
                        my $utrID = $name;
                        $utrID .= '.L.' . ($strand eq '+' ? $i+1 : @$utrStarts - $i) if defined $addIndex;
                        say join "\t", ($chr,
                                        $utrStarts->[$i],
                                        $utrEnds->[$i],
                                        $utrID,
                                        0,
                                        $strand,
                                        $name,
                                        $blockSizes[$i]
                                        );
                    }
                }
            }
            
            if($cdsEnd != $end){
                ($utrStarts, $utrEnds) = &gpeParser::getRightUTRExons($cdsEnd, \@exonStarts, \@exonEnds);
                my @blockSizes = &gpeParser::getSizes($utrStarts, $utrEnds);
                my $blockRelStarts = join ',', &bedParser::getRelStarts($utrStarts);
                if(defined $chain){
                    say join "\t", ($chr, $utrStarts->[0], $utrEnds->[-1], $name, 0, $strand,
                                    $utrEnds->[-1], $utrEnds->[-1], '0,0,0', $#blockSizes + 1,
                                    (join ',', @blockSizes),
                                    $blockRelStarts);
                }else{
                    for (my $i = 0; $i < @$utrStarts; $i++){
                        my $utrID = $name;
                        $utrID .= '.R.' . ($strand eq '+' ? $i+1 : @$utrStarts - $i) if defined $addIndex;
                        say join "\t", ($chr,
                                        $utrStarts->[$i],
                                        $utrEnds->[$i],
                                        $utrID,
                                        0,
                                        $strand,
                                        $name,
                                        $blockSizes[$i]
                                        );
                    }
                }
            }
        }
    }elsif(defined $upstream || defined $dnstream){
        die "Please specify --chrSize to avoid flank extended out of chromosome" unless $chrSizeFile;
        my ($intergenicStart, $intergenicEnd, $intergenicName);
        if(defined $upstream){
            if($strand eq '+'){
                $intergenicStart = $start - $upstream;
                $intergenicEnd = $intergenicStart + $upstream;
            }else{
                $intergenicStart = $end;
                $intergenicEnd = $intergenicStart + $upstream; 
            }
            $intergenicStart = 0 if $intergenicStart < 0;
            $intergenicEnd = $chrSize{$chr} if $intergenicEnd > $chrSize{$chr};
            $intergenicName = "$name.up". ($intergenicEnd - $intergenicStart) . "bp";
            say join "\t", ($chr, $intergenicStart, $intergenicEnd, $intergenicName, 0, $strand);
        }
        if(defined $dnstream){
            if($strand eq '+'){
                $intergenicStart = $end;
                $intergenicEnd = $intergenicStart + $dnstream;
            }else{
                $intergenicStart = $start - $dnstream;
                $intergenicEnd = $intergenicStart + $dnstream; 
            }
            $intergenicStart = 0 if $intergenicStart < 0;
            $intergenicEnd = $chrSize{$chr} if $intergenicEnd > $chrSize{$chr};
            $intergenicName = "$name.dn" . ($intergenicEnd - $intergenicStart) . "bp";
            say join "\t", ($chr, $intergenicStart, $intergenicEnd, $intergenicName, 0, $strand);
        }
    }
}