#!/usr/bin/perl -w

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;

my $extension;
GetOptions(
            'e|extension=i' => \$extension,
            'h|help'        => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";

my (%geneHash, %novelHash);
while(<IN>){
    chomp;
    my ($chr, $start, $end, $name, $score, $strand, $blockCount, $blockSizes, $relStarts, $type, $genesName) = (split "\t")[0..5, 9..13];
    my @relStarts = split ",", $relStarts;
    my @readSizes = split ",", $blockSizes;
    my ($readStarts, $readEnds) = &common::getAbsLoc($start, \@readSizes, \@relStarts);
    my $exonChain = join ';', map{"$readStarts->[$_]-$readEnds->[$_]"}0..(@$readStarts-1);
    if($type eq 'E'){
        for my $geneName(split ',', $genesName){
            $geneHash{$geneName}{read}{$name} = [$exonChain, $score];
            $geneHash{$geneName}{chr} = $chr;
            $geneHash{$geneName}{strand} = $strand;
        }
    }else{
        push @{$novelHash{$chr}{$strand}}, [$start, $end, $exonChain, $name, $score];
    }
}

for my $geneName(keys %geneHash){
    my $geneNameV = $geneHash{$geneName}{read};
    my @readNames = keys %$geneNameV;
    next if @readNames < 2;
    my $strand = $geneHash{$geneName}{strand};
    my %altHash;
    my $indexStart = $strand eq '+' ? 0 : 1;
    for(my $r1 = 0; $r1 < $#readNames; $r1++){
        my ($exons1, $score1) = @{$geneHash{$geneName}{read}{$readNames[$r1]}};
        my @exons1 = split ';', $exons1;
        my $indexEnd1 = $strand eq '+' ? $#exons1-1 : $#exons1;
        for(my $r2 = $r1 + 1; $r2 <= $#readNames; $r2++){
            my ($exons2, $score2) = @{$geneHash{$geneName}{read}{$readNames[$r2]}};
            my @exons2 = split ';', $exons2;
            my $indexEnd2 = $strand eq '+' ? $#exons2 - 1 : $#exons2;
            for(my $i = $indexStart; $i <= $indexEnd1; $i++){
                my ($exon1Start, $exon1End) = split '-', $exons1[$i];
                for(my $j = $indexStart; $j <= $indexEnd2; $j++){
                    my ($exon2Start, $exon2End) = split '-', $exons2[$j];
                    if($exon1Start < $exon2End + $extension && $exon1End + $extension > $exon2Start){
                        if($strand eq '+'){
                            if($exon1End < $exon2End && $exon2End < (split '-', $exons1[$i+1])[0]){
                                $altHash{"$exon1End-$exon2End"}{exc}{$readNames[$r1]} = $score1;
                                $altHash{"$exon1End-$exon2End"}{inc}{$readNames[$r2]} = $score2;
                            }elsif($exon1End > $exon2End && $exon1End < (split '-', $exons2[$j+1])[0]){
                                $altHash{"$exon2End-$exon1End"}{inc}{$readNames[$r1]} = $score1;
                                $altHash{"$exon2End-$exon1End"}{exc}{$readNames[$r2]} = $score2;
                            }
                        }else{
                            if($exon1Start < $exon2Start && $exon1Start > (split '-', $exons2[$j-1])[1]){
                                $altHash{"$exon1Start-$exon2Start"}{inc}{$readNames[$r1]} = $score1;
                                $altHash{"$exon1Start-$exon2Start"}{exc}{$readNames[$r2]} = $score2;
                            }elsif($exon1Start > $exon2Start && $exon2Start > (split '-', $exons1[$i-1])[1]){
                                $altHash{"$exon2Start-$exon1Start"}{exc}{$readNames[$r1]} = $score1;
                                $altHash{"$exon2Start-$exon1Start"}{inc}{$readNames[$r2]} = $score2;
                            }
                        }
                    }
                    last if $exon2Start > $exon1End;
                }
            }
        }
    }
    for my $alt(keys %altHash){
        my ($altStart, $altEnd) = split '-', $alt;
        my @incReads = keys %{$altHash{$alt}{inc}};
        my @excReads = keys %{$altHash{$alt}{exc}};
        my ($incReads, $excReads) = (0, 0);
        for my $read(@incReads){
            $incReads += $altHash{$alt}{inc}{$read};
        }
        for my $read(@excReads){
            $excReads += $altHash{$alt}{exc}{$read};
        }
        say join "\t", ($geneHash{$geneName}{chr}, $altStart, $altEnd, "$geneName:$alt", int($incReads/($incReads+$excReads)*1000), $strand, $geneName,
                        (join ',', @incReads), $incReads, (join ',', @excReads), $excReads);
    }
}

my $inc = 1;
while(my ($chr, $chrV) = each %novelHash){
    for my $strand (keys %$chrV){
        my @sortedStrandV = sort {$a->[0]<=>$b->[0]}@{$chrV->{$strand}};
        my @reads = ($sortedStrandV[0]);
        my $clusterEnd = $sortedStrandV[0]->[1];
        for(my $i = 1; $i <= $#sortedStrandV; $i++){
            my ($start, $end) = @{$sortedStrandV[$i]};
            if($start < $clusterEnd){
                push @reads, $sortedStrandV[$i];
                $clusterEnd = $end if $end > $clusterEnd;
            }else{
                &findA5SS($chr, $strand, @reads);
                $inc++;
                @reads = ($sortedStrandV[$i]);
                $clusterEnd = $end;
            }
        }
        &findA5SS($chr, $strand, @reads);
    }
}

sub findA5SS(){
    my ($chr, $strand, @reads) = @_;
    return if @reads < 2;
    my %altHash;
    my $indexStart = $strand eq '+' ? 0 : 1;
    for(my $r1 = 0; $r1 < $#reads; $r1++){
        my @exons1 = split ';', $reads[$r1]->[2];
        my ($readName1, $score1) = @{$reads[$r1]}[3, 4];
        my $indexEnd1 = $strand eq '+' ? $#exons1-1 : $#exons1;
        for(my $r2 = $r1 + 1; $r2 <= $#reads; $r2++){
            my @exons2 = split ';', $reads[$r2]->[2];
            my ($readName2, $score2) = @{$reads[$r2]}[3, 4];
            my $indexEnd2 = $strand eq '+' ? $#exons2-1 : $#exons2;
            for(my $i = $indexStart; $i <= $indexEnd1; $i++){
                my ($exon1Start, $exon1End) = split '-', $exons1[$i];
                for(my $j = $indexStart; $j <= $indexEnd2; $j++){
                    my ($exon2Start, $exon2End) = split '-', $exons2[$j];
                    if($exon1Start < $exon2End + $extension && $exon1End + $extension > $exon2Start){
                        if($strand eq '+'){
                            if($exon1End < $exon2End && $exon2End < (split '-', $exons1[$i+1])[0]){
                                $altHash{"$exon1End-$exon2End"}{exc}{$readName1} = $score1;
                                $altHash{"$exon1End-$exon2End"}{inc}{$readName2} = $score2;
                            }elsif($exon1End > $exon2End && $exon1End < (split '-', $exons2[$j+1])[0]){
                                $altHash{"$exon2End-$exon1End"}{inc}{$readName1} = $score1;
                                $altHash{"$exon2End-$exon1End"}{exc}{$readName2} = $score2;
                            }
                        }else{
                            if($exon1Start < $exon2Start && $exon1Start > (split '-', $exons2[$j-1])[1]){
                                $altHash{"$exon1Start-$exon2Start"}{inc}{$readName1} = $score1;
                                $altHash{"$exon1Start-$exon2Start"}{exc}{$readName2} = $score2;
                            }elsif($exon1Start > $exon2Start && $exon2Start > (split '-', $exons1[$i-1])[1]){
                                $altHash{"$exon2Start-$exon1Start"}{exc}{$readName1} = $score1;
                                $altHash{"$exon2Start-$exon1Start"}{inc}{$readName2} = $score2;
                            }
                        }
                    }
                    last if $exon2Start > $exon1End;
                }
            }
        }
    }
    for my $alt(keys %altHash){
        my ($altStart, $altEnd) = split '-', $alt;
        my @incReads = keys %{$altHash{$alt}{inc}};
        my @excReads = keys %{$altHash{$alt}{exc}};
        my ($incReads, $excReads) = (0, 0);
        for my $read(@incReads){
            $incReads += $altHash{$alt}{inc}{$read};
        }
        for my $read(@excReads){
            $excReads += $altHash{$alt}{exc}{$read};
        }
        say STDERR join "\t", ($chr, $altStart, $altEnd, "Novel$inc:$alt", int($incReads/($incReads+$excReads)*1000), $strand, '.',
                        (join ',', @incReads), $incReads, (join ',', @excReads), $excReads);
    }
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT.bed12+
    If INPUT isn't specified, input from STDIN
    OUTPUT.bed6+ contains 5 additional columns:
            gene (. for novel locus),
            inclusion reads, read count of previous column,
            exclusion reads, read count of previous column
Options:
    -e --extension  INT Extend junction block
    -h --help           Print this help information
HELP
    exit(-1);
}