#!/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Options:
    -j --junc   FILE    Junction file
    -t --type   INT     5 for A5SS, 3 for A3SS
    -h --help           Print this help information
HELP
    exit(-1);
}

my ($juncFile, $type);
GetOptions(
            'j|junc=s'  => \$juncFile,
            'type=i'    => \$type,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
open JUNC, "$juncFile" or die "Can't read file ($juncFile): $!";

my %juncHash;
while(<JUNC>){
    my ($chr, $start, $name, $score, $strand, $sizes, $relStarts) = (split "\t")[0, 1, 3..5, 10, 11];
    my @sizes = split ',', $sizes;
    my @relStarts = split ',', $relStarts;
    my ($blockStarts, $blockEnds) = &common::getAbsLoc($start, \@sizes, \@relStarts);
    my @sites;
    if($type == 5){
        @sites = $strand eq '+' ? @$blockEnds[0..@$blockStarts-2] : @$blockStarts[1..@$blockStarts-1];
    }elsif($type == 3){
        @sites = $strand eq '-' ? @$blockEnds[0..@$blockStarts-2] : @$blockStarts[1..@$blockStarts-1];
    }
    for my $site(@sites){
        push @{$juncHash{"$chr:$strand"}{$site}}, [$name, $score];
    }
}

while(<IN>){
    chomp;
    my ($chr, $start, $end, $name, $strand) = (split "\t")[0..3, 5];
    next if !exists $juncHash{"$chr:$strand"}{$start} || !exists $juncHash{"$chr:$strand"}{$end};
    $_ =~ s/Novel/NGSNovel/;
    say STDERR;
    my ($incCount, $excCount, @incJuncs, @excJuncs, @gene);
    if($type == 5){
        if($strand eq '+'){
            for my $junc(@{$juncHash{"$chr:$strand"}{$start}}){
                push @excJuncs, $junc->[0];
                $excCount += $junc->[1];
            }
            for my $junc(@{$juncHash{"$chr:$strand"}{$end}}){
                push @incJuncs, $junc->[0];
                $incCount += $junc->[1];
            }
        }else{
            for my $junc(@{$juncHash{"$chr:$strand"}{$start}}){
                push @incJuncs, $junc->[0];
                $incCount += $junc->[1];
            }
            for my $junc(@{$juncHash{"$chr:$strand"}{$end}}){
                push @excJuncs, $junc->[0];
                $excCount += $junc->[1];
            }
        }
    }elsif($type == 3){
        if($strand eq '-'){
            for my $junc(@{$juncHash{"$chr:$strand"}{$start}}){
                push @excJuncs, $junc->[0];
                $excCount += $junc->[1];
            }
            for my $junc(@{$juncHash{"$chr:$strand"}{$end}}){
                push @incJuncs, $junc->[0];
                $incCount += $junc->[1];
            }
        }else{
            for my $junc(@{$juncHash{"$chr:$strand"}{$start}}){
                push @incJuncs, $junc->[0];
                $incCount += $junc->[1];
            }
            for my $junc(@{$juncHash{"$chr:$strand"}{$end}}){
                push @excJuncs, $junc->[0];
                $excCount += $junc->[1];
            }
        }
    }
    my $genes = (split ':', $name)[0];
    if($genes =~ /^Novel/){
        $genes = '.';
        $name = "PB$name";
    }
    say join "\t", ($chr, $start, $end, $name, int($incCount/($incCount+$excCount)*1000), $strand,
                           $genes, (join ',', @incJuncs), $incCount, (join ',', @excJuncs), $excCount);
}

