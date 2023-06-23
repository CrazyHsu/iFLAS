#!/usr/bin/env perl

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

#my ($gpeFile, $bin);
GetOptions(
            #'g|gpe=s'   => \$gpeFile,
            #'b|bin'     => \$bin,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
#my $GPE;
#open $GPE, "$gpeFile" or die "Can't open $gpeFile: $!";
#
#my %gpeHash = &common::buildGpeHashByGeneNameWithTransAsArray($GPE, $bin);

my %gene2reads;
while(<IN>){
    chomp;
    my ($chr, $start, $end, $readName, $strand, $sizes, $relStarts, $gene) = (split "\t")[0..3, 5, 10..12];
    my @sizes = split ',', $sizes;
    my @relStarts = split ',', $relStarts;
    my $PA = $strand eq '+' ? $end : $start;
    my ($blockStarts, $blockEnds) = &common::getAbsLoc($start, \@sizes, \@relStarts);
    push @{$gene2reads{"$chr&$strand&$gene"}{$PA}}, [$blockStarts, $blockEnds, $readName];
}

for my $geneID(keys %gene2reads){
    my ($chr, $strand, $gene) = split '&', $geneID;
    my @PAs;
    for my $PA(sort{$a<=>$b}keys %{$gene2reads{$geneID}}){
        my @reads = @{$gene2reads{$geneID}{$PA}};
        my $constSite;
        if($strand eq '+'){
            $constSite = $reads[0]->[0]->[-1];
            for(my $i = 1; $i <= $#reads; $i++){
                my $nextReadLastExonUpSite = $reads[$i]->[0]->[-1];
                $constSite = $nextReadLastExonUpSite if $nextReadLastExonUpSite > $constSite;
            }
        }else{
            $constSite = $reads[0]->[1]->[0];
            for(my $i = 1; $i <= $#reads; $i++){
                my $nextReadLastExonUpSite = $reads[$i]->[1]->[0];
                $constSite = $nextReadLastExonUpSite if $nextReadLastExonUpSite < $constSite;
            }
        }
        push @PAs, [$PA, $constSite, scalar @reads];
        say STDERR join "\t", ($chr,
                               ($strand eq '+' ? ($constSite, $PA) : ($PA, $constSite)),
                               "$gene:$PA", scalar @reads, $strand);
    }
    my $PAs = join ',', map{$_->[0]}@PAs;
    my $constSites = join ',', map{$_->[1]}@PAs;
    my $supportReads = join ',', map{$_->[2]}@PAs;
    my $totoalReads = 0;
    $totoalReads += $_->[2] for(@PAs);
    my $freqs = join ',', map{$_->[2]/$totoalReads}(@PAs);
    say join "\t", ($chr, $strand, $gene, $PAs, $constSites, $supportReads, $freqs);
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Output:
    STDOUT: chr, strand, gene/locus ID, PAs and closest upstream sites
    STDERR: bed6 file composed of chr, paRegionStart, paRegionEnd, PA ID, supporting reads count and strand
Options:
    -h --help       Print this help information
HELP
    exit(-1);
}