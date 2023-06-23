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
    -r          Ref PA Group File
    -a          Annotated PA
    -n          Novel PA
    -h  --help  Print this help information
HELP
    exit(-1);
}

my ($refPaFile, $annoPaFile, $novelPaFile);
GetOptions(
            'r=s'       => \$refPaFile,
            'a=s'       => \$annoPaFile,
            'n=s'       => \$novelPaFile,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
open REF, "$refPaFile" or die "Can't read file ($refPaFile): $!";
open ANNO, ">$annoPaFile" or die "Can't read file ($annoPaFile): $!";
open NOVEL, ">$novelPaFile" or die "Can't read file ($novelPaFile): $!";

my %refPA;
while(<REF>){
    chomp;
    my ($chr, $strand, $gene, $PAs) = (split "\t")[0..3];
    for my $PA(split ",", $PAs){
        push @{$refPA{"$chr:$strand:$gene"}}, $PA;
    }
}

while(<IN>){
    chomp;
    my ($chr, $strand, $gene, $PAs) = (split "\t")[0..3];
    next if $gene =~ /^chr/;
    my $isNovelAPA = 0;
PA: for my $PA(split ",", $PAs){
        for my $refPA(@{$refPA{"$chr:$strand:$gene"}}){
            if(abs($refPA - $PA) <= 30){
                say ANNO join "\t", ($chr, $strand, $gene, $PA);
                next PA;
            }
        }
        say NOVEL join "\t", ($chr, $strand, $gene, $PA);
        $isNovelAPA = 1;
    }
    if($isNovelAPA == 1){
        say STDOUT $_; # novel
    }else{
        say STDERR $_; # anno
    }
}

