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
use pm::gpeParser;
use pm::bedParser;

my ($gpeFile, $bin);
GetOptions(
            'g|gpe=s'   => \$gpeFile,
            'b|bin'     => \$bin,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
my $gpeFH;
open $gpeFH, "$gpeFile" or die "Can't read file ($gpeFile): $!";
my %gpeHash = &gpeParser::buildHash2($gpeFH, $bin);

while(<IN>){
    chomp;
    my ($chr, $start, $end, $strand, $sizes, $relStarts) = (split "\t")[0, 1, 2, 5, 10, 11];
    my @sizes = split ',', $sizes;
    my @relStarts = split ',', $relStarts;
    my ($starts, $ends) = &bedParser::getAbsLoc($start, \@sizes, \@relStarts);
    my (%eGeneHash, %iGeneHash);
TRANS:
    for my $trans (@{$gpeHash{$chr}{$strand}}){
        my ($transStarts, $transEnds, $geneName) = (split "\t", $trans->[2])[8, 9, 11];
        if($start < $trans->[1] && $end > $trans->[0]){
            $iGeneHash{$geneName} = '';
            my @transStarts = split ',', $transStarts;
            my @transEnds = split ',', $transEnds;
            for(my $i = 0; $i < $#transStarts; $i++){
                if($transEnds[$i] == $ends->[0]){
                    $eGeneHash{$geneName} = '';
                    next TRANS;
                }
            }
            for(my $i = 1; $i < @transStarts; $i++){
                if($transStarts[$i] == $starts->[-1]){
                    $eGeneHash{$geneName} = '';
                    next TRANS;
                }
            }
            for(my $i = 0; $i < @transStarts; $i++){
                last if $transStarts[$i] > $ends->[-1];
                if($starts->[0] < $transEnds[$i] && $ends->[0] > $transStarts[$i]){
                    for(my $j = $i; $j < @transStarts; $j++){
                        next TRANS if $transStarts[$j] > $ends->[-1];
                        if($starts->[-1] < $transEnds[$j] && $ends->[-1] > $transStarts[$j]){
                            $eGeneHash{$geneName} = '';
                            next TRANS;
                        }
                    }
                }
            }
        }
        last if $trans->[0] > $end;
    }
    my @eGenesName = keys %eGeneHash;
    my @iGenesName = keys %iGeneHash;
    my ($type, $genesName);
    if(@eGenesName != 0){
        $type = 'E';
        $genesName = join ',', @eGenesName;
    }elsif(@iGenesName != 0){
        $type = 'I';
        $genesName = join ',', @iGenesName;
    }else{
        $type = 'IG';
        $genesName = 'NA';
    }
    say join "\t", ($_, $type, $genesName);
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Option:
    -g --gpe    FILE    The gpe file
    -b --bin            With bin in gpe
    -h --help           Print this help information
HELP
    exit(-1);
}