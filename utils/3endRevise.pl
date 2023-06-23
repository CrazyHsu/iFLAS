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

my $paFile;
GetOptions(
            'p|pa=s'    => \$paFile,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
open PA, $paFile or die "Can't read file ($paFile): $!";

my %read2PA;
while(<PA>){
    chomp;
    my ($pos, $reads) = (split "\t")[2, 3];
    for my $read(split ',', $reads){
        $read2PA{$read} = $pos;
    }
}

while(<IN>){
    chomp;
    my @fields = split "\t";
    
    my @blockSizes = split ',', $fields[10];
    my @relStarts = split ',', $fields[11];
    my ($blockStarts, $blockEnds) = &common::getAbsLoc($fields[1], \@blockSizes, \@relStarts);
    my $paPos = $read2PA{$fields[3]};
    if(defined $paPos){
        if($fields[5] eq '+'){
            $blockEnds->[-1] = $paPos if $paPos > $blockStarts->[-1];
        }else{
            $blockStarts->[0] = $paPos - 1 if $paPos <= $blockEnds->[0];
        }
    }
    @blockSizes = common::getSizes($blockStarts, $blockEnds);
    @relStarts = common::getRelStarts($blockStarts);
    say join "\t", ($fields[0], $blockStarts->[0], $blockEnds->[-1], @fields[3..5], $blockStarts->[0], $blockEnds->[-1],
                    @fields[8, 9], (join ',', @blockSizes), (join ',', @relStarts), @fields[12..$#fields]);
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Options:
    -p --pa     FILE    The PA file
    -h --help           Print this help information
HELP
    exit(-1);
}