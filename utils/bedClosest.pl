#!/usr/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::bedParser;

GetOptions(
            'h|help' => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
my $BED;
open $BED, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";
my %bedHash = &bedParser::buildHash2($BED);
for my $chr(keys %bedHash){
    my $chrV = $bedHash{$chr};
    for my $strand(keys %$chrV){
        my @lines = @{$chrV->{$strand}};
        next if @lines < 2;
        my ($start1, $end1, $line1) = @{$lines[0]};
        my ($start2, $end2, $line2) = @{$lines[1]};
        if($strand eq '+'){
            my $len1To2 = $end2 - $end1;
            say join "\t", ($line1, $line2, $len1To2);
            for(my $i = 2; $i <= $#lines; $i++){
                my ($start3, $end3, $line3) = @{$lines[$i]};
                my $len2To3 = $end3 - $end2;
                if($len1To2 <= $len2To3){
                    say join "\t", ($line2, $line1, $len1To2);
                }else{
                    say join "\t", ($line2, $line3, $len2To3);
                }
                $len1To2 = $len2To3;
                $end2 = $end3;
                $line1 = $line2;
                $line2 = $line3;
            }
            say join "\t", ($line2, $line1, $len1To2);
        }else{
            my $len1To2 = $start2 - $start1;
            say join "\t", ($line1, $line2, $len1To2);
            for(my $i = 2; $i <= $#lines; $i++){
                my ($start3, $end3, $line3) = @{$lines[$i]};
                my $len2To3 = $start3 - $start2;
                if($len1To2 <= $len2To3){
                    say join "\t", ($line2, $line1, $len1To2);
                }else{
                    say join "\t", ($line2, $line3, $len2To3);
                }
                $len1To2 = $len2To3;
                $start2 = $start3;
                $line1 = $line2;
                $line2 = $line3;
            }
            say join "\t", ($line2, $line1, $len1To2);
        }
    }
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.bed >OUTPUT.tsv
    If INPUT.bed isn't specified, input from STDIN
Option:
    -h --help       Print this help information
HELP
    exit(-1);
}