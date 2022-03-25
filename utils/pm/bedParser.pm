#!/bin/env perl
package bedParser;

use strict;
use 5.010;
require Exporter;
use List::Util qw/sum/;

#our @ISA = qw/Exporter/;
#our @EXPORT = qw/getCDSLength/;
1;

sub getAbsStarts{ # old name (getGpeBlockEnds) is deprecated
    my ($start, $blockStarts) = @_;
    my @blockStarts = @$blockStarts;
    map{$start + $_}@blockStarts;
}

sub getAbsEnds{
    my ($start, $blockSizes, $blockStarts) = @_;
    my @blockStarts = @$blockStarts;
    my @blockSizes = @$blockSizes;
    map{$start + $blockStarts[$_] + $blockSizes[$_]}0..$#blockStarts;
}

sub getAbsLoc{
    my ($start, $blockSizes, $blockStarts) = @_;
    my @blockStarts = getAbsStarts($start, $blockStarts);
    my @blockSizes = @$blockSizes;
    my @blockEnds = map{$blockStarts[$_] + $blockSizes[$_]}0..$#blockStarts;
    return (\@blockStarts, \@blockEnds);
}

sub getRelStarts{
    my ($blockStarts) = @_;
    map{$_ - $blockStarts->[0]}@$blockStarts;
}

####    Argument    Type    Description
#       fileHandle  FH      A file handle of bed file

####    Return      Type    Description
#       hash        hash    A hash with chr, strand as key and with array containing [start, end, line] as values
sub buildHash{
    my ($fh) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        my ($chr, $start, $end, $strand) = @fields[0..2, 5];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, $_];
        }else{
            $hash{$chr}{$strand} = [[$start, $end, $_]];
        }
    }
    while(my ($chr, $chrV) = each %hash){
        for my $strand (keys %$chrV){
            my @sortedStrandV = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]}@{$chrV->{$strand}};
            $chrV->{$strand} = \@sortedStrandV;
        }
    }
    return %hash;
}

####    Argument    Type    Description
#       fileHandle  FH      A file handle of bed file

####    Return      Type    Description
#       hash        hash    A hash with chr, strand as key and with array containing [start, end, line] as values
sub buildHash2{
    my ($fh) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        my ($chr, $start, $end, $strand) = @fields[0..2, 5];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, $_];
        }else{
            $hash{$chr}{$strand} = [[$start, $end, $_]];
        }
    }
    while(my ($chr, $chrV) = each %hash){
        for my $strand (keys %$chrV){
            my @sortedStrandV;
            if($strand eq '+'){
                @sortedStrandV = sort {$a->[1]<=>$b->[1] || $a->[0]<=>$b->[0]}@{$chrV->{$strand}};
            }else{
                @sortedStrandV = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]}@{$chrV->{$strand}};
            }
            $chrV->{$strand} = \@sortedStrandV;
        }
    }
    return %hash;
}
