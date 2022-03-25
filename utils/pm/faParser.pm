#!/bin/env perl
package faParser;

use strict;
use 5.010;
require Exporter;

1;

sub isPolyATailingSignal(){
    my ($seq) = @_;
    return 0 if length $seq < 6;
    my %signal = ('AATAAA' => '',
                  'AAATAA' => '',
                  'ATAAAA' => '',
                  'ATTAAA' => '',
                  'ATAAAT' => '',
                  'ATAAAG' => '',
                  'CAATAA' => '',
                  'TAATAA' => '',
                  'ATAAAC' => '',
                  'AAAATA' => '',
                  #'AAAAAA' => '',
                  'AAAAAG' => ''
                  );
    if(exists $signal{$seq}){
        return $seq;
    }else{
        return '';
    }
}

sub endWithPolyATailingSignal(){
    my ($seq) = @_;
    return 0 if length $seq < 6;
    my $tail = substr $seq, length ($seq) - 6;
    return &isPolyATailingSignal($tail);
}

sub reverseComplement(){
    my ($seq) = @_;
    $seq = join '', reverse (split '', $seq);
    $seq =~ tr/ATCG/TAGC/;
    $seq =~ tr/atcg/tagc/;
    $seq;
}