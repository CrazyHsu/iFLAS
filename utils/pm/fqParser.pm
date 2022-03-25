#!/bin/env perl
package fqParser;

use strict;
use 5.010;
require Exporter;

1;

####    Argument        Type    Description
#       seq             string  read sequence
#       base            string  base to count

####    Return          Type    Description
#       count           int     occurrence of base in seq
sub countBase{
    my ($seq, $base) = @_;
    my @bases = $seq =~ /([$base])/ig;
    return scalar @bases;
}

####    Argument        Type    Description
#       seq             string  read sequence
#       base            string  base to count

####    Return          Type    Description
#       fraction        double  fraction of base in seq
sub getBaseFraction{
    my ($seq, $base) = @_;
    return &countBase($seq, $base)/ length $seq;
}

sub getQualOffset{
    my ($platform) = @_;
    if($platform =~/^(S)|(Sanger)|(L)|(Illumina1.8+)$/i){
        return 33;
    }elsif($platform =~/^(X)|(Solexa)|(I)|(Illumina1.3+)|(J)|(Illumina1.5+)$/i){
        return 64;
    }else{
        die "Please specify the correct scoring platform
        (3 for S|Sanger or L|Illumina1.8+; 64 for X|Solexa or I|Illumina1.3+ or J|Illumina1.5+";
    }
}

####    Argument        Type    Description
#       qualSeq         string  quality string
#       platform        string  platform of quality scoring (   33 for S|Sanger or L|Illumina1.8+;
#                                                               64 for X|Solexa or I|Illumina1.3+ or J|Illumina1.5+
#                                                           )
#       qualCutoff      int     Low quality threshold

####    Return          Type    Description
#       count           int     occurrence of low quality bases in seq
sub countLowQualBase{
    my ($qualSeq, $platform, $qualCutoff) = @_;
    my $count = 0;
    my $qualOffset = &getQualOffset($platform);
    for (my $i=0; $i < length $qualSeq; $i++){
        $count++ if ord(substr $qualSeq, $i, 1) - $qualOffset < $qualCutoff;
    }
    $count;
}

####    Argument        Type    Description
#       qualSeq         string  quality string
#       platform        string  platform of quality scoring (   33 for S|Sanger or L|Illumina1.8+;
#                                                               64 for X|Solexa or I|Illumina1.3+ or J|Illumina1.5+
#                                                           )
#       qualCutoff      int     Low quality threshold

####    Return          Type    Description
#       fraction        double  fraction of low quality bases in seq
sub getLowQualBaseFraction{
    my ($qualSeq, $platform, $qualCutoff) = @_;
    &countLowQualBase(@_) / length $qualSeq;
}

####    Argument        Type    Description
#       qualSeq         string  quality string
#       platform        string  platform of quality scoring (   33 for S|Sanger or L|Illumina1.8+;
#                                                               64 for X|Solexa or I|Illumina1.3+ or J|Illumina1.5+
#                                                           )

####    Return          Type    Description
#       qualSum         int     quality sum
sub countQualSum{
    my ($qualSeq, $platform) = @_;
    my $qualSum = 0;
    my $qualOffset = &getQualOffset($platform);
    for (my $i=0; $i < length $qualSeq; $i++){
        $qualSum += ord(substr $qualSeq, $i, 1) - $qualOffset;
    }
    $qualSum;
}

####    Argument        Type    Description
#       qualSeq         string  quality string
#       platform        string  platform of quality scoring (   33 for S|Sanger or L|Illumina1.8+;
#                                                               64 for X|Solexa or I|Illumina1.3+ or J|Illumina1.5+
#                                                           )

####    Return          Type    Description
#       aveQual         int     average quality
sub countAveQual{
    my ($qualSeq, $platform) = @_;
    &countQualSum($qualSeq, $platform) / length $qualSeq;
}

####    Argument        Type    Description
#       tail            string  in the form of 'INT1,INT2', where INT1 is the tail length and INT2 is the cutoff quality scroe
#       qualSeq         string  quality string
#       platform        string  platform of quality scoring (   33 for S|Sanger or L|Illumina1.8+;
#                                                               64 for X|Solexa or I|Illumina1.3+ or J|Illumina1.5+
#                                                           )

####    Return          Type    Description
#       isBadTail       int     1 for bad tail, 0 for not bad tail
sub isBadTail{
    my ($tail, $qualSeq, $platform) = @_;
    my ($tailL, $tailV) = split /\D/, $tail;
    my @tailQual = split "", (substr $qualSeq, length ($qualSeq) - $tailL);
    my $qualOffset = &getQualOffset($platform);
    for my $qual (@tailQual){
        if (ord($qual) - $qualOffset < $tailV){
            return 1;
         }
    }
    0;
}

####    Argument    Type    Description
#       seq         STR     sequence line
#       qual        STR     quality line
#       trimLen     STR     length to trim
#       side                5' or 3' end to trim

####    Return      Type    Description
#       (seq, qual) ARRAY   sequence and quality line after trimming
sub trim{
    my ($seq, $qual, $trimLen, $side) = @_;
    my $len = length($seq);
    $trimLen = $len if $trimLen > $len;
    my $lenAfTrim = $len - $trimLen;
    if($side == 5){
        $seq = substr( $seq, $trimLen);
        $qual = substr( $qual, $trimLen);
    }elsif($side ==3 ){        
        $seq = substr( $seq, 0, $lenAfTrim );
        $qual = substr( $qual, 0, $lenAfTrim);
    }
    ($seq, $qual);
}

sub trimByRegExp{
    my ($seq, $qual, $regExp, $side) = @_;
    my $oriLen = length($seq);
    if($side == 5){
        $seq =~ s/^$regExp//;
        my $lenTrimmed = $oriLen - length $seq;
        $qual = substr($qual, $lenTrimmed)
    }elsif($side == 3){
        $seq =~ s/$regExp$//;
        $qual = substr($qual, 0, length $seq)
    }
    ($seq, $qual);
}

sub toQualValue(){
    my ($qual, $platform) = @_;
    my $offset = getQualOffset($platform);
    my @BQs;
    $qual =~ s/(.)/push @BQs, ord($1)-$offset/eg;
    return @BQs;
}


