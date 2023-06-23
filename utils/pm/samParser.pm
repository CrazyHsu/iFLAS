#!/bin/env perl
package samParser;

use strict;
use 5.010;
require Exporter;

1;

sub isPaired{
    my ($flag) = @_;
    return ($flag & 0b1) == 0b1 ? 1 : 0;
}

sub isProperPair{
    my ($flag) = @_;
    return ($flag & 0b11) == 0b11 ? 1 : 0;
}

sub isUnmapped{
    my ($flag) = @_;
    return ($flag & 0b100) == 0b100 ? 1 : 0;
}

sub isMateUnmapped{
    my ($flag) = @_;
    return ($flag & 0b1000) == 0b1000 ? 1 : 0;
}

sub isReverseStrand{
    my ($flag) = @_;
    return ($flag & 0b10000) == 0b10000 ? 1 : 0;
}

sub isMateReverseStrand{
    my ($flag) = @_;
    return ($flag & 0b100000) == 0b100000 ? 1 : 0;
}

sub isFirstMate{
    my ($flag) = @_;
    return ($flag & 0b1000000) == 0b1000000 ? 1 : 0;
}

sub isSecondMate{
    my ($flag) = @_;
    return ($flag & 0b10000000) == 0b10000000 ? 1 : 0;
}

sub isSecondaryAlignment{
    my ($flag) = @_;
    return ($flag & 0b100000000) == 0b100000000 ? 1 : 0;
}

sub isFailedQC{
    my ($flag) = @_;
    return ($flag & 0b1000000000) == 0b1000000000 ? 1 : 0;
}

sub isDuplicate{
    my ($flag) = @_;
    return ($flag & 0b10000000000) == 0b10000000000 ? 1 : 0;
}

sub flag2Char{
    my ($flag) = @_;
    my $char = '';
    $char .= 'p' if isPaired($flag) == 1;
    $char .= 'P' if isProperPair($flag) == 1;
    $char .= 'u' if isUnmapped($flag) == 1;
    $char .= 'U' if isMateUnmapped($flag) == 1;
    $char .= 'r' if isReverseStrand($flag) == 1;
    $char .= 'R' if isMateReverseStrand($flag) == 1;
    $char .= '1' if isFirstMate($flag) == 1;
    $char .= '2' if isSecondMate($flag) == 1;
    $char .= 's' if isSecondaryAlignment($flag) == 1;
    $char .= 'f' if isFailedQC($flag) == 1;
    $char .= 'd' if isDuplicate($flag) == 1;
    return $char;
}

sub determineCodingStrand{
    my ($libType, $flag) = @_;
    if($libType eq "fr-unstranded"){
        return '.';
    }elsif($libType eq 'fr-firststrand'){
        my $charFlag = &flag2Char($flag);
        if($charFlag !~/[12]/){# single end
            return '+' if $charFlag =~ /r/;
            return '-' if $charFlag !~ /[ur]/;
        }else{
            return '+' if( ($charFlag =~ /1/ && $charFlag =~ /r/ && $charFlag !~ /R/) ||
            ($charFlag =~ /2/ && $charFlag !~ /[ur]/ && $charFlag =~ /[RU]/)  );
            return '-' if( ($charFlag =~ /1/ && $charFlag !~ /[ur]/ && $charFlag =~ /[RU]/) ||
           ($charFlag =~ /2/ && $charFlag =~ /r/ && $charFlag !~ /R/)    );
        }            
        return '';
    }elsif($libType eq 'fr-secondstrand'){
        my $charFlag = &flag2Char($flag);
        if($charFlag !~/[12]/){# single end
            return '+' if $charFlag !~ /[ur]/;
            return '-' if $charFlag =~ /r/;
        }else{
            return '+' if( ($charFlag =~ /1/ && $charFlag !~ /[ur]/ && $charFlag =~ /[RU]/) ||
           ($charFlag =~ /2/ && $charFlag =~ /r/ && $charFlag !~ /R/)    );
            return '-' if( ($charFlag =~ /1/ && $charFlag =~ /r/ && $charFlag !~ /R/) ||
           ($charFlag =~ /2/ && $charFlag !~ /[ur]/ && $charFlag =~ /[RU]/)    );
        }
        return '';
    }
}

sub cigar2multiJunction{
    my ($readStart, $cigar) = @_;
    my ($juncStart, $juncEnd);
    my $blockStart = $readStart;
    my (@juncStarts, @juncEnds);
    while($cigar =~ s/([^N]+?)(\d+)N//){
        my ($match1, $match2) = ($1, $2);
        my $blockSize = 0;
        $blockSize += $_ for($match1 =~ /(\d+)[MD=X]/g);
        $juncStart = $blockStart + $blockSize;
        #$blockStart += $1 if $match1 =~ /^(\d+)D/;
        push @juncEnds, $blockStart;
        $blockStart = $juncStart + $match2;
        #$juncStart -= $1 if $match1 =~ /(\d+)D$/;
        push @juncStarts, $juncStart;
    }
    shift @juncEnds;
    my $lastBlockSize = 0;
    $lastBlockSize += $_ for($cigar =~ /(\d+)[MD=X]/g);
    my $readEnd = $blockStart + $lastBlockSize;
    #$blockStart += $1 if $cigar =~ /^(\d+)D/;
    push @juncEnds, $blockStart;
    return (\@juncStarts, \@juncEnds, $readEnd);
}

sub cigar2eachJunction{
    my ($readStart, $cigar) = @_;
    my ($juncStarts, $juncEnds, $readEnd) = &cigar2multiJunction(@_);
    my @junctions;
    my $start = $readStart;
    for(my $i = 0; $i < @$juncStarts - 1; $i++){
    my @junction = ($start, $juncStarts->[$i], $juncEnds->[$i], $juncStarts->[$i+1]);
    push @junctions, \@junction;
    $start = $juncEnds->[$i];
    }
    my @junction = ($start, $juncStarts->[-1], $juncEnds->[-1], $readEnd);
    push @junctions, \@junction;
    return @junctions;
}

sub cigar2Blocks{ # $readStart is 0-based
    my ($readStart, $cigar) = @_;
    my $blockStart = $readStart;
    my (@blockStarts, @blockEnds);
    while($cigar =~ s/([^N]+?)(\d+)N//){
    my ($match1, $match2) = ($1, $2);
    my $blockEnd = $blockStart;
        #$blockStart += $1 if $match1 =~ /^(\d+)D/;
        push @blockStarts, $blockStart;
    $blockEnd += $_ for($match1 =~ /(\d+)[MD=X]/g);
        $blockStart = $blockEnd + $match2;
        #$blockEnd -= $1 if $match1 =~ /(\d+)D$/;
        push @blockEnds, $blockEnd;
    }
    my $blockEnd = $blockStart;
    #$blockStart += $1 if $cigar =~ /^(\d+)D/;
    push @blockStarts, $blockStart;
    $blockEnd += $_ for($cigar =~ /(\d+)[MD=X]/g);
    push @blockEnds, $blockEnd;
    return (\@blockStarts, \@blockEnds);
}

sub cigar2BlocksWithClip{
    my ($readStart, $cigar) = @_;
    my ($blockStarts, $blockEnds) = &cigar2Blocks($readStart, $cigar);
    my ($cdsStart, $cdsEnd) = ($blockStarts->[0], $blockEnds->[-1]);
    $blockStarts->[0] -= $1 if $cigar =~ /^(\d+)[SH]/;
    $blockEnds->[-1] += $1 if $cigar =~ /(\d+)[SH]$/;
    return ($blockStarts, $blockEnds, $cdsStart, $cdsEnd);
    
#    my $cdsStart = $readStart;
#    my $blockStart = $cigar =~ /^(\d+)[SH]/ ? ($cdsStart - $1) : $cdsStart;
#    my (@blockStarts, @blockEnds);
#    push @blockStarts, $blockStart;
#    $blockStart = $cdsStart;
#    while($cigar =~ s/([^N]+?)(\d+)N//){
#    my ($match1, $match2) = ($1, $2);
#    my $blockEnd = $blockStart;
#    $blockEnd += $_ for($match1 =~ /(\d+)[MD=X]/g);
#        $blockStart = $blockEnd + $match2;
#        push @blockStarts, $blockStart;
#        $blockEnd -= $1 if $match1 =~ /(\d+)D^/;
#        push @blockEnds, $blockEnd;
#    }
#    my $cdsEnd = $blockStart;
#    $cdsEnd += $_ for($cigar =~ /(\d+)[MD=X]/g);
#    my $blockEnd = $cigar =~ /(\d+)[SH]$/ ? ($cdsEnd + $1) : $cdsEnd;
#    push @blockEnds, $blockEnd;
#    return (\@blockStarts, \@blockEnds, $cdsStart, $cdsEnd);
}

####    Argument    Type    Description
#       samLine     string  A tab-separated line of a sam read record

####    Return      Type    Description
#       isMapped    int     1 when read is mappable
#                           0 when unmapped
sub isMapped{ #deprecated
    my ($samLine) = @_;
    my $cigar = (split "\t", $samLine)[5];
    if($cigar eq "*"){
        return 0;
    }else{
        return 1;
    }
}

####    Argument    Type    Description
#       samLine     string  A tab-separated line of a sam read record

####    Return      Type    Description
#       strand      char    The strand on to the read is mapped. Null string when inferable                      
sub getStrandOrientation{#deprecated
    my ($samLine) = @_;
    my $flag = (split "\t", $samLine)[1];
    if(
            ($flag &  0b1111101) ==  0b1010001 ||#read1 map - (infer read2 map +)
            ($flag & 0b11111101) == 0b10100001   #read2 map + (infer read1 map -)
        ){
            return '+';
        }elsif(
                ( $flag &  0b1111101) ==  0b1100001 ||#read1 map to +
                ( $flag & 0b11111101) == 0b10010001   #read2 map to -
              ){
            return '-';
        }else{
            return '';
        }
}

####    Argument    Type    Description
#       samLine     string  A tab-separated line of a sam read record

####    Return      Type    Description
#       bestHitsN   int     Number of best hits
sub getBestHitsN{
    my ($samLine) = @_;
    chomp $samLine;
    if($samLine =~ /X0:i:(\d+)/){
        return $1;
    }else{
        return;
    }
}

####    Argument        Type    Description
#       samLine         string  A tab-separated line of a sam read record

####    Return          Type    Description
#       subBestHitsN    int     Number of best hits
sub getSubHitsN{
    my ($samLine) = @_;
    chomp $samLine;
    if($samLine =~ /X1:i:(\d+)/){
        return $1;
    }else{
        return;
    }
}

####    Argument    Type    Description
#       samLine     string  A tab-separated line of a sam read record

####    Return      Type    Description
#       minNM       int     NM value of the alternative hit with minimal NM value
sub getAltHitsMinNM{
    my ($samLine) = @_;
    return unless $samLine=~/XA:Z:(\S+)+/;
    my @altHitsNM = $1=~/(\d+);/g;
    my $minNM = $altHitsNM[0];
    for my $altHitNM (@altHitsNM){
        $minNM = $altHitNM if $altHitNM < $minNM;
    }
    $minNM;
}

sub getTagsString{
    my ($line) = @_;
    my @fields = split "\t", $line;
    return join "\t", @fields[11..$#fields]
}

sub getTagValue{
    my ($line, $tag) = @_;
    my $tags = &getTagsString($line);
    if($tags =~ /\b$tag:\w:([^\t]+)/){
        return $1;
    }else{
        return;
    }
}

sub getReadLength{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MISH=X]/g);
    $length;
}

sub getReadHardClippedLength{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MIS=X]/g);
    $length;
}

sub getAlignedLength{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MI=X]/g);
    $length;
}

sub getMatchLength{
    my ($md) = @_;
    my $length = 0;
    $length += $_ for($md =~ /(\d+)/g);
    $length;
}

sub getMisMatchSeq{ # including N bases
    my ($md) = @_;
    my @misMatchSeq;
    for my $seq ($md =~ /([\^ATCGN]+)/g){
        push @misMatchSeq, $seq unless $seq =~ /^\^/;
    }
    @misMatchSeq;
}

sub getNBaseCount{
    my ($md) = @_;
    my $length = 0;
    for my $seq ($md =~ /([\^ATCGN]+)/g){
        next if $seq =~ /^\^/;
        $length++ for($seq =~ /(N)/g);
    }
    $length
}

sub getMisMatchLength{
    my ($md) = @_;
    my @misMatchSeq = &getMisMatchSeq($md);
    my $length = 0;
    $length += length($_) for(@misMatchSeq);
    $length;
}

sub getIndelLength{
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[ID]/g);
    $length;
}

sub reverseStrand{
    my ($flag) = @_;
    if(&isReverseStrand($flag) == 0){
        $flag += 0x10;
    }else{
        $flag -= 0x10;
    }
    $flag;
}

sub reverseMateStrand{
    my ($flag) = @_;
    if(&isMateReverseStrand($flag) == 0){
        $flag += 0x20;
    }else{
        $flag -= 0x20;
    }
    $flag;
}
