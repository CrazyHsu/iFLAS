#!/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;

my ($scoreFile);
my ($clusterLength) = (15);
GetOptions(
            's|score=s'     => \$scoreFile,
            'l|length=i'    => \$clusterLength,
            'h|help'        => sub{usage()}
        ) || usage();

open SCO, $scoreFile or die "Can't open $scoreFile: $!";
$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";

my %juncs;
while(<SCO>){
    chomp;
    my ($chr, $strand, $start, $end, $lScore, $jScore, $rScore) = split "\t";
    if(exists $juncs{$chr}{$strand}){
        push @{$juncs{$chr}{$strand}}, [$start, $end, $lScore, $jScore, $rScore];
    }else{
        $juncs{$chr}{$strand} = [[$start, $end, $lScore, $jScore, $rScore]];
    }
}

my %alignHash;
while(<IN>){
    chomp;
    next if /^#/;
    my ($chr, $start, $name, $strand, $blockSizes, $blockRelStarts) = (split "\t")[0, 1, 3, 5, 10, 11];
    my @blockSizes = split ",", $blockSizes;
    my @blockRelStarts = split ",", $blockRelStarts;
    my ($juncStarts, $juncEnds) = &common::getIntrons(&common::getAbsLoc($start, \@blockSizes, \@blockRelStarts));
    if(@$juncStarts == 0){ # single-exon read
        print "$_\n";
        next;
    }
    my %juncs = map{("$juncStarts->[$_]-$juncEnds->[$_]", '')}0..(@$juncStarts-1);
    $alignHash{$chr}{$strand}{$name} = {    line    => $_,
                                            juncs   => \%juncs};
}

for my $chr (keys %juncs){
    my $chrV = $juncs{$chr};
    for my $strand (keys %$chrV){
        my $strandV = $chrV->{$strand};
        my @sorted = sort {$a->[0]<=>$b->[0]}@$strandV;
        my @juncs = ($sorted[0]);
        for(my $i = 1; $i < @sorted; $i++){
            if($sorted[$i]->[0] - $sorted[$i-1]->[0] <= $clusterLength){
                push @juncs, $sorted[$i];
            }else{
                &endClusterAndRefine(\@juncs, $chr, $strand);
                @juncs = ($sorted[$i]);
            }
        }
        &endClusterAndRefine(\@juncs, $chr, $strand);
    }
}

for my $chr (keys %alignHash){  # output refined alignment
    my $chrV = $alignHash{$chr};
    for my $strand (keys %$chrV){
        my $strandV = $chrV->{$strand};
        for my $name (keys %$strandV){
            my ($line, $juncs) = @{$strandV->{$name}}{qw/line juncs/};
            my @fields = split "\t", $line;
            my ($start, $end) = @fields[1, 2];
            my (@juncStarts, @juncEnds);
            for my $junc(keys %$juncs){
                my ($juncStart, $juncEnd) = split "-", $junc;
                push @juncStarts, $juncStart;
                push @juncEnds, $juncEnd;
            }
            @juncStarts = sort{$a<=>$b}@juncStarts;
            @juncEnds = sort{$a<=>$b}@juncEnds;
            if($juncStarts[0] <= $start){ # junction start is revised across the read start
                shift @juncStarts;
                next if @juncStarts == 0;
                $start = shift @juncEnds;
            }
            if($juncEnds[-1] >= $end){
                pop @juncEnds;
                next if @juncEnds == 0;
                $end = pop @juncStarts;
            }
            my @blockStarts = ($start, @juncEnds);
            my @blockEnds = (@juncStarts, $end);
            my $blockSizes = join ",", &common::getSizes(\@blockStarts, \@blockEnds);
            my $blockRelStarts = join ",", &common::getRelStarts(\@blockStarts);
            @fields[1, 2, 6, 7, 9, 10, 11] = ($start, $end, $start, $end, $#blockStarts + 1, $blockSizes, $blockRelStarts);
            print join "\t", @fields;
            print "\n";
        }
    }
}

sub endClusterAndRefine{
    my ($juncs, $chr, $strand) = @_;
    my @sorted = sort {$a->[1]<=>$b->[1]}@$juncs;   # sort by end
    my @juncs = ($sorted[0]);
    for(my $i = 1; $i < @sorted; $i++){
        if($sorted[$i]->[1] - $sorted[$i-1]->[1] <= $clusterLength){
            push @juncs, $sorted[$i];
        }else{
            &refine(\@juncs, $chr, $strand);
            @juncs = ($sorted[$i]);
        }
    }
    &refine(\@juncs, $chr, $strand);
}

sub refine{
    my ($juncs, $chr, $strand) = @_;
    my $highJScore = 0;
    my (@realJuncStarts, @realJuncEnds);
    for my $junc (@$juncs){
        my ($start, $end, $lScore, $jScore, $rScore) = @$junc;
        if($jScore ne ''){
            if($jScore > $highJScore){ # run into the higher jScore
                @realJuncStarts = ($start);
                @realJuncEnds = ($end);
                $highJScore = $jScore;
            }elsif($jScore == $highJScore){
                push @realJuncStarts, $start;
                push @realJuncEnds, $end;
            }
        }
    }
    if($highJScore == 0 || @realJuncStarts > 1){ # no junction supprots or multi equal jScore junctions
        my $highLScore = 0;
        @realJuncStarts = ();
        @realJuncEnds = ();
        for my $junc (@$juncs){
            my ($start, $end, $lScore, $jScore, $rScore) = @$junc;
            if($lScore ne ''){
                if($lScore > $highLScore){
                    @realJuncStarts = ($start);
                    @realJuncEnds = ($end);
                    $highLScore = $lScore;
                }elsif($lScore == $highLScore){
                    push @realJuncStarts, $start;
                    push @realJuncEnds, $end;
                }
            }
        }
        if(@realJuncStarts > 1){ # has multiple high score junction in long reads, try the reference supporting (rScore)
            @realJuncStarts = ();
            @realJuncEnds = ();
            my @refJuncs = grep{$_->[4] ne ''}@$juncs;
            if(@refJuncs == 1){ # junction in annotaion is uniq
                @realJuncStarts = ($refJuncs[0]->[0]);
                @realJuncEnds = ($refJuncs[0]->[1]);
            }
        }
    }
    if(@realJuncStarts == 1){ # only one consensus result in three types of supporting
    JUNC:for my $junc (@$juncs){ # refine junctions
            my ($start, $end, $lScore, $jScore, $rScore) = @$junc;
            next if $jScore ne '' || $rScore ne ''; # skip junction with j or r support
            next if $start == $realJuncStarts[0] && $end == $realJuncEnds[0];
            my $strandV = $alignHash{$chr}{$strand};
            for my $name (keys %$strandV){
                my ($line, $juncs) = @{$strandV->{$name}}{qw/line juncs/};
                if(exists $juncs->{"$start-$end"}){
                    next JUNC if &isExonBroken($juncs, $realJuncStarts[0], $realJuncEnds[0]) == 1;
                    delete $juncs->{"$start-$end"};
                    $juncs->{"$realJuncStarts[0]-$realJuncEnds[0]"} = '';
                }
            }
        }
    }
}

sub isExonBroken(){
    my ($juncs, $newJuncStart, $newJuncEnd) = @_;
    my (@juncStarts, @juncEnds);
    for my $junc(keys %$juncs){
        my ($juncStart, $juncEnd) = split '-', $junc;
        push @juncStarts, $juncStart;
        push @juncEnds, $juncEnd;
    }
    @juncStarts = sort{$a<=>$b}@juncStarts;
    @juncEnds = sort{$a<=>$b}@juncEnds;
    for(my $i = 1; $i < @juncStarts; $i++){ # test whether is start out-of-exon_boundary
        if( $newJuncStart <= $juncEnds[$i-1] && $newJuncEnd > $juncStarts[$i]){
            return 1
        }
    }
    for(my $i = 0; $i < @juncStarts - 1; $i++){ # test whether is end out-of-exon_boundary
         if( $newJuncEnd >= $juncStarts[$i+1] && $newJuncStart < $juncEnds[$i]){
            return 1
         }
    }
    return 0;
}
sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName gmap.bed12+ >refined.bed12+
    If INPUT isn't specified, input from STDIN
Options:
    -s --score  FILE    The score file generated by juncScoring.pl
    -l --length INT     The max length between two adjacent splicing sites to be clustered[15]
    -h --help           Print this help information
HELP
    exit(-1);
}