package common;

use strict;
use 5.010;
require Exporter;
use List::Util qw/sum/;
#use DBI;

#our @ISA = qw/Exporter/;
#our @EXPORT = qw/getCDSLength/;
1;

sub getExonsLength(){#debugged
    my($exonStarts, $exonEnds) = @_;
    my @exonLengths = map{$exonEnds->[$_]-$exonStarts->[$_]}(0..scalar(@$exonStarts)-1);
    return sum(@exonLengths);
}

sub getOverlap(){
    my ($blockStarts1, $blockEnds1, $blockStarts2, $blockEnds2) = @_;
    my @blockStarts1 = @$blockStarts1;
    my @blockEnds1 = @$blockEnds1;
    my @blockStarts2 = @$blockStarts2;
    my @blockEnds2 = @$blockEnds2;
    my $j = 0;
    my (@overlapStarts, @overlapEnds);
I:  for(my $i = 0; $i <= $#blockStarts1; $i++){
	my ($blockStart1, $blockEnd1) = ($blockStarts1[$i], $blockEnds1[$i]);
	for(;$j <= $#blockStarts2; $j++){
	    my ($blockStart2, $blockEnd2) = ($blockStarts2[$j], $blockEnds2[$j]);
	    if($blockStart2 >= $blockEnd1){ # 2 is at right of 1
		next I;
	    }elsif($blockEnd2 <= $blockStart1){
		next;
	    }else{ # intersect
		if($blockStart2 <= $blockStart1){
		    if($blockEnd2 <= $blockEnd1){
			push @overlapStarts, $blockStarts1[$i];
			push @overlapEnds, $blockEnds2[$j];
		    }else{
			push @overlapStarts, $blockStarts1[$i];
			push @overlapEnds, $blockEnds1[$i];
			next I;
		    }
		}else{
		    if($blockEnd2 <= $blockEnd1){
			push @overlapStarts, $blockStarts2[$j];
			push @overlapEnds, $blockEnds2[$j];
		    }else{
			push @overlapStarts, $blockStarts2[$j];
			push @overlapEnds, $blockEnds1[$i];
			next I;
		    }
		}
	    }
	}
    }
    return (\@overlapStarts, \@overlapEnds);
}

sub getOverlapLength(){
    &getExonsLength(&getOverlap(@_));
}

sub getIntrons{
    my ($exonStarts, $exonEnds) = @_;
    my @exonStarts = @$exonStarts;
    my @exonEnds = @$exonEnds;
    shift @exonStarts;	#change to be intron ends
    pop @exonEnds;	#change to be intron starts
    return (\@exonEnds, \@exonStarts);
}

sub isCanonicalSite{
    my ($chr, $start, $strand, $type, $seqHash) = @_;
    return unless exists $seqHash->{$chr};
    return if $start + 2 > length($seqHash->{$chr});
    my $siteSeq = substr ($seqHash->{$chr}, $start, 2);
    $siteSeq = uc $siteSeq;
    if( $strand eq '+'){
        if($type eq 'donor'){
            return 1 if $siteSeq =~ /(G[TC])|(AT)/g;
        }else{ # acceptor
            return 1 if $siteSeq =~ /A[GC]/g;
        }
    }else{
        if($type eq 'donor'){
            return 1 if $siteSeq =~ /([AG]C)|(AT)/g;
        }else{ # acceptor
	    
            return 1 if $siteSeq =~ /[CG]T/g;
        }
    }
    return 0;
}

sub uniqArray{
    my ($array) = @_;
    my %uniqArray;
    for my $value (@$array){
        $uniqArray{"$value"} = '';
    }
    keys %uniqArray;
}

sub removeUndefElement(){
    my ($array) = @_;
    my @array = @$array;
    my @returnArray;
    for(my $i = 0; $i <= $#array; $i++){
	push @returnArray, $array[$i] if defined $array[$i];
    }
    return @returnArray;
}

####    Argument    Type    Description
#       fileHandle  FH      A file handle of gpe file

####    Return      Type    Description
#       hash        hash    A hash with chr, strand as key and with array containing [start, end, name] as values
sub buildHash{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($name, $chr, $strand, $start, $end) = @fields[0..4];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, $name];
        }else{
            $hash{$chr}{$strand} = [ [$start, $end, $name] ];
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
#       fileHandle  FH      A file handle of gpe file

####    Return      Type    Description
#       hash        hash    A hash with chr, strand as key and with array containing [start, end, line] as values
sub buildHash2{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($chr, $strand, $start, $end) = @fields[1..4];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, (join "\t", @fields)];
        }else{
            $hash{$chr}{$strand} = [ [$start, $end, (join "\t", @fields)] ];
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

sub buildHash3{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($name, $chr, $start, $end) = @fields[0, 1, 3, 4];
	if(exists $hash{"$chr"}){
	    push @{$hash{"$chr"}}, [$start, $end, (join "\t", @fields)];
	}else{
	    $hash{"$chr"} = [ [$start, $end, (join "\t", @fields)] ];
	}        
    }
    for my $chr(keys %hash){
	my @sorted = sort{$a->[0]<=>$b->[0]}(@{$hash{$chr}});
	$hash{"$chr"} = \@sorted;
    }
    return %hash;
}

sub buildHash4{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($chr, $strand) = @fields[1, 2];
	$hash{$chr}{$strand} = $_;
    }
    return %hash;
}

sub buildHash5{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($chr, $strand, $start, $end) = @fields[1..4];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, \@fields];
        }else{
            $hash{$chr}{$strand} = [ [$start, $end, \@fields] ];
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

sub buildGpeHashByTransName{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my $name = $fields[0];
	$hash{"$name"} = $_;
    }
    return %hash;
}

sub buildGpeHashByGeneNameWithTransAsArray{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my $gene = $fields[11];
	push @{$hash{"$gene"}}, [$_];
    }
    return %hash;
}

sub buildGpeHashByGeneNameWithTransAsHash{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($name, $gene) = @fields[0, 11];
	$hash{"$gene"}{"$name"} = $_;
    }
    return %hash;
}

sub buildBedHash{
    my ($fh) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my ($chr, $start, $end, $strand, $blockSizes, $blockRelStarts) = (split "\t")[0..2, 5, 10, 11];
	my @blockSizes = split ",", $blockSizes;
	my @blockRelStarts = split ",", $blockRelStarts;
	my ($blockStarts, $blockEnds) = getAbsLoc($start, \@blockSizes, \@blockRelStarts);
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, $blockStarts, $blockEnds, $_];
        }else{
            $hash{$chr}{$strand} = [ [$start, $end, $blockStarts, $blockEnds, $_] ];
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

sub buildBedHash2{
    my ($fh) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my ($chr, $start, $end, $name, $strand) = (split "\t")[0..3, 5];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, $name];
        }else{
            $hash{$chr}{$strand} = [ [$start, $end, $name] ];
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

sub buildExonsHash{
    my ($fh, $bin) = @_;
    my (%hash, %auxHash);
    while(<$fh>){
	chomp;
	my @fields = split "\t";
	shift @fields if defined $bin;
	my ($chr, $strand, $blockStarts, $blockEnds) = @fields[1..4, 8, 9];
	my @blockStarts = split ',', $blockStarts;
	my @blockEnds = split ',', $blockEnds;
	for(my $i = 0; $i <= $#blockStarts; $i++){
	    next if exists $auxHash{"$chr:$strand:$blockStarts[$i]-$blockEnds[$i]"};
	    $auxHash{"$chr:$strand:$blockStarts[$i]-$blockEnds[$i]"} = '';
	    if(exists $hash{$chr}{$strand}){
		push @{$hash{$chr}{$strand}}, [$blockStarts[$i], $blockEnds[$i]];
	    }else{
		$hash{$chr}{$strand} = [ [$blockStarts[$i], $blockEnds[$i]] ];
	    }
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

sub buildBedExonsHash{
    my ($fh) = @_;
    my (%hash, %auxHash);
    while(<$fh>){
	chomp;
	my ($chr, $start, $strand, $blockSizes, $blockRelStarts) = (split "\t")[0, 1, 5, 10, 11];
	my ($blockStarts, $blockEnds) = getAbsLoc($start, $blockSizes, $blockRelStarts);
	for(my $i = 0; $i < @$blockStarts; $i++){
	    next if exists $auxHash{"$chr:$strand:$blockStarts->[$i]-$blockEnds->[$i]"};
	    $auxHash{"$chr:$strand:$blockStarts->[$i]-$blockEnds->[$i]"} = '';
	    if(exists $hash{$chr}{$strand}){
		push @{$hash{$chr}{$strand}}, [$blockStarts->[$i], $blockEnds->[$i]];
	    }else{
		$hash{$chr}{$strand} = [ [$blockStarts->[$i], $blockEnds->[$i]] ];
	    }
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

sub getConsensusIntronN(){
    my ($exonStarts1, $exonEnds1, $exonStarts2, $exonEnds2, $offset) = @_;
    $offset = 0 unless defined $offset;
    my ($intronStarts1, $intronEnds1) = getIntrons($exonStarts1, $exonEnds1);
    my ($intronStarts2, $intronEnds2) = getIntrons($exonStarts2, $exonEnds2);
    my ($j, $consensusN) = (0, 0);
I:  for(my $i = 0; $i < @$intronStarts1; $i++){
	for(; $j < @$intronStarts2; $j++){
	    next I if($intronStarts2->[$j] > $intronEnds1->[$i]);
	    if($intronStarts1->[$i] - $offset <= $intronStarts2->[$j] && $intronStarts2->[$j] <= $intronStarts1->[$i] + $offset &&
		$intronEnds1->[$i] - $offset <= $intronEnds2->[$j] && $intronEnds2->[$j] <= $intronEnds1->[$i] + $offset){
		$consensusN++;
	    }
	}
    }
    return $consensusN;
}

sub getCnsSite{
    my ($sites1, $sites2, $offset) = @_;
    $offset = 0 unless defined $offset;
    my ($s2, $cnsN) = (0, 0);
S:  for my $site1(@$sites1){
	for(; $s2 < @$sites2; $s2++){
	    my $site2 = $sites2->[$s2];
	    next S if $site2 > $site1 + $offset;
	    $cnsN++ if $site1 - $offset <= $site2 && $site2 <= $site1 + $offset;
	}
    }
    return $cnsN;
}

sub getAbsStarts{ # old name (getGpeBlockEnds) is deprecated
    my ($start, $blockRelStarts) = @_;
    my @blockRelStarts = @$blockRelStarts;
    map{$start + $_}@blockRelStarts;
}

sub getAbsLoc{
    my ($start, $blockSizes, $blockRelStarts) = @_;
    my @blockStarts = getAbsStarts($start, $blockRelStarts);
    my @blockSizes = @$blockSizes;
    my @blockEnds = map{$blockStarts[$_] + $blockSizes[$_]}0..$#blockStarts;
    return (\@blockStarts, \@blockEnds);
}

sub getSizes{
    my ($blockStarts, $blockEnds) = @_;
    map{$blockEnds->[$_] - $blockStarts->[$_]}0..(@$blockStarts-1);
}

sub getRelStarts{
    my ($blockStarts) = @_;
    map{$_ - $blockStarts->[0]}@$blockStarts;
}

sub getMergedTrans(){
    my ($blockStarts1, $blockEnd1, $blockStarts2, $blockEnd2) = @_;
    my @blockStarts1 = @$blockStarts1;
    my @blockEnds1 = @$blockEnd1;
    my @blockStarts2 = @$blockStarts2;
    my @blockEnds2 = @$blockEnd2;
    my (@mergedStarts, @mergedEnds);
    my %hash_mergedTran;
    my ($i,$j,$c,$d);
    my (@a, @b);
    $j=0;
    foreach (@blockStarts1){
	if($j>0&&$blockStarts1[$j]<$blockEnds1[$j-1]){
	    if($blockEnds1[$j]<$blockEnds1[$j-1])
	    {
		$blockStarts1[$j]=$blockStarts1[$j-1];
		$blockEnds1[$j]=$blockEnds1[$j-1];
	    }
	    else{
		$blockStarts1[$j]=$blockStarts1[$j-1];
		$hash_mergedTran{$blockStarts1[$j-1]}=$blockEnds1[$j];
	    }
	    $j++;
	    next;
	}
	elsif($j>0&&$blockStarts1[$j]==$blockEnds1[$j-1]){
	    $blockStarts1[$j]=$blockStarts1[$j-1];
	}
	my $mergedStart=$blockStarts1[$j];
	$i=0;
	if(@blockStarts2){
	    foreach (@blockStarts2){
		if($_<$mergedStart){
		    if($blockEnds2[$i]<$mergedStart){
			$hash_mergedTran{$blockStarts2[$i]}=$blockEnds2[$i];
			my $t=scalar @blockStarts2;
			if($blockEnds2[$t-1]<$blockStarts1[$j]){
			    $hash_mergedTran{$blockStarts1[$j]}=$blockEnds1[$j];
			}
		    }
		    elsif($blockEnds2[$i]==$mergedStart){
			$hash_mergedTran{$blockStarts2[$i]}=$blockEnds1[$j];
			$blockStarts1[$j]=$blockStarts2[$i];
		    }
		    else{
			if($blockEnds2[$i]<=$blockEnds1[$j]){
			    $hash_mergedTran{$blockStarts2[$i]}=$blockEnds1[$j];
			    $blockStarts1[$j]=$blockStarts2[$i];
			}
			else{
			    $hash_mergedTran{$blockStarts2[$i]}=$blockEnds2[$i];
			    $blockStarts1[$j]=$blockStarts2[$i];
			    $blockEnds1[$j]=$blockEnds2[$i];
			}
		    }
		}elsif($_>=$mergedStart){
		    if($_>$blockEnds1[$j]){
			$hash_mergedTran{$blockStarts1[$j]}=$blockEnds1[$j];
			last;
		    }
		    elsif($blockEnds2[$i]>=$blockEnds1[$j]){
			$hash_mergedTran{$blockStarts1[$j]} = $blockEnds2[$i];
			$blockEnds1[$j]=$blockEnds2[$i];
		    }
		    else{
			$hash_mergedTran{$blockStarts1[$j]} = $blockEnds1[$j];
		    }
		}
		$i++;											
	    }
	    if($i > 0){
		for(0..$i-1){
		    shift @blockStarts2;
		    shift @blockEnds2;
		}
	    }
	}else{
	    $hash_mergedTran{$blockStarts1[$j]} = $blockEnds1[$j];
	}							
	$j++;	
    }
    if(@blockStarts2){
	my $i = scalar @blockStarts2;
	foreach(0..$i-1){
	    $hash_mergedTran{$blockStarts2[$_]} = $blockEnds2[$_];
	}
    }
    @a= sort {$a <=> $b } keys %hash_mergedTran;
    @b= sort {$a <=> $b } values %hash_mergedTran;
    return (\@a,\@b);    
}

sub parseBlockChain{
    my ($chain) = @_;
    my (@blockStarts, @blockEnds);
    for my $block(split ';', $chain){
	my ($blockStart, $blockEnd) = split '-', $block;
	push @blockStarts, $blockStart;
	push @blockEnds, $blockEnd;
    }
    return (\@blockStarts, \@blockEnds)
}

sub is1Contain2{
    my ($exonStarts1, $exonEnds1, $exonStarts2, $exonEnds2) = @_;
    for(my $i1 = 0; $i1 < @$exonStarts1; $i1++){
	if($exonEnds2->[0] == $exonEnds1->[$i1]){
	    return 0 if $i1 != 0 && $exonStarts2->[0] < $exonStarts1->[$i1];
	    $i1++;
	    for(my $i2 = 1;$i2 < @$exonStarts2 - 1; $i2++){
		return 0 if $i1 >= @$exonStarts1; # out of exon index of trans1
		return 0 if $exonStarts2->[$i2] != $exonStarts1->[$i1] || $exonEnds2->[$i2] != $exonEnds1->[$i1];
		$i1++;
	    }
	    return 0 if $exonStarts2->[-1] != $exonStarts1->[$i1];
	    return 0 if $i1 != @$exonStarts1 - 1 && $exonEnds2->[-1] > $exonEnds1->[$i1];
	    return 1;
	}
    }
    return 0;
}

sub getTagsString(){
    my ($line) = @_;
    my @fields = split "\t", $line;
    return join "\t", @fields[11..$#fields]
}

sub getTagValue(){
    my ($line, $tag) = @_;
    my $tags = &getTagsString($line);
    if($tags =~ /\b$tag:\w:([^\t]+)/){
	return $1;
    }else{
	return;
    }
}

sub getMatchLength(){
    my ($md) = @_;
    my $length = 0;
    $length += $_ for($md =~ /(\d+)/g);
    $length;
}

sub getReadLength(){
    my ($cigar) = @_;
    my $length = 0;
    $length += $_ for($cigar =~ /(\d+)[MISH=X]/g);
    $length;
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

sub toQualValue(){
    my ($qual, $platform) = @_;
    my $offset = getQualOffset($platform);
    my @BQs;
    $qual =~ s/(.)/push @BQs, ord($1)-$offset/eg;
    return @BQs;
}

sub reverseComplement(){
    my ($seq) =@_;
    $seq = join '', reverse (split '', $seq);
    $seq =~ tr/ATCGRYKMBVDH/TAGCYRMKVBHD/;
    $seq =~ tr/atcgrykmbvdh/tagcyrmkvbhd/;
    $seq;
}

sub getLeftUTRExons(){ # debugged
    my ($CDSStart, $exonStarts, $exonEnds) = @_;
    my $exonCount = @$exonStarts;
    for( my $exonSIndex = 0; $exonSIndex < $exonCount; $exonSIndex++){
	if($exonStarts->[$exonSIndex] <= $CDSStart && $CDSStart <= $exonEnds->[$exonSIndex]){
	    my (@leftUTRExonStarts, @leftUTRExonEnds);
	    if($exonStarts->[$exonSIndex] == $CDSStart){
		$exonSIndex--;
		@leftUTRExonStarts = @{$exonStarts}[0..$exonSIndex];
		@leftUTRExonEnds = @{$exonEnds}[0..$exonSIndex];
	    }else{
		@leftUTRExonStarts = @{$exonStarts}[0..$exonSIndex];
		@leftUTRExonEnds = @{$exonEnds}[0..$exonSIndex];
		$leftUTRExonEnds[-1] = $CDSStart;
	    }
	    return (\@leftUTRExonStarts, \@leftUTRExonEnds);
	}elsif($CDSStart < $exonStarts->[$exonSIndex]){
	    return;
	}
    }
}

sub getRightUTRExons(){ #debugged
    my ($CDSEnd, $exonStarts, $exonEnds) = @_;
    my $maxIndex = @$exonStarts -1;
    for( my $i = $maxIndex; $i >= 0; $i--){
	if($exonStarts->[$i] < $CDSEnd && $CDSEnd <= $exonEnds->[$i]){
	    my (@UTRExonStarts, @UTRExonEnds);
	    if($exonEnds->[$i] == $CDSEnd){
		$i++;
		@UTRExonStarts = @{$exonStarts}[$i..$maxIndex];
		@UTRExonEnds = @{$exonEnds}[$i..$maxIndex];
	    }else{
		@UTRExonStarts = @{$exonStarts}[$i..$maxIndex];
		@UTRExonEnds = @{$exonEnds}[$i..$maxIndex];
		$UTRExonStarts[0] = $CDSEnd;
	    }
	    return (\@UTRExonStarts, \@UTRExonEnds);
	}elsif($CDSEnd > $exonEnds->[$i]){
	    return;
	}
    }
}

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
    if($libType == "fr-unstranded"){
        return '.';
    }elsif($libType == 'fr-firststrand'){
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
    }elsif($libType == 'fr-secondstrand'){
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

sub dichotomy{
    my ($qPos, $tPoss, $toSortRegion) = @_; # pos is in 1-based
    my @tPoss = @$tPoss;
    return ('n') if @tPoss == 0;
    if(defined $toSortRegion && $toSortRegion == 1){
        @tPoss = sort {$a <=> $b}@tPoss;
    }
    return ('l') if $qPos < $tPoss[0];
    return ('r') if $qPos > $tPoss[-1];
    
    my ($l, $r) = (0, $#tPoss);
    while($r - $l > 1){
        my $m = int (($l + $r)/2);
        if($qPos < $tPoss[$m]){
            $r = $m;
        }elsif($qPos > $tPoss[$m]){
            $l = $m;
        }else{
	    my @indexes;
	    for(my $i = $m; $i >= 0 && $tPoss[$i] == $qPos; $i--){
		push @indexes, $i;
	    }
	    @indexes = reverse @indexes;
	    for(my $j = $m + 1; $j <= $#tPoss && $tPoss[$j] == $qPos ; $j++){
		push @indexes, $j;
	    }
	    return ('e', @indexes);
	}
    }
    if($tPoss[$l] == $qPos){
	return ('e', $l);
    }elsif($tPoss[$r] == $qPos){
	return ('e', $r);
    }else{
	return ('b', $l, $r);
    }
}

sub getOvlRegs{
    my ($qStart, $qEnd, $tRegions, $toSortRegion) = @_;
    return if !defined $tRegions;
    my @tRegions = @$tRegions;
    return if @tRegions == 0;
    if(defined $toSortRegion && $toSortRegion == 1){
        @tRegions = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]}@tRegions;
    }
    my @targetStarts = map{$tRegions[$_]->[0]}0..$#tRegions;
    my ($status, @indexes) = &dichotomy($qEnd, \@targetStarts);
    my $rightIndex;
    if($status eq 'e'){ # e: equal
	$rightIndex = $indexes[0] - 1;
    }elsif($status eq 'b'){ # b: between
	$rightIndex = $indexes[0];
    }elsif($status eq 'r'){ # r: right
	$rightIndex = $#tRegions;
    }else{
	return;
    }
    my @ovlRegions;
    my $qLen = $qEnd - $qStart;
    for(my $i = 0; $i <= $rightIndex; $i++){
	my ($tStart, $tEnd) = @{$tRegions[$i]};
	if($qStart < $tEnd && $qEnd > $tStart){
	    my $ovlLen = ($qEnd < $tEnd ? $qEnd : $tEnd) - ($qStart > $tStart ? $qStart : $tStart);
	    my $ovlRatioInQuery = $ovlLen / $qLen;
	    my $ovlRatioInTarget = $ovlLen / ($tEnd - $tStart);
	    push @ovlRegions, [$tRegions[$i], $ovlLen, $ovlRatioInQuery, $ovlRatioInTarget];
	}
    }
    return @ovlRegions;
}

sub arraySplice{
    my ($arrayI,$splices)=@_;
    my @splicesA=split ",",$splices;
    my @indexs;
    for my $splice (@splicesA){
        if($splice=~/\D+/){#-f 1- or -7 or 1-4 or -f 4-1 or -f1-1
            if( $splice=~/\D+$/){#-f 1-
                my $from=(split /\D+/,$splice)[0];
                die "'$splice' in the fields your specify ($splices) isn't in correct form" unless defined $from;
                push @indexs,$_ for( $from-1..$#$arrayI );
            }elsif( $splice=~/^\D+/ ){#-f -7
                die "'$splice' in the fields your specify ($splices) isn't in correct form" if split /\D+/,$splice !=2;
                my $from=(split /\D+/,$splice)[1];
                push @indexs,$_ for reverse( $from-1..$#$arrayI );
            }else{#-f 1-4 or -f 4-1 or -f 1-1
                my ($from,$to)=split /\D+/,$splice;
                if($from<$to){#-f 1-4
                   if ($to>@$arrayI){
                        say STDERR "Warnning: no $to columns in: ".join "\t",@$arrayI;
                        $to=@$arrayI;
                    }
                    push @indexs,$_ for( ($from-1)..($to-1) ); 
                }else{#-f 4-1 or -f 1-1
                    ($from,$to)=($to,$from);
                   if ($to>@$arrayI){
                        say STDERR "Warnning: no $to columns in: ".join "\t",@$arrayI;
                        $to=$#$arrayI;                        
                    }
                    push @indexs,$_ for reverse( ($from-1)..($to-1) ); 
                }    
            }
        }else{#-f 1
            if ($splice>@$arrayI){
                say STDERR "Warnning: no $splice columns in: ".join "\t",@$arrayI;
            }else{
                push @indexs,($splice-1);
            }            
        }
    }
    my @arraySpliced=map{ $arrayI->[$_] }@indexs;
    return {
                "index" => \@indexs,
                "array" => \@arraySpliced
            };
}

