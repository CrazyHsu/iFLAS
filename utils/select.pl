#!/usr/bin/env perl

use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;

my ($sep) = ("\t");
my ($indexes, $names, $header, $value);
sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.tsv >OUTPUT.tsv
    If INPUT.tsv isn't specified, input from STDIN
Option:
    -i --index  STRs    Comma-separated list specifying the fields to output
                        The element of the list can be a single column number or a range with nonnumeric char as separator
                        To specify the last column left the range right margin blank
                        If continuous range specified like '1-3-6', the first range '1-3' will be output
                        eg.:
                            -i 1,4          output columns 1,4
                            -i 1-4,6..8     output columns 1,2,3,4,6,7,8
                            -i 1,4,6-       output columns 1,4,6,7,... last column
                            -i 1-3          output columns 1,2,3
    -n --name   STRs    Comma-separated names to be select columns
                        Columns of the first line in INPUT.tsv will bed treated as names
       --header         Print the header as the first line (commented) to OUTPUT.tsv
    -s --sep    STR     The seperator to join fields[\t]
    -v --value  STR     The value to perform the select on the value joined by the selected columns
    -h --help           Print this help information
HELP
    exit(-1);
}

GetOptions(
            'i|index=s'     => \$indexes,
            'n|name=s'      => \$names,
            'header'        => \$header,
            's|sep=s'       => \$sep,
            'v|value=s'     => \$value,
            'h|help'        => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";

my @indexes;
if(defined $names){
    my $headers = <IN>;
    chomp $headers;
    my @headers = split "\t", $headers;
    my %header2index;
    $header2index{$headers[$_]} = $_+1 for(0..$#headers);
    my @names = split ",", $names;
    for my $name(@names){
        die "No column named as $name" if !exists $header2index{$name};
        push @indexes, $header2index{$name};
    }
    if(defined $header){
        print "#";
        say join "\t", @names;
    }
}elsif(!defined $indexes){
    die "--index or --name must be specified!\n";
}

while(<IN>){
    my @fields = split "\t";
    $fields[$#fields] =~ s/\n$//;
    if(!defined $names){
        @indexes = map{$_+1}(@{&common::arraySplice(\@fields, $indexes)->{index}});
    }
    my $output = join $sep, map{$fields[$_-1]}(@indexes);
    if(defined $value){
        next if $output ne $value;
    }
    say $output;
}
