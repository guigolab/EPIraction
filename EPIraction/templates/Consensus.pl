#!/usr/bin/env perl

use strict;
use warnings;

my $data_folder = "!{data_folder}";
my $regions     = "!{regions}";
my $tissues     = "!{tissues}";
my $chromosome  = "!{chromosome}";

die "wrong folder $data_folder\n" unless -d $data_folder;
$tissues =~s/^\[//;
$tissues =~s/\]$//;

if(-f "$data_folder/temp/normalize/Consensus.$chromosome.raw.data.gz")
{
    print "already Consensus.$chromosome.raw.data.gz\n";
    exit;
}

open FILE_IN,"$data_folder/files/$regions" or die $!;
open REGIONS,">Current_regions.bed";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next unless $Data[0] eq $chromosome;
    print REGIONS $string,"\n";
}
close FILE_IN;
close REGIONS;

my (@List,%Signal) = ();
foreach my $tissue (split ", ", $tissues)
{
    my $bigWig = "$data_folder/BigWig/$tissue.H3K27ac.consensus.bigWig";
    open FILE_IN, "bigWigAverageOverBed $bigWig Current_regions.bed -minMax stdout |";
    my %Maximum = ();
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	my $max  = sprintf("%.6f",$Data[7]);
	my $mean = sprintf("%.6f",$Data[4]);
	$Signal{$Data[0]}.="\t$max";
    }
    close FILE_IN;
    push @List,$tissue;
}
print "Load all\n";

open FILE_OUT,"| gzip > Consensus.$chromosome.data.gz";
print FILE_OUT "region\t",(join "\t", @List),"\n";

open FILE_IN,"Current_regions.bed";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    print FILE_OUT $Data[3],$Signal{$Data[3]},"\n";
}
close FILE_OUT;
system "mv Consensus.$chromosome.data.gz $data_folder/temp/normalize/Consensus.$chromosome.raw.data.gz";
system "rm -f Current_big.bed Current_regions.bed";
print "Done matrix\n";

sub get_nearest
{
    my $value = $_[0];
    my $times = int(($value+25)/50);
    return $times*50;
}
