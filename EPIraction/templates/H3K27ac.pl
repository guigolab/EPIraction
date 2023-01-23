#!/usr/bin/env perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $samples_file = "!{samples_index}";
my $regions      = "!{regions}";
my $tissue       = "!{tissue}";
my $chromosomes  = "!{chromosomes}";

die "wrong folder $data_folder\n" unless -d $data_folder;
$chromosomes =~s/^\[//;
$chromosomes =~s/\]$//;

my ($samples,undef) = load_data("$data_folder/files/$samples_file");

my $cur_chrom = "Primer";
open FILE_IN,"$data_folder/files/$regions";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    unless($Data[0] eq $cur_chrom)
    {
	close REGIONS;
	print "Done split $cur_chrom\n";
	$cur_chrom = $Data[0];
	open REGIONS,">Regions.$cur_chrom.bed";
    }
    print REGIONS $string,"\n";
}
print "Done split $cur_chrom\n\n";
close REGIONS;
close FILE_IN;

my $prefix = "$data_folder/temp/signal/$tissue";
system "mkdir $prefix" unless -d $prefix;

foreach my $chromosome (split ", ", $chromosomes)
{
    if(-f "$prefix/Raw.$chromosome.data.gz")
    {
	print "already $prefix/Raw.$chromosome.data.gz";
	system "rm -f Binned.$chromosome.bed.gz Regions.$chromosome.bed";
	next;
    }

    my (@List,%Signal) = ();
    foreach my $sample_id (grep { $samples->{$_}->{"tissue"} eq $tissue } sort { $a <=> $b } keys %{$samples})
    {
	my $data      = $samples->{$sample_id};
	my $accession = $data->{"accession"};
	my $file      = $data->{"file"};
	my $bigWig    = "$data_folder/BigWig/$file.H3K27ac.bigWig";
	open FILE_IN, "bigWigAverageOverBed $bigWig Regions.$chromosome.bed -minMax stdout |";
	while(my $string = <FILE_IN>)
	{
	    chomp $string;
	    my @Data     = split "\t", $string;
	    my $max      = sprintf("%.6f",$Data[7]);
	    my $mean     = sprintf("%.6f",$Data[4]);
	    $Data[0] =~ s/\.\d+?$//;
	    $Signal{$Data[0]} .= "\t$max";
	}
	close FILE_IN;
	push @List,$accession;
    }
    print "Load $chromosome samples\n";

    open FILE_OUT,"| gzip > $tissue.$chromosome.data.gz";
    print FILE_OUT "region\t",(join "\t", @List),"\n";

    open FILE_IN,"Regions.$chromosome.bed";
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	print FILE_OUT $Data[3],$Signal{$Data[3]},"\n";
    }
    close FILE_OUT;
    system "mv $tissue.$chromosome.data.gz $prefix/Raw.$chromosome.data.gz";
    system "rm -f Regions.$chromosome.bed";
    print "Done $tissue.$chromosome.data.gz\n";
}
print "Done all\n";

sub load_data
{
    open FILE_L,$_[0] or die "$! $_[0]";
    my $string = <FILE_L>;chomp $string;$string=~s/\r//;
    my @Head = split "\t",$string;
    my $num_rows=$#Head;
    my %Hash_data=();
    while($string = <FILE_L>)
    {
        chomp $string;$string=~s/\r//;
        $string.="NA" if $string=~/\t$/;
        my @Content = split "\t",$string;
        %{$Hash_data{$Content[0]}}=();
        my $hash_ref = $Hash_data{$Content[0]};
        map { $hash_ref->{$Head[$_]}= $Content[$_]} (1..$num_rows);
    }
    close FILE_L;
    shift @Head;
    return (\%Hash_data,\@Head);
}
