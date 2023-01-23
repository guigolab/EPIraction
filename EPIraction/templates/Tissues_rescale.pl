#!/usr/bin/env perl

use strict;

my $data_folder = "!{data_folder}";
my $chromosome  = "!{chromosome}";

die "wrong folder $data_folder\n" unless -d $data_folder;

open FILE_RAW,"zcat $data_folder/temp/normalize/Consensus.$chromosome.raw.data.gz |";
open FILE_NORM,"zcat $data_folder/temp/normalize/Consensus.$chromosome.norm.data.gz |";
my $raw_string  = <FILE_RAW>;
my $norm_string = <FILE_NORM>;
unless($raw_string eq $norm_string)
{
    print "Wrong order:\n$raw_string$norm_string";
    die(127)
}
open FILE_OUT,"| gzip > Current.scale.data.gz";
print FILE_OUT $raw_string;

while()
{
    $raw_string  = <FILE_RAW>;
    last unless $raw_string;
    $norm_string = <FILE_NORM>;
    chomp $raw_string;
    chomp $norm_string;
    my @Raw_data  = split "\t", $raw_string;
    my @Norm_data = split "\t", $norm_string;
    my $region = shift @Raw_data;shift @Norm_data;

    my @Report = ();
    push @Report, $region;
    foreach my $index (0..$#Raw_data)
    {
	if($Raw_data[$index] < 0.1 or $Norm_data[$index] < 0.1)
	{
	    push @Report,"1.00";
	}
	else
	{
	    my $scale = $Norm_data[$index]/$Raw_data[$index];
	    if($scale < 0.25)
	    {
		push @Report,"0.25";
	    }
	    elsif($scale > 5)
	    {
		push @Report,"5.00";
	    }
	    else
	    {
		push @Report,sprintf("%.6f",$scale);
	    }
#	    print "Raw\t",$Raw_data[$index],"\nNorm\t",$Norm_data[$index],"\nScale\t$scale\n";exit;
	}
    }
    print FILE_OUT join "\t", @Report;
    print FILE_OUT "\n";
}
close FILE_RAW;
close FILE_NORM;
close FILE_OUT;
system "mv Current.scale.data.gz $data_folder/temp/normalize/Consensus.$chromosome.scale.data.gz";
print "Done $chromosome\n";
