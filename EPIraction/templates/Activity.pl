#!/usr/bin/perl -w

use strict;
use warnings;

my $data_folder   = "!{data_folder}";
my $tissue        = "!{tissue}";

die "wrong folder $data_folder\n" unless -d $data_folder;

my @Regions = ();
my ($open,$cofactor,$regions) = (0,0,0);
open FILE_IN,"$data_folder/temp/tissues/$tissue.regions.data" or die "absent $tissue\n";
my $header = <FILE_IN>;
chomp $header;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    push @Regions,\@Data;
    if($Data[6] eq "enhancer")
    {
	$open     += $Data[8];
	$cofactor += $Data[9];
	$regions++;
    }
}
close FILE_IN;
$open     /= $regions;
$cofactor /= $regions;

open FILE_OUT,">tissue.regions.data";
print FILE_OUT $header,"\tOpen.activity\tCof.activity\n";

foreach my $region (@Regions)
{
    print FILE_OUT join "\t",@{$region};
    print FILE_OUT "\t",sprintf("%.6f\t%.6f\n",$region->[8]/$open, $region->[9]/$cofactor);
}
close FILE_OUT;
system "mv tissue.regions.data $data_folder/data/$tissue.regions.data";
