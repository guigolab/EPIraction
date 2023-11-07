#!/usr/bin/perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $tissue       = "!{tissue}";
my $regions      = "!{regions}";
my $HiC_samples  = "!{HiC_samples}";
my $chromosomes  = "!{chromosomes}";
my $baseline_hic = "!{HiC_baseline}";
my $upstream     = "!{upstream}";
my $downstream   = "!{downstream}";
my $min_rescaled = 0.05;

die "wrong folder $data_folder\n" unless -d $data_folder;

my $HiC_baseline     = load_contacts($baseline_hic);
my $Upstream_limit   = load_contacts($upstream);
my $Downstream_limit = load_contacts($downstream);
$chromosomes     =~ s/\[|\]//g;

my %Regions = ();
open FILE_IN,"$data_folder/data/$tissue.regions.data";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;

while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my ($start,$end,$symbol,$activity,$strand) = @Data[1,2,11,13,5];
    $Regions{"$Data[3]"} = ["$start","$end","$activity","$symbol","$strand"];
}
close FILE_IN;

open FILE_OUT,"| gzip > Tissue.contacts.data.gz";
print FILE_OUT "chrom\tsymbol\tpromoter\tenhancer\tdistance\tactivity\tcontact\trescaled\tHiC-fold\n";

foreach my $chrom (sort split ", ", $chromosomes)
{
    my %Contacts = ();
    foreach my $HiC_sample (split ";", $HiC_samples)
    {
	open FILE_IN,"zcat $data_folder/temp/HiC/$HiC_sample.$chrom.overlap.gz |";
	while($string = <FILE_IN>)
	{
	    chomp $string;
	    next if $string eq "";
	    my ($promoter,$enhancer,$distance,$reads,$contact) = split "\t", $string;
	    unless(defined $Contacts{"$promoter<>$enhancer"})
	    {
		$Contacts{"$promoter<>$enhancer"} = $contact;
	    }
	    else
	    {
		$Contacts{"$promoter<>$enhancer"} = $contact if $Contacts{"$promoter<>$enhancer"} < $contact;
	    }
	}
	close FILE_IN;
    }
    print "Load tissues HiC for $chrom\n";

    my %Consensus = ();
    open FILE_IN,"zcat $data_folder/temp/HiC/Consensus.tissues.Encode.intact.$chrom.overlap.gz |";
    while($string = <FILE_IN>)
    {
	chomp $string;
	next if $string eq "";
	my ($promoter,$enhancer,$distance,$reads,$contact) = split "\t", $string;
	$Consensus{"$promoter<>$enhancer"} = $contact;
    }
    close FILE_IN;
    print "Load consensus HiC for $chrom\n";

    open FILE_IN,"zcat $data_folder/data/zzz.Pairs.$chrom.gz |";
    $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my ($promoter,undef,$list) = split "\t", $string;
	next unless defined $Regions{$promoter};
	my ($p_start,$p_end,$p_activity,$symbol,$strand) = @{$Regions{$promoter}};

	foreach my $enhancer (split ";", $list)
	{
	    next unless defined $Regions{$enhancer};
	    my ($e_start,$e_end,$e_activity,undef,undef) = @{$Regions{$enhancer}};
	    my @Sorted       = sort { $a <=> $b } ($p_start,$p_end,$e_start,$e_end);
	    my $distance     = $Sorted[2] - $Sorted[1];
	    my $distance_200 = get_200_distance($distance);
	    my $baseline     = $HiC_baseline->{$distance_200};
	    my $contact      = $baseline;
	    $contact         = $Consensus{"$promoter<>$enhancer"} if defined $Consensus{"$promoter<>$enhancer"} and $Consensus{"$promoter<>$enhancer"} >= $contact;
	    $contact         = $Contacts{"$promoter<>$enhancer"}  if defined $Contacts{"$promoter<>$enhancer"}  and $Contacts{"$promoter<>$enhancer"}  >= $contact;
	    my $hic_fold     = sprintf("%.3f",$contact/$baseline);

	    my $maximum      = $Upstream_limit->{$distance_200};
	    if($e_start > $p_end)
	    {
		if($strand eq "-")
		{
		    $maximum = $Upstream_limit->{$distance_200};
		}
		else
		{
		    $maximum = $Downstream_limit->{$distance_200};
		}
	    }
	    else
	    {
		if($strand eq "-")
		{
		    $maximum = $Downstream_limit->{$distance_200};
		}
		else
		{
		    $maximum = $Upstream_limit->{$distance_200};
		}
	    }
	    $contact         = $maximum if $contact > $maximum;
	    $contact         = sprintf("%.8f",$contact);
	    my $rescaled     = sprintf("%.6f",$e_activity*$contact);
	    next if $rescaled < $min_rescaled;
	    print FILE_OUT "$chrom\t$symbol\t$promoter\t$enhancer\t$distance\t$e_activity\t$contact\t$rescaled\t$hic_fold\n";
	}
    }
    close FILE_IN;
    print "Done $tissue.$chrom\n";
}
close FILE_OUT;
system "mv -f Tissue.contacts.data.gz $data_folder/data/$tissue.contacts.data.gz";
print "Done all\n";

sub load_contacts
{
    my $file  = $_[0];
    my %ret_h = ();
    open FILE_IN,$file or die "Absent $file";
    my $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[0]} = $Data[1];
    }
    close FILE_IN;
    return \%ret_h;
}

sub get_200_distance
{
    my $value = $_[0];
    return 1000000 if $value >= 999900;

    my $times = int(($value+100)/200);
    return $times*200;
}
