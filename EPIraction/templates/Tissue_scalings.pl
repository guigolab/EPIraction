#!/usr/bin/perl

use strict;
use warnings;

my $data_folder   = "!{data_folder}";
my $tissue        = "!{tissue}";
my $regions       = "!{regions}";
my $HiC_samples   = "!{HiC_samples}";
my $chromosomes   = "!{chromosomes}";

die "wrong folder $data_folder\n" unless -d $data_folder;

if(-f "$data_folder/data/$tissue.scalings.data.gz")
{
    print "Already $data_folder/data/$tissue.scalings.data.gz\n";
    exit;
}

my $HiC_baseline     = load_baseline();
my $distance_penalty = load_distance_penalty();
my $genes_data       = load_gene_data();

$chromosomes =~ s/\[|\]//g;

my %Active_regions = ();
open FILE_IN,"$data_folder/data/$tissue.regions.data";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;

while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my ($start,$end,$H3K27ac,$symbol) = @Data[1,2,7,11];
    my $open_activity = $Data[13];
    $open_activity    = 2 if $open_activity > 2;
    my $cofactor_activity = $Data[14];
    my $activity = sprintf("%.8f",$open_activity*$cofactor_activity);
    $activity = "3.0000000" if $activity > 3;
    $Active_regions{"$Data[3]"} = ["$start","$end","$activity","$H3K27ac","$symbol"];
#    map { print $_,"\t",$Header[$_],"\t",$Data[$_],"\n" } (0..$#Data);exit;
}
close FILE_IN;

open FILE_OUT,"| gzip > Tissue.scalings.data.gz";
print FILE_OUT "chrom\tsymbol\tpromoter\tenhancer\tdistance\tactivity\tcontact\tdistance_penalty\tscaling\trescaled\toriginal\n";

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
	    my ($promoter,$enhancer,undef,$contact,undef,undef) = split "\t", $string;
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

    open FILE_IN,"zcat $data_folder/data/zzz.Pairs.$chrom.gz |";
    $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my ($promoter,undef,$list) = split "\t", $string;
	next unless defined $genes_data->{$promoter};
	my ($p_start,$p_end,$symbol,$strand) = @{$genes_data->{$promoter}};

	my @Enhancers = ();
	foreach my $enhancer (split ";", $list)
	{
	    next unless defined $Active_regions{$enhancer};
	    my ($e_start,$e_end,$activity,$H3K27ac,undef) = @{$Active_regions{$enhancer}};
	    if($e_start > $p_end)
	    {
		my $distance = $e_start - $p_end;
		my $HiC_distance = get_HiC_distance($distance);
		my $contact      = $HiC_baseline->{$HiC_distance};
		$contact         = $Contacts{"$promoter<>$enhancer"} if defined $Contacts{"$promoter<>$enhancer"} and $Contacts{"$promoter<>$enhancer"} >= $contact;

		my $distance_200 = get_200_distance($distance);
		my $dist_penalty = "NA";

		if($strand eq "-")	{
		    $dist_penalty = $distance_penalty->{$distance_200}->[0];
		}
		else	{
		    $dist_penalty = $distance_penalty->{$distance_200}->[1];
		}
		my $scalings     = sprintf("%.6f",$contact*$activity*$dist_penalty);
		$scalings        = "2.000000" if $scalings > 2;
		my $rescaled     = sprintf("%.6f",$scalings*$H3K27ac);
		push @Enhancers,"$chrom\t$symbol\t$promoter\t$enhancer\t$distance\t$activity\t$contact\t$dist_penalty\t$scalings\t$rescaled\t$H3K27ac\n";
		next;
	    }
	    elsif($e_end < $p_start)
	    {
		my $distance = $p_start - $e_end;
		my $HiC_distance = get_HiC_distance($distance);
		my $contact      = $HiC_baseline->{$HiC_distance};
		$contact         = $Contacts{"$promoter<>$enhancer"} if defined $Contacts{"$promoter<>$enhancer"} and $Contacts{"$promoter<>$enhancer"} >= $contact;

		my $distance_200 = get_200_distance($distance);
		my $dist_penalty = "NA";

		if($strand eq "+")
		{
		    $dist_penalty = $distance_penalty->{$distance_200}->[0];
		}
		else
		{
		    $dist_penalty = $distance_penalty->{$distance_200}->[1];
		}
		my $scalings     = sprintf("%.6f",$contact*$activity*$dist_penalty);
		$scalings        = "2.000000" if $scalings > 2;
		my $rescaled     = sprintf("%.6f",$scalings*$H3K27ac);
		push @Enhancers,"$chrom\t$symbol\t$promoter\t$enhancer\t$distance\t$activity\t$contact\t$dist_penalty\t$scalings\t$rescaled\t$H3K27ac\n";
		next;
	    }
	    die "$p_start,$p_end,$symbol,$strand\n$e_start,$e_end,$activity,$H3K27ac\n";
	}
	next if $#Enhancers == -1;
	print FILE_OUT join "",@Enhancers;
	print FILE_OUT "\n";
    }
    close FILE_IN;
    print "Done $tissue.$chrom\n";
}
system "mv Tissue.scalings.data.gz $data_folder/data/$tissue.scalings.data.gz";
print "Done all\n";

sub load_baseline
{
    my %ret_h = ();
    open FILE_IN,"$data_folder/files/HiC.baseline";
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

sub load_distance_penalty
{
    my %ret_h = ();
    open FILE_IN,"$data_folder/files/Upstream.penalty";
    my $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[0]} = [$Data[1],$Data[1]];
    }
    close FILE_IN;

    open FILE_IN,"$data_folder/files/Downstream.penalty";
    $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[0]}->[1] = $Data[1];
    }
    close FILE_IN;
    return \%ret_h;
}

sub get_HiC_distance
{
    my $value = $_[0];
    return 1000000 if $value >= 997500;
    return 0       if $value <=   2500;

    my $times = int(($value+2500)/5000);
    return $times*5000;
}

sub get_200_distance
{
    my $value = $_[0];
    return 1000000 if $value >= 999900;

    my $times = int(($value+100)/200);
    return $times*200;
}

sub load_gene_data
{
    my %ret_h = ();
    open FILE_IN,"$data_folder/files/Gencode.v40.promoters.data";
    my $string = <FILE_IN>;
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[0]} = [0,0,$Data[5],$Data[2]];
    }
    close FILE_IN;
    open FILE_IN,"$data_folder/files/$regions";
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	next unless $Data[3] =~ /^ENSG/;
	next unless $ret_h{$Data[3]};
	$ret_h{$Data[3]}->[0] = $Data[1];
	$ret_h{$Data[3]}->[1] = $Data[2];
    }
    close FILE_IN;
    return \%ret_h;
}
