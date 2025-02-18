#!/usr/bin/perl

use strict;

my $data_folder = "!{data_folder}";
my $tissues     = "!{tissues}";
my $hic_samples = "!{hic_samples}";
my $chromosomes = "!{chromosomes}";

die "wrong folder $data_folder\n" unless -d $data_folder;
$tissues     =~ s/\[|\]//g;
$hic_samples =~ s/\[|\]//g;
$chromosomes =~ s/\[|\]//g;
my @Chroms = sort split ", ", $chromosomes;

open FILE_OUT,">First_QC.log";
print "==============================";
print FILE_OUT "==============================";
foreach my $tissue (sort split ", ", $tissues)
{
    unless(-f "$data_folder/data/$tissue.regions.data.gz")
    {
	print "Absent $data_folder/data/$tissue.regions.data.gz";
	exit(127);
    }
    open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.regions.data.gz |";
    my $string = <FILE_IN>;
    my ($promoters,$enhancers) = (0,0);
    while($string = <FILE_IN>)
    {
	$string =~ /\tENSG\d+\t/ ? $promoters++ : $enhancers++;
    }
    close FILE_IN;
    print "\ntissue    = $tissue\npromoters = $promoters\nenhancers = $enhancers\n";
    print FILE_OUT "\ntissue    = $tissue\npromoters = $promoters\nenhancers = $enhancers\n";
}
print "\n==============================";
print FILE_OUT "\n==============================";

foreach my $chrom (@Chroms)
{
    unless(-f "$data_folder/data/zzz.Pairs.$chrom.gz")
    {
	print "Absent $data_folder/data/zzz.Pairs.$chrom.gz";
	exit(127);
    }
    open FILE_IN,"zcat $data_folder/data/zzz.Pairs.$chrom.gz |";
    my $string = <FILE_IN>;
    my ($promoters,$pairs) = (0,0);
    while($string = <FILE_IN>)
    {
	my @Data = split "\t", $string;
	$promoters++;
	$pairs+=$Data[1];
    }
    close FILE_IN;
    print "\nchrom     = $chrom\npromoters = $promoters\npairs     = $pairs\n";
    print FILE_OUT "\nchrom     = $chrom\npromoters = $promoters\npairs     = $pairs\n";
}
print "\n==============================";
print FILE_OUT "\n==============================";

foreach my $hic_sample (split ", ", $hic_samples)
{
    foreach my $chrom (@Chroms)
    {
	my $prefix = "$hic_sample.$chrom";
	unless(-f "$data_folder/temp/contacts/$prefix.overlap.gz")
	{
	    print "Absent $data_folder/temp/contacts/$prefix.overlap.gz";
	    exit(127);
	}
    }
    print "\nHiC: $hic_sample is good";
    print FILE_OUT "\nHiC: $hic_sample is good";
}
print "\n\n\nAll good\n";
print FILE_OUT "\n\nAll good\n";
close FILE_OUT;
print "Done First_qc\n";
