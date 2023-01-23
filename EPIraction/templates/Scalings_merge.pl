#!/usr/bin/perl -w
use strict;

my $data_folder = "!{data_folder}";
my $chromosome  = "!{chrom}";
my $tissues     = "!{tissues}";
my $regions     = "!{regions}";
die "wrong folder $data_folder\n" unless -d $data_folder;

if(-f "$data_folder/data/xxx.Scalings.$chromosome.data.gz")
{
    print "Already $data_folder/data/xxx.Scalings.$chromosome.data.gz\n";
    exit;
}

$tissues =~ s/\[|\]//g;
my @Tissues = sort split ", ", $tissues;
my $positions   = load_regions("$data_folder/files/$regions");

my %Scalings = ();
foreach my $tissue (@Tissues)
{
    open FILE_IN,"zcat $data_folder/data/$tissue.scalings.data.gz |";
    my $string = <FILE_IN>;
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	next if $string eq "";
	my ($chrom,$symbol,$promoter,$enhancer,$distance,$activity,$contact,$d_penalty,$scaling,$rescaled,$original) = split "\t", $string;
	next unless $chrom eq $chromosome;
	$Scalings{"$tissue<>$promoter<>$enhancer"} = sprintf("%.8f",$scaling);
    }
    close FILE_IN;
}

open FILE_OUT,"| gzip > Scalings.data.gz";
print FILE_OUT "promoter\tenhancer\tdistance\t",(join "\t", @Tissues),"\n";

open FILE_IN,"zcat $data_folder/data/zzz.Pairs.$chromosome.gz |";
my $string = <FILE_IN>;
while($string = <FILE_IN>)
{
    chomp $string;
    my ($promoter,undef,$list) = split "\t", $string;
    my ($p_start,$p_end) = @{$positions->{$promoter}};
    foreach my $enhancer (split ";", $list)
    {
	my ($e_start,$e_end) = @{$positions->{$enhancer}};
	my @Sorted       = sort { $a <=> $b } ($e_start,$e_end,$p_start,$p_end);
	my $distance     = $Sorted[2] - $Sorted[1];
	my @Report = ();
	foreach my $tissue (@Tissues)
	{
	    if(defined $Scalings{"$tissue<>$promoter<>$enhancer"})
	    {
		push @Report,$Scalings{"$tissue<>$promoter<>$enhancer"};
	    }
	    else
	    {
		push @Report,"0.000000";
	    }
	}
	print FILE_OUT "$promoter\t$enhancer\t$distance\t",(join "\t", @Report),"\n";
    }
    print FILE_OUT "\n";
}
close FILE_IN;
close FILE_OUT;
system "mv Scalings.data.gz $data_folder/data/xxx.Scalings.$chromosome.data.gz";
print "Done $chromosome\n";

sub load_regions
{
    my %ret_h = ();
    open FILE_IN,"$_[0]";
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[3]} = [$Data[1],$Data[2]];
    }
    close FILE_IN;
    return\%ret_h;
}
