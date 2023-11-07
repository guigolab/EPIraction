#!/usr/bin/env perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $samples_file = "!{samples_index}";
my $tissue       = "!{tissue}";
my $chromosomes  = "!{chromosomes}";

die "wrong folder $data_folder\n" unless -d $data_folder;
srand(127);

my ($samples,undef) = load_data("$data_folder/files/$samples_file");
my $prefix = "$data_folder/temp/signal/$tissue";

foreach my $chromosome (split ", ", $chromosomes)
{
    my %Activity = ();
    open FILE_OUT,">Regions.$chromosome.bed";
    open FILE_IN,"$data_folder/data/$tissue.regions.data";
    my $string = <FILE_IN>;
    chomp $string;
    my @Header = split "\t", $string;
    my ($elevated,$total) = (0,0);
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	next unless $Data[0] eq $chromosome;
	$total++;
	print FILE_OUT "$Data[0]\t$Data[1]\t$Data[2]\t$Data[3]\n";
	$Data[7] = 1 if $Data[7] < 1;
	$Activity{$Data[3]} = ["$Data[7]","$Data[13]","$Data[6]"];
    }
    close FILE_IN;
    close FILE_OUT;
    print "Load activies and wrote bed file for $chromosome $tissue\n";

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
	    my ($region,undef,undef,undef,undef,undef,undef,$signal) = split "\t", $string;
	    $signal = 90 if $signal > 90;
	    if($Activity{$region}->[2] eq "enhancer")
	    {
		if($Activity{$region}->[1] == 0)
		{
		    $Signal{$region} .= "\t0.000";
		}
		else
		{
		    my $rescale = ($signal/$Activity{$region}->[0])**0.5;
		    my $report  = $rescale*$Activity{$region}->[1];
		    $report     = 110 if $report > 110;
		    $report     = sprintf("%.6f",$report);
		    $Signal{$region} .= "\t$report";
		}
	    }
	    else
	    {
		$signal     = 1 if $signal < 1;
		my $rescale = ($signal/$Activity{$region}->[0])**0.5;
		my $report  = $rescale*$Activity{$region}->[1];
		$report     = 110 if $report > 110;
		$report     = sprintf("%.6f",$report);
		$Signal{$region} .= "\t$report";
	    }
	}
	close FILE_IN;
	push @List,$accession;
    }
    print "Load $chromosome samples\n";

    open FILE_OUT,"| gzip > $tissue.$chromosome.data.gz";
    print FILE_OUT "region\tActivity.norm\tActivity.mean\tActivity.error\t",(join "\t", @List),"\n";

    open FILE_IN,"Regions.$chromosome.bed";
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$Signal{$Data[3]} =~s/^\t//;
	my ($sum,$count) = (0,0);
	map { $sum+=$_;$count++ } split "\t", $Signal{$Data[3]};
	my $Activity_norm  = $Activity{$Data[3]}->[1];
	$Activity_norm     = 0.1 if $Activity_norm < 0.1;
	my $Activity_mean  = sprintf("%.6f",$sum/$count);
	my $Activity_error = sprintf("%.6f",abs($Activity_norm - $Activity_mean)/$Activity_norm);
	$Activity_error = sprintf("%.6f",abs($Activity_norm - $Activity_mean)) if $Activity_norm < 1;
	
	print FILE_OUT "$Data[3]\t$Activity_norm\t$Activity_mean\t$Activity_error\t",$Signal{$Data[3]},"\n";
    }
    close FILE_OUT;
    system "mv $tissue.$chromosome.data.gz $prefix.$chromosome.data.gz";
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

