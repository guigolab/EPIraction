#!/usr/bin/env perl

use strict;
use warnings;

my $data_folder  = "!{data_folder}";
my $samples_file = "!{samples_index}";
my $tissue       = "!{tissue}";
my $chromosomes  = "!{chromosomes}";
my $ctcf_impact  = 0.5;
srand(127);

die "wrong folder $data_folder\n" unless -d $data_folder;

my ($samples,undef)    = load_data("$data_folder/files/$samples_file");
my ($expression,undef) = load_data("$data_folder/data/$tissue.expression.data");
my $prefix = "$data_folder/temp/samples/$tissue";

foreach my $chromosome (split ", ", $chromosomes)
{
    my %Activity = ();
    if(-f "$prefix.$chromosome.data.done")
    {
	print "Already $prefix.$chromosome.data.done\n";
    }

    open FILE_OUT,">Regions.$chromosome.bed";
    open FILE_IN,"unpigz -p3 -c $data_folder/data/$tissue.regions.data.gz |";
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
	$Activity{$Data[3]} = ["$Data[8]","$Data[9]","$Data[10]","$Data[11]"];
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
	    if($region =~/^ENSG/ and $expression->{$region}->{"expression"} >= 1)
	    {
######----------For promoters of expressed protein-coding or lncRNA genes I set minimum H3K27ac and Open since I did it calculating tissue-consensus activity ------
		my $Open     = $Activity{$region}->[0];
		my $Cofactor = $Activity{$region}->[1];
		my $CTCF     = $Activity{$region}->[2];
		$signal      = 5 if $signal < 5;
		$Open        = 5 if $Open < 5;
		my $report   = sprintf("%.8f",sqrt($Open*($signal + $Cofactor)));
		$Signal{$region} .= "\t$report";
	    }
	    elsif($Activity{$region}->[3] == 0)
	    {
######----------For non-active regions I set activity to zero------
		my $report = "0.00000000";
		$Signal{$region} .= "\t$report";
		next;
	    }
	    else
	    {
######----------For the rest I calculate Activity conventually ------
		my $Open     = $Activity{$region}->[0];
		my $Cofactor = $Activity{$region}->[1];
		my $CTCF     = $Activity{$region}->[2];
		$Open        = 1 if $Open < 1;
		my $report   = sprintf("%.8f",sqrt($Open*($signal + $Cofactor + $ctcf_impact*$CTCF)));
		$Signal{$region} .= "\t$report";
	    }
	}
	close FILE_IN;
	push @List,$accession;
    }
    print "Load $chromosome samples\n";

    open FILE_OUT,"| gzip > $tissue.$chromosome.data.gz";
    print FILE_OUT "region\tConsensus\tActivity.mean\tShift\t",(join "\t", @List),"\n";

    open FILE_IN,"Regions.$chromosome.bed";
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data   = split "\t", $string;
	$Signal{$Data[3]} =~s/^\t//;
	my @Values = split "\t", $Signal{$Data[3]};
	my ($mean,$count) = (0,0);
	map { $mean+=$_;$count++ } split "\t", $Signal{$Data[3]};
	$mean = sprintf("%.8f",$mean/$count);
	my $Consensus = $Activity{$Data[3]}->[3];
	my $difference = sprintf("%.8f",$Consensus - $mean);
	print FILE_OUT "$Data[3]\t$Consensus\t$mean\t$difference\t",$Signal{$Data[3]},"\n";
    }
    close FILE_OUT;
    system "mv $tissue.$chromosome.data.gz $prefix.$chromosome.data.gz";
    system "touch $prefix.$chromosome.data.done";
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

